document("~/dev/sleuth")
install("~/dev/sleuth")

library("cowplot")
library("ggplot2")
library("ggvis")
library("nlme")

base_dir <- "../data/de_4_4"

de_truth <- fread(file.path(base_dir, "de.txt"), header = TRUE, data.table = FALSE)

de_counts <- fread(file.path(base_dir, "counts.txt"), header = TRUE,
  data.table = FALSE)

de_truth <- de_truth %>%
  mutate(size = est_counts / 3) %>%
  mutate(dispersion = 3 / est_counts)

kal_fnames <- Sys.glob(file.path(base_dir, "*", "kallisto", "expression.h5"))

kal_obj <- lapply(kal_fnames, read_kallisto_h5, read_bootstrap = TRUE)

ids <- kal_fnames %>%
  sub(paste0(base_dir, "/"), "", .) %>%
  sub(file.path("", "kallisto", "expression.h5"), "", .) %>%
  paste0("sample_", .)

id_to_condition <- data.frame(sample = paste0("sample_", 1:8),
  condition = c(rep("a", 4), rep("b", 4)))

################################################################################
s_o <- new_sleuth(kal_obj, id_to_condition, ~ condition, norm_bootstraps = FALSE)
s_o$obs_norm <- s_o$obs_raw

s_o <- prep_bs_plots(s_o)

de_truth %>%
  filter(!is_de, tpm > 1 & tpm < 2)

plot_transcript_variability(s_o, 31, "condition")

hi <- lmm_design(s_o)
bs_est_counts <- dcast_bootstrap(s_o, "est_counts")
fit_lmm(hi, log(expression + 0.5) ~ condition, list(sample = ~ 1),
  bs_est_counts['ENST00000217233',])

pass_filt <- s_o$obs_norm %>%
  group_by(target_id) %>%
  summarise( count_filt = sum(est_counts > 5) > 0.8 * n())

save(s_o, all_lmm, pass_filt, file = "../data/lmm.RData")

debugonce(lmm_by_row)

system.time(all_lmm <- lmm_by_row(s_o, log(expression + 0.5) ~ condition,
    list(sample = ~ 1), pass_filt))


system.time(all_lmm_ml <- lmm_by_row(s_o, log(expression + 0.5) ~ condition,
    list(sample = ~ 1), pass_filt, "ML"))

################################################################################

load("../data/lmm.RData", verbose = TRUE)

names(all_lmm) <- s_o$kal[[1]]$bootstrap[[1]]$target_id
filt_lmm <- all_lmm %>%
  Filter(function(x) is(x, "lme"), .)

names(all_lmm_ml) <- s_o$kal[[1]]$bootstrap[[1]]$target_id
filt_lmm_ml <- all_lmm_ml %>%
  Filter(function(x) is(x, "lme"), .)

summary_lmm <- lapply(filt_lmm, summary)

pvals <- summary_lmm %>%
  seq_along() %>%
  lapply(
    function(i)
    {
      id <- names(summary_lmm)[i]
      s <- summary_lmm[[i]]
      data.frame(target_id = id, b = s$tTable["conditionb", "Value"],
        pval = s$tTable["conditionb", "p-value"],
        var_b = s$varFix[1,1],
        sigma = summary_lmm[[i]]$sigma,
        t = s$tTable["conditionb", "t-value"])
    }) %>%
      rbind_all()

summary_lmm_ml <- filt_lmm_ml %>%
  seq_along() %>%
  lapply(
    function(i)
    {
      id <- names(filt_lmm_ml)[i]
      s <- summary(filt_lmm_ml[[i]])
      #s_summary <- summary(s)
      data.frame(target_id = id, b = s$coefficients$fixed[2],
        pval = s$tTable["conditionb", "p-value"],
        sigma = s$sigma,
        var_b = s$varFix[1,1],
        t = s$tTable["conditionb", "t-value"])
    })
summary_lmm_ml_df <- rbind_all(summary_lmm_ml)

pvals_truth <- de_truth %>%
  select(target_id, fold_change, valid_trans, est_counts, is_de) %>%
  inner_join(pvals, by = c("target_id"))

################################################################################
# mv relationship
################################################################################

s_o <- sleuth_summarize_bootstrap(s_o)

cv_summary <- s_o$bootstrap_summary %>%
  group_by(target_id) %>%
  summarise(cv_est_counts = mean(bs_cv_est_counts, na.rm = TRUE))

LWR <- 0.5
UPR <- 20

s_o$obs_norm %>%
  filter(est_counts > 0) %>%
  ggplot(aes(log(est_counts))) +
    geom_histogram(aes(y = ..density..)) +
    xlim(0, 10) +
    facet_wrap(~ sample)


mv <- s_o %>%
  null_mean_var(transform = function(x) log(x), 0.0) %>%
  filter(counts_mean > 0) %>%
  inner_join(data.table(cv_summary), by = c("target_id")) %>%
  arrange(counts_mean) %>%
  mutate(grp = cut(ecdf(cv_est_counts)(cv_est_counts), 10))


smooth <- with(mv, lowess(counts_mean, sqrt(sqrt(counts_var)), f = 0.5))
#smooth <- with(mv, lowess(counts_mean, sqrt(sqrt(counts_var))))
mv <- mv %>%
  mutate(smooth_var = smooth$y)

debugonce(get_quantile)

alpha <- 0.10
eps <- 0.25
tmp <- s_o %>%
  null_mean_var(transform = function(x) log(x), 0.0) %>%
  filter(counts_mean > 0) %>%
  # compute the ECDF of the counts based on mean counts
  mutate(exp_ecdf = ecdf(counts_mean)(counts_mean)) %>%
  # break the expression into 100 different bins
  mutate(grp = cut(exp_ecdf, 100)) %>%
  data.frame() %>%
  # remove the lowly expressed and highly expressed stuff
  filter(0 + alpha <= exp_ecdf & exp_ecdf <= 1 - alpha) %>%
  group_by(grp) %>%
  # compute the sliding window using the "center" of the data
  do(get_quantile(., "counts_var", eps, 1 - eps))

lfit <- with(filter(tmp, cdf),
  lowess(counts_mean, sqrt(sqrt(counts_var)))) %>%
  as.data.frame()

ggplot(tmp, aes(counts_mean, sqrt(sqrt(counts_var)))) +
  geom_point(aes(colour = cdf), alpha = 0.2) +
  geom_line(aes(x, y), data = lfit, colour = "blue") +
  xlim(0, 10) +
  ylim(0, 4)


  arrange(counts_mean) %>%
  mutate(smooth_x = smooth$x, smooth_var = smooth$y) %>%
  ggplot(aes(counts_mean, sqrt(sqrt(counts_var))))  +
    geom_point(aes(colour = obs_grp), alpha = 0.2) +
    geom_line(aes(counts_mean, smooth_var), colour = "blue") +
    ylim(0, 2) +
    xlim(0, 15)

debugonce(null_mean_var)

mv <- s_o %>%
  null_mean_var(transform = function(x) log(x))

mv %>%
  as.data.frame() %>%
  filter(counts_mean > 0.1) %>%
  ggvis(x = ~counts_mean, y = ~sqrt(sqrt(counts_var))) %>%
    layer_points(fill := "black", opacity := 0.2) %>%
    add_tooltip(function(dat){
      print(dat)
      print(class(dat))
      }, "hover") %>%
    scale_numeric('x', domain = c(LWR, UPR), nice = FALSE)


s_o %>%
  null_mean_var(transform = function(x) log(x), 0.0) %>%
  filter(counts_mean > 0) %>%
  inner_join(data.table(cv_summary), by = c("target_id")) %>%
  mutate(grp = cut(ecdf(cv_est_counts)(cv_est_counts), 10)) %>%
  arrange(counts_mean) %>%
  mutate(smooth_x = smooth$x, smooth_var = smooth$y) %>%
  ggplot(aes(counts_mean, cv_est_counts))  +
    geom_point(aes(colour = grp), alpha = 0.2)

################################################################################
# shrink sigma
################################################################################


mv <- s_o %>%
  null_mean_var(transform = function(x) log(x))

mv_pvals <- inner_join(mv, data.table(pvals), by = "target_id")

mv_pvals <- with(mv_pvals,
  lowess(counts_mean, sigma)) %>%
  as.data.frame() %>%
  bind_cols(mv_pvals)

plt_lmm <- ggplot(mv_pvals, aes(counts_mean, sigma)) +
  geom_point(alpha = 0.2) +
  geom_line(aes(x, y), colour = "blue")

plt_sigma <- ggplot(mv_pvals, aes(counts_mean, sqrt(counts_var))) +
  geom_point(alpha = 0.2)

YMAX <- 1
plt <- plot_grid(plt_lmm + ylim(0, YMAX), plt_sigma + ylim(0, YMAX))
dev.print(pdf, "../img/sigma_lmm_vs_raw.pdf")

mv_pvals <- mv_pvals %>%
  mutate(sigma_shrink = pmax(sigma, y)) %>%
  mutate(sigma_shrink_2 = pmin(sigma, y))


summary_lmm_shrink <- Map(
  function(lme_obj, sigma)
  {
    summary_lme(lme_obj, sigma)
  }, filt_lmm[mv_pvals$target_id], mv_pvals$sigma_shrink_2)

pvals_shrink <- summary_lmm_shrink %>%
  seq_along() %>%
  lapply(
    function(i)
    {
      id <- names(summary_lmm_shrink)[i]
      s <- summary_lmm_shrink[[i]]
      data.frame(target_id = id, b = s$tTable["conditionb", "Value"],
        pval = s$tTable["conditionb", "p-value"],
        sigma = summary_lmm_shrink[[i]]$sigma)
    }) %>%
      rbind_all()
pvals_shrink <- pvals_shrink %>%
  mutate(qval = p.adjust(pval, method = "BH"))

debugonce(nlme:::summary.lme)

debugonce(nlme:::MEestimate)
summary(filt_lmm[[1]])

################################################################################
# here, we want to get a working copy of lmm so we can figure out what it's
# doing
################################################################################

dcast_res <- dcast_bootstrap(s_o, "est_counts")
tmp_trans <- "ENST00000361739"
tmp_design <- model.matrix(~ condition, lmm_design(s_o))

debugonce(lmm)
ah <- sleuth::lmm(tmp_design, log(dcast_res[tmp_trans,] + 0.5), 8)
sqrt( ah$sigma_sq )
filt_lmm[[tmp_trans]]$sigma
filt_lmm[[tmp_trans]]

# let's see how the closed form lmm correlates with nlme::lme
system.time(hp_filt_lmm <- lapply(names(filt_lmm),
  function(target_id)
  {
    lmm(tmp_design, log(dcast_res[target_id,] + 0.5), 8)
  }))

hp_filt_lmm_df <- data.frame(target_id = names(filt_lmm),
  hp_sigma = sapply(hp_filt_lmm, function(x) sqrt(x$sigma)),
  hp_b = sapply(hp_filt_lmm, function(x) x$coef[2]),
  hp_var_b1 = sapply(hp_filt_lmm, function(x) x$cov_b[1,1]),
  hp_t = sapply(hp_filt_lmm, function(x) x$t_value[2]),
  stringsAsFactors = FALSE)

lme_lmm_sigma <- inner_join(hp_filt_lmm_df, select(pvals, target_id, var_b, sigma, t),
  by = "target_id")

lme_ml_lmm_sigma <- inner_join(hp_filt_lmm_df, summary_lmm_ml_df,
  by = "target_id")

ggplot(lme_lmm_sigma, aes(hp_sigma, sigma)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  geom_smooth(method = "lm", colour = "blue")
with(lme_lmm_sigma, cor(hp_sigma, sigma, use = "complete.obs"))

ggsave("../img/lme_lmm_sigma.png")

p_reml_b <- ggplot(lme_lmm_sigma, aes(hp_var_b1, var_b)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  geom_smooth(method = "lm", colour = "blue") +
  xlim(0, 1.5) +
  ylim(0, 1.5) +
  ylab("var_b using REML")

p_ml_b <- ggplot(lme_ml_lmm_sigma, aes(hp_var_b1, var_b)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  geom_smooth(method = "lm", colour = "blue") +
  xlim(0, 1.5) +
  ylim(0, 1.5) +
  ylab("var_b using ML")

plt <- plot_grid(p_reml_b, p_ml_b)
plt
dev.print(pdf, "../img/lme_ml_lmm_vbeta.pdf")
dev.off()

p_reml_t <- ggplot(lme_lmm_sigma, aes(abs(hp_t), abs(t))) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  geom_smooth(method = "lm", colour = "blue") +
  # xlim(0, 1.5) +
  # ylim(0, 1.5) +
  ylab("t using REML")

p_ml_t <- ggplot(lme_ml_lmm_sigma, aes(abs(hp_t), abs(t))) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  geom_smooth(method = "lm", colour = "blue") +
  # xlim(0, 1.5) +
  # ylim(0, 1.5) +
  ylab("t using ML")

plt <- plot_grid(p_reml_t, p_ml_t)
plt
dev.print(pdf, "../img/lme_lmm_t.pdf")
dev.off()

save_plot("../img/lme_ml_lmm_vbeta.png", plt)

################################################################################
# fit my own lmm
################################################################################

debugonce(compute_lmms)
system.time(s_o <- compute_lmms(s_o, function(x) log(x + 0.5), pass_filt))
s_o$sigma$bs_mean <- s_o$bs_means

debugonce(naive_shrink_lmm)
s_o <- naive_shrink_lmm(s_o)

plot_mean_var(s_o) +
  geom_line(aes(bs_mean, naive_locfit), colour = "red")

s_o$sigma %>%
  head()
s_o$adjustment <- s_o$adjustment[1,1]

all_ts <- compute_t(s_o, "naive_shrink_sigma")

all_ts_raw <- compute_t(s_o, "raw_sigma")

naive_shrink_pvals <- all_ts$p_vals %>%
  data.frame() %>%
  select(pval = conditionb) %>%
  mutate(target_id = rownames(all_ts$p_vals)) %>%
  mutate(qval = p.adjust(pval, method = "BH"))

raw_pvals <- all_ts_raw$p_vals %>%
  data.frame() %>%
  select(pval = conditionb) %>%
  mutate(target_id = rownames(all_ts_raw$p_vals)) %>%
  mutate(qval = p.adjust(pval, method = "BH"))

colnames(all_ts$p_vals)[2] <- "shrink_pval"

################################################################################
# looking at the p-values output from my method vs the R lme method
shrink_pvals <- all_ts$p_vals %>%
  as.data.frame() %>%
  mutate(target_id = rownames(all_ts$p_vals)) %>%
  select(-`(Intercept)`)

tmp <- inner_join(shrink_pvals, summary_lmm_ml_df, by = "target_id")


ggplot(tmp, aes(pval, shrink_pval)) +
  geom_point(alpha = 0.04) +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  geom_smooth(method = "lm")

# put t-values and p-values in the same table for 'conditionb'
raw_pvals <- all_ts_raw$p_vals %>%
  as.data.frame() %>%
  mutate(target_id = rownames(all_ts_raw$p_vals)) %>%
  select(-`(Intercept)`)
raw_pvals <- rename(raw_pvals, raw_pval = conditionb)
raw_pvals <- all_ts_raw$t_stats %>%
  as.data.frame() %>%
  mutate(target_id = rownames(all_ts_raw$t_stats)) %>%
  select(-`(Intercept)`, raw_t_stat = conditionb) %>%
  inner_join(raw_pvals, by = "target_id")

tmp2 <- inner_join(raw_pvals, summary_lmm_ml_df, by = "target_id")

ggplot(tmp2, aes(pval, raw_pval)) +
  geom_point(alpha = 0.04) +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  geom_smooth(method = "lm")

ggplot(tmp2, aes(abs(t), abs(raw_t_stat))) +
  geom_point(alpha = 0.04) +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  geom_smooth(method = "lm") +
  xlim(0, 100) +
  ylim(0, 100)

##########
# there seems to be an issue with the p-values and the t-statistics. a lot of
# stuff that has a smaller p-value using my method vs the lme method. let's
# investigate that...
##########

# first, let's compare raw sigmas again...
raw_sigma_df <- s_o$sigma[c(1)] %>%
  mutate(target_id = rownames(s_o$sigma)) %>%
  inner_join(summary_lmm_ml_df, by = "target_id")

# correlation is very high... looks good.
with(raw_sigma_df, cor(raw_sigma, sigma))

# looks almost perfect. issue is NOT in sigma. must be in the way var(\beta) is computed
ggplot(raw_sigma_df, aes(sigma, raw_sigma)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = "lm")

fixed_sd <- s_o$fixed_sd %>%
  data.frame() %>%
  select(hp_sd_b = conditionb) %>%
  mutate(target_id = rownames(s_o$fixed_sd))

fixed_sd <- summary_lmm_ml_df %>%
  mutate(sd_b = sqrt(var_b * s_o$adjustment)) %>%
  inner_join(fixed_sd, by = "target_id")

# this plot shows that my calculation of the SD is much greater than the
# R version (lme).  this issue doesn't explain why there is a strange variation
# in the calculation of p-values (e.g. p-values that are less in our version)
ggplot(fixed_sd, aes(sd_b, hp_sd_b)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  geom_smooth(method = "lm", colour = "blue") +
  xlim(0, 2) +
  ylim(0, 2)

all_ts_raw <- compute_t(s_o, "raw_sigma")

# let's compare the p-values using compute_t versus doing them from my closed
# form output
t_values <- sapply(s_o$lmms, function(x) x$t_value[2]) %>%
  data_frame(t_value = ., target_id = names(s_o$lmms))

compute_t_values <- data_frame(compute_t_value = all_ts_raw$t_stats[,2],
  target_id = rownames(all_ts_raw$t_stats))

# after finding bug, they now seem consistent -- woohoo!
ggplot(inner_join(t_values, compute_t_values, by = "target_id"),
  aes(t_value, compute_t_value)) +
    geom_point(alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1)

################################################################################

