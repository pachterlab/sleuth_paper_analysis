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

################################################################################

load("../data/lmm.RData", verbose = TRUE)

names(all_lmm) <- s_o$kal[[1]]$bootstrap[[1]]$target_id
filt_lmm <- all_lmm %>%
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
        sigma = summary_lmm[[i]]$sigma)
    }) %>%
      rbind_all()

pvals <- pvals %>%
  mutate(qval = p.adjust(pval, method = "BH"))

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



