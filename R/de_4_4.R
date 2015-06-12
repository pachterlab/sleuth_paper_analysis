library("sleuth")
library("data.table")

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

s_o <- new_sleuth(kal_obj, id_to_condition, ~ condition)

################################################################################
# looking at dispersion
################################################################################

de_truth %>%
  filter(valid_trans) %>%
  ggplot(aes(est_counts, dispersion)) +
    geom_point(alpha = 0.2) +
    xlim(0, 250)

################################################################################
# exploring GLMs
################################################################################

obs_norm <- spread_abundance_by(s_o$obs_norm, "est_counts")

de_truth %>%
  filter(is_de) %>%
  arrange(desc(tpm)) %>%
  head(20)
# 'ENST00000329251' looks good

cur_trans <- 'ENST00000329251'
cur_resp <- obs_norm[cur_trans,]

cur_glm <- fit_gamma_glm(s_o$design_mat, cur_resp) %>%
  structure(class = c("glm", "lm"))

samp_bs <- sample_bootstrap(s_o, 100L)

# ~ 50s
system.time(obs_glms <- glm_by_rows(obs_norm, s_o$design_matrix))
is_valid_obs_glms <- unlist(lapply(obs_glms, is, "glm"))
valid_obs_glms <- obs_glms[ is_valid_obs_glms ]
valid_obs_summary <- lapply(valid_obs_glms, summary)

# ~35s per it
system.time(samp_glm <- lapply(seq_along(samp_bs)[1:30], function(x)
    {
      if (x %% 10 == 0) {
        cat("it:",x,"\n")
      }
      glm_coef_by_rows(samp_bs[[x]], s_o$design_matrix)
    }))

samp_glm_merge <- lapply(samp_glm,
  function(samp)
  {
    samp <- data.frame(samp) %>%
      mutate(target_id = rownames(samp))
    samp[complete.cases(samp),]
  }) %>%
    rbind_all()
setnames(samp_glm_merge, c("intercept", "conditionb", "target_id"))

samp_glm_merge %>%
  group_by(target_id) %>%
  mutate()

################################################################################

# currently only takes one coefficient for the p-value calculation
glm_summary_to_df <- function(glm_summary) {
  coef_data <- as.data.frame((glm_summary)$coef)[2,]
  data.frame(
    conditionb = coef_data[1,"Estimate"],
    se_conditionb = coef_data[1, "Std. Error"],
    dispersion = glm_summary$dispersion,
    pval = coef_data[1, "Pr(>|t|)"]
    )
}

valid_obs_dat <- lapply(valid_obs_summary, glm_summary_to_df) %>%
  rbind_all() %>%
  mutate(target_id = names(valid_obs_summary))

valid_obs_dat <- valid_obs_dat %>%
  mutate(qval = p.adjust(pval, method = "BH"))

# compute SE on bootstraps
samp_glm_se <- samp_glm_merge %>%
  group_by(target_id) %>%
 summarise(
    bs_mean_conditionb = mean(conditionb),
    bs_se_conditionb = sd(conditionb)
    )

# merge the two tables
valid_obs <- inner_join(
  data.table(valid_obs_dat),
  data.table(samp_glm_se),
  by = c("target_id"))

valid_obs <- valid_obs %>%
  mutate(qval = p.adjust(pval, method = "BH"))

valid_obs <- de_truth %>%
  select(-size, -dispersion, fc_a, fc_b) %>%
  inner_join(valid_obs, by = c("target_id")) %>%
  mutate(est_fold_change = exp(conditionb))

pe <- function(est, true) {
  (est - true) / true
}

save.image(file = "de_4.RData")

load("de_4.Rdata", verbose = TRUE)

################################################################################
# fit the null GLMs
################################################################################


system.time(null_obs_glms <- glm_by_rows(obs_norm,
    model.matrix( ~ 1, s_o$sample_to_condition)))

valid_null_obs_glms <- Filter(function(x) is(x, "glm"), null_obs_glms)
valid_null_summary <- Map(summary, valid_null_obs_glms)
valid_null_mean <- apply(obs_norm[names(valid_null_summary),], 1, mean)
valid_null_phi <- Map(function(x) x$dispersion, valid_null_summary)

valid_null_df <- data.frame(
  target_id = names(valid_null_phi),
  null_mean = valid_null_mean,
  null_phi = unlist(valid_null_phi),
  stringsAsFactors = FALSE
  )

valid_null_df <- valid_null_df %>%
  inner_join(select(valid_obs, target_id, is_de), by = c("target_id"))

ggplot(valid_null_df, aes(null_mean, null_phi)) +
  geom_point(aes(colour = is_de), alpha = 0.1) +
  scale_color_manual(values = c("black", "red")) +
  xlim(0, 1000) +
  ylim(0, 1)

s_o <- fit_null_glms(s_o)
lfp <- locfit_phi(s_o)
lfr <- locfit::locfit.raw(log(lfp[[1]]), lfp[[2]], family = "gamma", maxk = 1000)

lf_df <- data.frame(target_id = names(s_o$null_glm_summary),
  null_mean = lfp[[1]], null_phi = lfp[[2]])

lf_df <- lf_df %>%
  mutate(smooth_phi = predict(lfr, log(null_mean)))

ggplot(lf_df, aes(null_mean, null_phi)) +
  geom_point(alpha = 0.2) +
  geom_line(aes(null_mean, smooth_phi), colour = "red") +
  xlim(0, 1000) +
  ylim(0, 1)

document("~/dev/sleuth")
install("~/dev/sleuth")

hi <- dcast_bootstrap(s_o$kal[[1]], "est_counts")

hi <- dcast_bootstrap(s_o, "est_counts", 30)

debugonce(bs_null_obs_glm_fit)
all_bs_glm <- bs_null_obs_glm_fit(s_o)
x <- apply(hi, 1, function(row) all(row > 0))
y <- which(!x)
z <- intersect(names(null_phis), names(y))
zp <- intersect(names(null_phis), names(which(x))

bs_glm <- bs_null_obs_glm_fit(s_o, 30) %>%
    Filter(function(x) is(x, "glm"), .) %>%
    Map(summary, .)

tech_phis <- bs_glm %>%
  lapply(function(x) x$dispersion) %>%
  unlist()

null_phis <- s_o$null_glm_summary %>%
  lapply(function(x) x$dispersion) %>%
  unlist()

length(null_phis)
length(tech_disp)
length(intersect(names(null_phis), names(tech_disp)))

isect_phis <- intersect(names(null_phis), names(tech_disp))
data.frame(null_phi = null_phis[isect_phis], tech_phi = tech_phis[isect_phis]) %>%
  ggplot(aes(null_phi, tech_phi)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_log10() +
  scale_y_log10()

true_dat <- data.frame(target_id = isect_phis, stringsAsFactors = FALSE) %>%
  inner_join(valid_obs, by = "target_id")


s_o <- sleuth_summarize_bootstrap(s_o)
bs_means <- spread_abundance_by(s_o$bootstrap_summary, "bs_mean_est_counts")

bs_means_null_glm <- glm_by_rows(bs_means, model.matrix( ~1, s_o$sample_to_condition )) %>%
  Filter(function(x) is(x, "glm"), .)
bs_means_null_glm <- Map(summary, bs_means_null_glm)

bs_null_phis <- bs_means_null_glm %>%
  lapply(function(x) x$dispersion) %>%
  unlist()


isect_phis <- intersect(names(bs_null_phis), names(tech_disp))
data.frame(bs_null_phi = bs_null_phis[isect_phis], tech_phi = tech_phis[isect_phis]) %>%
  ggplot(aes(bs_null_phi, tech_phi)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_log10() +
  scale_y_log10()
