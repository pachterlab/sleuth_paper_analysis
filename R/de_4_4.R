library("sleuth")
library("data.table")

base_dir <- "../data/de_4_4"

de_truth <- fread(file.path(base_dir, "de.txt"), header = TRUE, data.table = FALSE)

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

