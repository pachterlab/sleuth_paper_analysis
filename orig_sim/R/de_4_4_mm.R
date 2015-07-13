document("~/dev/sleuth")
install("~/dev/sleuth")

library("sleuth")
library("data.table")
library("ggplot2")

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
# fit a fixed effect for each sample
################################################################################

fen <- fixed_null_glm_all_bs(s_o, 30)
fen_summary <- fen %>%
  Filter(function(x) is(x, "glm"), .) %>%
  Map(summary, .)
fen_dispersion <- fen_summary %>%
  lapply(function(x) x$dispersion) %>%
  unlist()
fen_df <- data.frame(target_id = names(fen_dispersion), disp_fen = fen_dispersion)

# only intercept, so totally pooled
fen_oi <- fixed_null_glm_all_bs_no_grouping(s_o, 30) %>%
  Filter(function(x) is(x, "glm"), .) %>%
  Map(summary, .)

fen_oi_disp <- fen_oi %>%
  sapply(function(x) x$dispersion)
fen_oi_df <- data.frame(target_id = names(fen_oi_disp),
  disp_fen_oi = fen_oi_disp)

ggplot(inner_join(fen_df, fen_oi_df, by = "target_id"),
  aes(disp_fen, disp_fen_oi)) +
    geom_point(alpha = 0.2) +
    geom_abline(slope = 1, intercept = 0) +
    scale_x_log10(limits = c(0.0001, 2)) +
    scale_y_log10(limits = c(0.0001, 2))
