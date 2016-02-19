devtools::document("~/dev/sleuth")

devtools::install("~/dev/sleuth")

library("cowplot")
library("sleuth")

base_dir <- "../results/3_3_1_1"

# # de_info contains the 'true' fold change as well as which transcripts are DE
# de_info <- read.table(gzfile(file.path(base_dir, "de_info.tsv.gz")),
#   header = TRUE, stringsAsFactors = FALSE)

kal_dirs <- file.path(base_dir, "sample_1", 1:6, "kallisto")

s2c <- data.frame(sample = paste0("sample_", 1:6),
  condition = rep(c("A", "B"), each = 3), stringsAsFactors = FALSE)

sobj <- sleuth_prep(kal_dirs, s2c$sample, s2c, ~ condition)

sobj <- sleuth_fit(sobj)

models(sobj)

sobj <- sleuth_test(sobj, 'conditionB')
models(sobj)

sleuth_results(sobj, 'conditionB') %>% head

plot_mean_var(sobj)

plot_pca(sobj, pc_x = 3, color_by = 'condition', text_labels = TRUE)
plot_pca(sobj, pc_x = 3, color_by = 'condition')

ggplot(tmp, aes_string)

sleuth_interact(sobj)

debugonce(plot_scatter)

plot_scatter(sobj, offset = 0, xlim = c(0, 12), ylim = c(0, 12))

ah <- 'log2'
as.name(ah)
debugonce(plot_scatter)
plot_scatter(sobj, offset = 0, xlim = c(0, 12), ylim = c(0, 12), xtrans = eval(as.name(ah)))

plot_scatter(sobj, offset = 0, xlim = c(0, 12), ylim = c(0, 12), trans = ah)

plot_scatter(sobj, xtrans = NULL, ytrans = NULL, offset = 0, xlim = c(0, 12), ylim = c(0, 12))

hi <- spread_abundance_by(sobj$obs_norm, 'est_counts')
hi <- as.data.frame(hi)
hi$target_id <- rownames(hi)

ggplot(hi, aes(sample_1 + 1, sample_2 + 1)) +
  geom_point() +
  coord_trans(xtrans = 'log', ytrans = 'log')


