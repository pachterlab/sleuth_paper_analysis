
suppressPackageStartupMessages({
library('cowplot')
library('data.table')
library('dplyr')
library('mamabear')
library('parallel')
})
options(mc.cores = n_cpu)

base_dir <- '../results/final_figures'
default_extension <- '.pdf'

################################################################################
# independent isoform simulation
################################################################################

sim_name <- 'isoform_3_3_20_1_1'

each_filter_benchmark <- readRDS(paste0('../results/', sim_name,
  '/isoform_benchmarks_filter.rds'))

each_filter_benchmark[[1]]$color_mapping <- method_colors

# because some of the true things will inevitably get filtered, remove the ones we think
# will be removed from the truth
suppressMessages(current_fdr <- get_fdr(each_filter_benchmark, sim_filter = TRUE)$pvals)
tmp <- fdr_efdr_power_plot(current_fdr, start = 100, jump = 100, rank_fdr = NULL)

p <- tmp +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_cowplot(25)
p

filename <- file.path(base_dir, paste0('each_filter_nozoom_', sim_name,
  default_extension))
save_plot(filename, p, base_aspect_ratio = 1.6, base_height = 15)

################################################################################
# fold changes learned from reference (gfr)
################################################################################

sim_name <- 'gfr_3_3_20_42_2'

###
# isoform level
###

each_filter_benchmark <- readRDS(paste0('../results/', sim_name,
  '/isoform_benchmarks_filter.rds'))

# because some of the true things will inevitably get filtered, remove the ones we think
# will be removed from the truth
suppressMessages(current_fdr <- get_fdr(each_filter_benchmark, sim_filter = TRUE)$pvals)

tmp <- fdr_efdr_power_plot(current_fdr, start = 500, jump = 500, rank_fdr = NULL,
  de_bench = each_filter_benchmark[[1]])

p <- tmp +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_cowplot(25)
p

filename <- file.path(base_dir, paste0('isoform.each_filter_nozoom_', sim_name,
  default_extension))
save_plot(filename, p, base_aspect_ratio = 1.6, base_height = 15)

#
tmp <- fdr_efdr_power_plot(current_fdr, start = 100, jump = 100, rank_fdr = 0.10,
  de_bench = each_filter_benchmark[[1]])

p <- tmp +
  coord_cartesian(xlim = c(-0.01, 0.25), ylim = c(-0.01, 0.075), expand = FALSE) +
  theme_cowplot(25)
p

###
# gene level
###

each_filter_benchmark <- readRDS(paste0('../results/', sim_name,
  '/gene_benchmarks_filter.rds'))

suppressMessages(current_fdr <- get_fdr(each_filter_benchmark,
  sim_filter = FALSE)$pvals)

tmp <- fdr_efdr_power_plot(current_fdr, start = 500, jump = 500, rank_fdr = 0.10,
  de_bench = each_filter_benchmark[[1]])

p <- tmp +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_cowplot(25)
p

filename <- file.path(base_dir, paste0('gene.each_filter_nozoom_', sim_name,
  default_extension))
save_plot(filename, p, base_aspect_ratio = 1.6, base_height = 15)

#
tmp <- fdr_efdr_power_plot(current_fdr, start = 100, jump = 100, rank_fdr = 0.10,
  de_bench = each_filter_benchmark[[1]])

p <- tmp +
  coord_cartesian(xlim = c(-0.01, 0.25), ylim = c(-0.01, 0.25), expand = FALSE) +
  theme_cowplot(25)
p
