# args <- commandArgs(trailingOnly = TRUE)

# if (length(args) != 2) {
#   stop('Usage: Rscript bench_isoform.R N_CPU SIM_NAME')
# }

# temporary for debugging
# n_cpu <- 20
# sim_name <- 'isoform_3_3_20_1_1'
sim_name <- 'gfr_3_3_20_42_2'

# n_cpu <- args[1]
# sim_name <- args[2]
base_dir <- '../results/final_figures'
default_extension <- '.pdf'

source("benchmark_methods.R")
source("gene_common.R")

library("parallel")
library("mamabear")

transcript_gene_mapping <- get_human_gene_names()

# keep this as a placeholder and replace the path later
sim <- parse_simulation(sim_name)

sample_info <- lapply(1:20,
  function(i) {
    n <- sim$a + sim$b
    kal_dirs <- file.path('..', 'sims', sim_name, paste0("exp_", i),
      1:n, "kallisto")
    get_sample_to_condition(sim$a, sim$b, kal_dirs)
  })

# this stuff was run by Nick and deposited here
rsem_base <- '/home/hjp/sleuth_paper_analysis/geuvadis/sims/gfr_3_3_20_42_2/exp_1'

rsem_h5 <- file.path(rsem_base, 1:6, 'rsem.h5')

sample_info[[1]]$path <- rsem_h5
sample_info <- sample_info[[1]]

res <- run_sleuth(sample_info)
res$so <- NULL

saveRDS(res, file = '../results/gfr_3_3_20_42_2/1/rsem_sleuth.rds')


each_filter_benchmark <- readRDS(paste0('../results/', sim_name,
  '/isoform_benchmarks_filter_lfc.rds'))

sr_1 <- each_filter_benchmark[[1]]$original_data$sleuth

library('mamabear')

rsem <- readRDS('../results/gfr_3_3_20_42_2/1/rsem_sleuth.rds')

# add BitSeq to put it on the same scale as the other plot

bitseq <- readRDS(file.path('..', 'results', sim_name, 1,
  'BitSeq_de.rds'))

# length differential expression as recommended in the Bioconductor vignette
bitseq$pplr <- bitseq$pplr[ order(abs(0.5 - bitseq$pplr$pplr), decreasing=TRUE ), ]

bs_1 <- bitseq$pplr
bs_1$target_id <- rownames(bs_1)
bs_1 <- as.data.frame(bs_1)
# make a dummy column that looks like a pval
bs_1 <- dplyr::mutate(bs_1, rank = 1:length(target_id))
bs_1 <- dplyr::mutate(bs_1, pval = rank / length(target_id))
bs_1 <- dplyr::mutate(bs_1, qval = pval)

bs_filtered_1 <- dplyr::semi_join(bs_1, sr_1, by = 'target_id')
bs_filtered_1 <- dplyr::select(bs_filtered_1, target_id, pval, qval)

bs_1 <- dplyr::select(bs_1, target_id, pval, qval)

filter_benchmark_rsem <- new_de_benchmark(
  list(sr_1, rsem$sleuth.lrt, bs_1, bs_filtered_1),
  # list(sr_1, rsem$sleuth.lrt),
  # c(each_filter_benchmark[[1]]$original_data, list(res$sleuth.lrt)),
  c('sleuth', 'sleuth (RSEM)', 'BitSeq', 'BitSeq (filtered)'),
  # c(each_filter_benchmark[[1]]$labels, 'sleuth (RSEM)'),
  oracle = each_filter_benchmark[[1]]$oracle, join_mode = 'union')

# turn into a list so that we can call the list version of `get_fdr`
# filter_benchmark_rsem <- rename_benchmark(filter_benchmark_rsem,
#   c('Cuffdiff2', 'limmaVoom'),
#   c('Cuffdiff 2', 'voom'), join_mode = 'union')
# filter_benchmark_rsem <- list(filter_benchmark_rsem)

fdr_rsem <- suppressMessages(get_fdr(list(filter_benchmark_rsem),
  sim_filter = TRUE)$pvals)

# remove BitSeq because we don't care about it here (point of comparison is
# kallisto-sleuth versus rsem-sleuth)
fdr_rsem <- dplyr::filter(fdr_rsem, !grepl('BitSeq', method))

# produce a zoomed in version of the plot
tmp <- fdr_efdr_power_plot(fdr_rsem, start = 100, jump = 100, rank_fdr = 0.10,
  # method_colors = c(method_colors$sleuth, 'sleuth (RSEM)' = 'black'),
  method_colors = c(method_colors, 'LFC' = 'gray', 'GLFC' = 'lightgray',
  'sleuth (RSEM)' = 'black'),
  fdr_level_position = -0.005)

simulation_mode <- 'reference'
current_limits <- switch(simulation_mode,
  independent = list(x = c(-0.01, 0.25), y = c(-0.01, 0.28)),
  common = list(x = c(-0.01, 0.25), y = c(-0.01, 0.20)),
  reference = list(x = c(-0.01, 0.25), y = c(-0.01, 0.075))
  )

library('cowplot')
theme_hp <- function() {
  theme_cowplot(25) +
    theme(legend.key.size = unit(2, "lines"))
}

p <- tmp + theme_hp()
p <- p + coord_cartesian(xlim = current_limits$x, ylim = current_limits$y,
  expand = FALSE)
p

filename <- file.path(base_dir, paste0('rsem_isoform.each_filter_', sim_name,
  default_extension))
save_plot(filename, p, base_aspect_ratio = 1.6, base_height = 15)

# produce the zoomed out version of the plot
tmp <- fdr_efdr_power_plot(fdr_rsem, start = 500, jump = 500, rank_fdr = 0.10,
  method_colors = c(method_colors, 'LFC' = 'gray', 'GLFC' = 'lightgray',
  'sleuth (RSEM)' = 'black'),
  fdr_level_position = -0.02)

p <- tmp +
  coord_cartesian(xlim = c(-0.05, 1), ylim = c(-0.05, 1), expand = FALSE) +
  theme_hp()
p <- p +
  geom_polygon(aes(x, y), alpha = 0.20,
    data = data.frame(
    x = c(0, 0, current_limits$x[2], current_limits$x[2]),
    y = c(0, current_limits$y[2], current_limits$y[2], 0)))
p

filename <- file.path(base_dir, paste0('rsem_isoform.each_filter_nozoom_',
  sim_name, default_extension))
save_plot(filename, p, base_aspect_ratio = 1.6, base_height = 15)
