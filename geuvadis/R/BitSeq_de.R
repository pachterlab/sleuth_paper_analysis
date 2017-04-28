# args <- commandArgs(trailingOnly = TRUE)
#
# if (length(args) != 3) {
#   stop('Usage: Rscript BitSeq_de.R SIM_NAME TRANSCRIPTOME_FA OUTPUT_PREFIX')
# }



sim_name <- 'gfr_3_3_20_42_2'
sim_number <- 1

# sim_name <- args[1]

source('benchmark_methods.R')
source('gene_common.R')

base_dir <- '../results/final_figures'
default_extension <- '.pdf'

sim_info <- parse_simulation(sim_name)

n <- sim_info$a + sim_info$b

expression_file <- list()
expression_file$a <- file.path('..', 'sims', sim_name,
  paste0('exp_', sim_number), 1:(sim_info$a))
expression_file$b <- file.path('..', 'sims', sim_name,
  paste0('exp_', sim_number), (sim_info$a + 1):(sim_info$a + sim_info$b))

expression_file <- lapply(expression_file,
  function(file_names) {
    file.path(file_names, 'BitSeq_cli.rpkm')
  })

library('BitSeq')

bs_timing <- system.time(res <- getDE(expression_file, samples = TRUE, seed = 42))
saveRDS(bs_timing, file = file.path(destination, 'bs_timing.rds'))

destination <- file.path('..', 'results', sim_name, sim_number)
suppressWarnings(dir.create(destination))

saveRDS(res, file = file.path(destination, 'BitSeq_de.rds'))

###
# eventually export this into the analysis
###
sim_name <- args[1]


sim_name <- 'gfr_3_3_20_42_2'
sim_number <- 1


source('benchmark_methods.R')
source('gene_common.R')

base_dir <- '../results/final_figures'
default_extension <- '.pdf'

sim_info <- parse_simulation(sim_name)

n <- sim_info$a + sim_info$b


bitseq <- readRDS(file.path('..', 'results', sim_name, sim_number,
  'BitSeq_de.rds'))

# length differential expression as recommended in the Bioconductor vignette
bitseq$pplr <- bitseq$pplr[ order(abs(0.5 - bitseq$pplr$pplr), decreasing=TRUE ), ]

bs_1 <- bitseq$pplr
bs_1$target_id <- rownames(bs_1)
bs_1 <- as.data.frame(bs_1)

# load sleuth result
each_filter_benchmark <- readRDS(paste0('../results/', sim_name,
  '/isoform_benchmarks_filter_lfc.rds'))

sr_1 <- each_filter_benchmark[[1]]$original_data$sleuth

# make a dummy column that looks like a pval
bs_1 <- dplyr::mutate(bs_1, rank = 1:length(target_id))
bs_1 <- dplyr::mutate(bs_1, pval = rank / length(target_id))
bs_1 <- dplyr::mutate(bs_1, qval = pval)

bs_filtered_1 <- dplyr::semi_join(bs_1, sr_1, by = 'target_id')
bs_filtered_1 <- dplyr::select(bs_filtered_1, target_id, pval, qval)

bs_1 <- dplyr::select(bs_1, target_id, pval, qval)

library('mamabear')

filter_benchmark_bitseq <- new_de_benchmark(
  c(each_filter_benchmark[[1]]$original_data, list(bs_filtered_1, bs_1)),
  c(each_filter_benchmark[[1]]$labels, 'BitSeq filtered', 'BitSeq'),
  oracle = each_filter_benchmark[[1]]$oracle, join_mode = 'union')
# turn into a list so that we can call the list version of `get_fdr`
filter_benchmark_bitseq <- rename_benchmark(filter_benchmark_bitseq,
  c('Cuffdiff2', 'limmaVoom'),
  c('Cuffdiff 2', 'voom'), join_mode = 'union')
filter_benchmark_bitseq <- list(filter_benchmark_bitseq)

fdr_bitseq <- suppressMessages(get_fdr(filter_benchmark_bitseq,
  sim_filter = TRUE)$pvals)

#
fdr_bitseq <- dplyr::filter(fdr_bitseq, !grepl('LFC', method))
method_colors_bs <- c(method_colors, 'BitSeq filtered' = '#4393c3', 'BitSeq' = '#fddbc7')

# produce a zoomed in version of the plot
tmp <- fdr_efdr_power_plot(fdr_bitseq, start = 100, jump = 100, rank_fdr = 0.10,
  method_colors = method_colors_bs, fdr_level_position = -0.005,
  ignore_estimated_fdr = c('BitSeq filtered', 'BitSeq'))

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

filename <- file.path(base_dir, paste0('BitSeq_isoform.each_filter_', sim_name,
  default_extension))
save_plot(filename, p, base_aspect_ratio = 1.6, base_height = 15)

# produce the zoomed out version of the plot
tmp <- fdr_efdr_power_plot(fdr_bitseq, start = 500, jump = 500, rank_fdr = 0.10,
  method_colors = method_colors_bs, fdr_level_position = -0.02,
  ignore_estimated_fdr = c('BitSeq filtered', 'BitSeq'))

p <- tmp +
  coord_cartesian(xlim = c(-0.05, 1), ylim = c(-0.05, 1), expand = FALSE) +
  theme_hp()
p <- p +
  geom_polygon(aes(x, y), alpha = 0.20,
    data = data.frame(
    x = c(0, 0, current_limits$x[2], current_limits$x[2]),
    y = c(0, current_limits$y[2], current_limits$y[2], 0)))
p

filename <- file.path(base_dir, paste0('BitSeq_isoform.each_filter_nozoom_',
  sim_name, default_extension))
save_plot(filename, p, base_aspect_ratio = 1.6, base_height = 15)
