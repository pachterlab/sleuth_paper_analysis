# library('devtools')
# install('~/sleuth_revision')

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop('Usage: Rscript bench_isoform.R N_CPU SIM_NAME')
}

# temporary for debugging
n_cpu <- 20
# sim_name <- 'isoform_3_3_20_1_1'
sim_name <- 'gfr_3_3_20_42_2'

n_cpu <- args[1]
sim_name <- args[2]


source("benchmark_methods.R")
source("gene_common.R")

library("parallel")
library("mamabear")
options(mc.cores = n_cpu)

transcript_gene_mapping <- get_human_gene_names()

sim <- parse_simulation(sim_name)

sample_info <- lapply(1:20,
  function(i) {
    n <- sim$a + sim$b
    kal_dirs <- file.path('..', 'sims', sim_name, paste0("exp_", i),
      1:n, "kallisto")
    get_sample_to_condition(sim$a, sim$b, kal_dirs)
  })

zero_variance_filter <- list()

# zero_variance_filter$sleuth_zero_variance <- lapply(1,
zero_variance_filter$sleuth_zero_variance <- mclapply(1:sim$n,
  function(i) {
    si <- sample_info[[i]]
    res <- run_sleuth(si, zero_technical_variance = TRUE)
    res$so <- NULL
    res
  })

zero_variance_filter$sleuth_zero_variance <- lapply(zero_variance_filter$sleuth_zero_variance,
  function(x) {
    x$sleuth.lrt
  })

each_filter_benchmark <- readRDS(paste0('../results/', sim_name,
  '/isoform_benchmarks_filter_lfc.rds'))

each_filter <- lapply(each_filter_benchmark[[1]]$labels,
  function(label) {
    lapply(each_filter_benchmark, function(bench) {
      bench$original_data[[label]]
    })
  })
names(each_filter) <- each_filter_benchmark[[1]]$labels

each_filter_zero_variance <- c(each_filter, zero_variance_filter)
each_filter_benchmark_zero_variance <- lapply(1:sim$n,
  function(i) {
    current <- lapply(each_filter_zero_variance, '[[', i)
    current_names <- names(each_filter_zero_variance)
    sim_info <- get_de_info(sim_name, i, transcript_gene_mapping)
    new_de_benchmark(current, current_names, sim_info$de_info,
      join_mode = 'union')
  })

saveRDS(each_filter_benchmark_zero_variance, file = paste0('../results/', sim_name,
  '/isoform_benchmarks_filter_zero_variance.rds'))
