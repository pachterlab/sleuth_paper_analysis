args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop('Usage: Rscript bench_isoform.R N_CPU SIM_NAME')
}

# temporary for debugging
n_cpu <- 20
sim_name <- 'isoform_3_3_20_1_1'
# sim_name <- 'gfr_3_3_20_42_2'

n_cpu <- args[1]
sim_name <- args[2]


source("benchmark_methods.R")
source("gene_common.R")

library("parallel")
library("mamabear")
options(mc.cores = n_cpu)

transcript_gene_mapping <- get_human_gene_names()

# put the isoform information for EBSeq in global variables
NG_LIST <- GetNg(transcript_gene_mapping$target_id,
  transcript_gene_mapping$ens_gene)

###
# run each method on their own filter
###

sim <- parse_simulation(sim_name)

sample_info <- lapply(1:20,
  function(i) {
    n <- sim$a + sim$b
    kal_dirs <- file.path('..', 'sims', sim_name, paste0("exp_", i),
      1:n, "kallisto")
    get_sample_to_condition(sim$a, sim$b, kal_dirs)
  })

smooth_variance_filter <- list()

smooth_variance_filter$sleuth_smooth_variance <- mclapply(1:sim$n,
  function(i) {
    si <- sample_info[[i]]
    res <- run_sleuth(si, which_bio_var = 'smooth_sigma_sq')
    res$so <- NULL
    res
  })

smooth_variance_filter$sleuth_smooth_variance <- lapply(smooth_variance_filter$sleuth_smooth_variance,
  function(x) {
    x$sleuth.lrt
  })

each_filter_benchmark <- readRDS(paste0('../results/', sim_name,
  '/isoform_benchmarks_filter_poisson_variance.rds'))

each_filter <- lapply(each_filter_benchmark[[1]]$labels,
  function(label) {
    lapply(each_filter_benchmark, function(bench) {
      bench$original_data[[label]]
    })
  })
names(each_filter) <- each_filter_benchmark[[1]]$labels

each_filter_zero_variance <- c(each_filter, smooth_variance_filter)
each_filter_benchmark_zero_variance <- lapply(1:sim$n,
  function(i) {
    current <- lapply(each_filter_zero_variance, '[[', i)
    current_names <- names(each_filter_zero_variance)
    sim_info <- get_de_info(sim_name, i, transcript_gene_mapping)
    new_de_benchmark(current, current_names, sim_info$de_info,
      join_mode = 'union')
  })

saveRDS(each_filter_benchmark_zero_variance, file = paste0('../results/', sim_name,
  '/isoform_benchmarks_filter_smooth_variance.rds'))
