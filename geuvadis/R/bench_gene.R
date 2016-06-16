args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop('Usage: Rscript bench_gene.R N_CPU SIM_NAME')
}

# temporary for debugging
n_cpu <- 20
sim_name <- 'gfr_3_3_20_42_2'

n_cpu <- args[1]
sim_name <- args[2]


source("benchmark_methods.R")
source("gene_common.R")

library("parallel")
library("mamabear")
options(mc.cores = n_cpu)

transcript_gene_mapping <- get_human_gene_names()

###
# run the each method and intersect it with sleuth
###

all_results <- list()

all_results$DESeq <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_gene_results_intersect(sim_name, i, 'DESeq',
      DESeq_filter_and_run)
  }, mc.cores = n_cpu)

all_results$DESeq2 <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_gene_results_intersect(sim_name, i, 'DESeq2',
      DESeq2_filter_and_run_intersect)
  }, mc.cores = n_cpu)

all_results$limmaVoom <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_gene_results_intersect(sim_name, i, 'limmaVoom',
      limma_filter_and_run)
  }, mc.cores = n_cpu)

all_results$Cuffdiff2 <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_gene_results_intersect(sim_name, i, 'Cuffdiff2',
      cuffdiff_filter_and_run_gene)
  }, mc.cores = n_cpu)

all_results$edgeR <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_gene_results_intersect(sim_name, i, 'edgeR',
      edgeR_filter_and_run)
  }, mc.cores = n_cpu)

all_results$EBSeq <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_gene_results_intersect(sim_name, i, 'EBSeq',
      EBSeq_gene_filter_and_run)
  }, mc.cores = n_cpu)

all_benchmarks <- lapply(all_results,
  function(result) {
    mclapply(seq_along(result),
      function(i) {
        x <- result[[i]]
        sim_info <- get_de_info(sim_name, i, transcript_gene_mapping)
        new_de_benchmark(x, names(x), sim_info$de_genes)
      }, mc.cores = n_cpu)
  })

saveRDS(all_benchmarks, file = paste0('../results/', sim_name,
  '/gene_benchmarks.rds'))
