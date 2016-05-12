args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop('Usage: Rscript bench_isoform.R N_CPU SIM_NAME')
}

# temporary for debugging
# n_cpu <- 20
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

# put the isoform information for EBSeq in global variables
NG_LIST <- GetNg(transcript_gene_mapping$target_id,
  transcript_gene_mapping$ens_gene)

###
# new refactoring
###

all_results <- list()

all_results$limmaVoom <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_isoform_results_intersect(sim_name, i, 'limmaVoom',
      limma_filter_and_run)
  }, mc.cores = n_cpu)

all_results$DESeq <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_isoform_results_intersect(sim_name, i, 'DESeq', DESeq_filter_and_run)
  }, mc.cores = n_cpu)

all_results$DESeq2 <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_isoform_results_intersect(sim_name, i, 'DESeq2',
      DESeq2_filter_and_run_intersect)
  }, mc.cores = n_cpu)

all_results$edgeR <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_isoform_results_intersect(sim_name, i, 'edgeR',
      edgeR_filter_and_run)
  }, mc.cores = n_cpu)

all_results$EBSeq <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_isoform_results_intersect(sim_name, i, 'EBSeq',
      EBSeq_isoform_filter_and_run)
  }, mc.cores = n_cpu)

all_results$Cuffdiff2 <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_isoform_results_intersect(sim_name, i, 'Cuffdiff2',
      cuffdiff_filter_and_run)
  }, mc.cores = n_cpu)

all_benchmarks <- lapply(all_results,
  function(result) {
    mclapply(seq_along(result),
      function(i) {
        x <- result[[i]]
        sim_info <- get_de_info(sim_name, i, transcript_gene_mapping)
        new_de_benchmark(x, names(x), sim_info$de_info)
      }, mc.cores = n_cpu)
  })

saveRDS(all_benchmarks, file = paste0('../results/', sim_name,
  '/isoform_benchmarks.rds'))


#
# saveRDS(all_nde_plot, file = '../sims/new_all_nde_plot.rds')
#
# all_nde_plot <- readRDS('../sims/new_all_nde_plot.rds')
#
#
# saveRDS(all_nde_plot, file = '../sims/temporary_all_nde_plot.rds')
#
# all_nde_plot <- readRDS('../sims/temporary_all_nde_plot.rds')
#
# ###
# # TODO: remove this... for debugging
# ###
#
# cat('Sample: ', i, '\n')
#
# sim_info <- get_de_info(sim_name, 1, transcript_gene_mapping)
#
# tmp <- load_isoform_results_intersect(sim_name, 1, 'limmaVoom', limma_filter_and_run,
#   include_empirical_bayes = TRUE)
#
# tmp_bench <- new_de_benchmark(tmp, names(tmp), sim_info$de_info)
#
# fdr_nde_plot(tmp_bench) +
#   xlim(0, 3000) +
#   ylim(0, 0.10) +
#   theme_cowplot(25) +
#   theme(legend.position = c(0.15, 0.80))
