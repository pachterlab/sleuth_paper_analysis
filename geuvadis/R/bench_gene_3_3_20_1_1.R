source("benchmark_methods.R")
source("gene_common.R")

library("parallel")
library("cowplot")
library("mamabear")

sim_name <- 'isoform_3_3_20_1_1'

n_cpu <- 20

n_cpu <- 2

transcript_gene_mapping <- get_human_gene_names()

# debugonce(load_gene_results_intersect)
debugonce(edgeR_filter_and_run)
tmp <- load_gene_results_intersect(sim_name, 1, 'edgeR', edgeR_filter_and_run)

de_info <- get_de_info(sim_name, 1, transcript_gene_mapping)
de_genes <- de_info$de_genes

tmp_bench <- new_de_benchmark(tmp, names(tmp), de_genes)

fdr_nde_plot(tmp_bench) +
  xlim(0, 4500) +
  ylim(0, 0.10) +
  theme_cowplot(25) +
  theme(legend.position = c(0.15, 0.80))

###
# DESeq2
###

debugonce(limma_filter_and_run)

tmp <- load_gene_results_intersect(sim_name, 1, 'limmaVoom', limma_filter_and_run,
  include_empirical_bayes = TRUE)

tmp <- load_gene_results_intersect(sim_name, 1, 'DESeq2',
  DESeq2_filter_and_run_intersect)

# debugonce(EBSeq_gene_filter_and_run)
debugonce(runEBSeq)

tmp <- load_gene_results_intersect(sim_name, 1, 'EBSeq',
  EBSeq_gene_filter_and_run)

de_info <- get_de_info(sim_name, 1, transcript_gene_mapping)
de_genes <- de_info$de_genes

tmp_bench <- new_de_benchmark(tmp, names(tmp), de_genes)

fdr_nde_plot(tmp_bench) +
  xlim(0, 3000) +
  ylim(0, 0.10) +
  theme_cowplot(25) +
  theme(legend.position = c(0.15, 0.80))

#
###
# refactor
###

n_cpu <- 1

all_results <- list()

all_results$limmaVoom <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_gene_results_intersect(sim_name, i, 'limmaVoom',
      limma_filter_and_run, include_empirical_bayes = TRUE)
      # limma_filter_and_run, include_empirical_bayes = TRUE)
  }, mc.cores = n_cpu)

all_results$DESeq2 <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_gene_results_intersect(sim_name, i, 'DESeq2',
      DESeq2_filter_and_run_intersect, include_empirical_bayes = TRUE)
      # DESeq2_filter_and_run_intersect)
  }, mc.cores = n_cpu)

all_results$Cuffdiff2 <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_gene_results_intersect(sim_name, i, 'Cuffdiff2',
      cuffdiff_filter_and_run)
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
        new_de_benchmark(x, names(x), sim_info$de_info)
      }, mc.cores = n_cpu)
  })

all_nde_plot <- mclapply(all_benchmarks,
  function(bench) {
    fdr_nde_plot(bench) +
      xlim(0, 3000) +
      ylim(0, 0.10) +
      theme_cowplot(25) +
      theme(legend.position = c(0.15, 0.80))
  })

###
# debugging
###

all_results$limmaVoom <- lapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_gene_results_intersect(sim_name, i, 'limmaVoom',
      limma_filter_and_run)
  })

all_benchmarks <- lapply(all_results,
  function(result) {
    lapply(seq_along(result),
      function(i) {
        x <- result[[i]]
        sim_info <- get_de_info(sim_name, i, transcript_gene_mapping)
        new_de_benchmark(x, names(x), sim_info$de_genes)
      })
  })

sim_info <- get_de_info(sim_name, 1, transcript_gene_mapping)
tmp <- all_benchmarks[['limmaVoom']][[1]]
fdr_nde_plot(tmp)
temporary_bench <- new_de_benchmark(tmp, names(tmp), sim_info$de_info)

all_nde_plot <- lapply(all_benchmarks,
  function(bench) {
    fdr_nde_plot(bench) +
      xlim(0, 3000) +
      ylim(0, 0.10) +
      theme_cowplot(25) +
      theme(legend.position = c(0.15, 0.80))
  })
