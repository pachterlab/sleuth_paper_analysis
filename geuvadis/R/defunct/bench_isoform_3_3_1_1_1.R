source("benchmark_methods.R")
source("gene_common.R")

library("cowplot")
library("mamabear")

###
# simulation specific constants
###
sim_name <- 'isoform_3_3_1_1_1'
sims_dir <- file.path('../sims', sim_name)
base_dir <- file.path('../results', sim_name)
n_a <- 3
n_b <- 3
n <- n_a + n_b

sim_info <- get_de_info(sim_name, transcript_gene_mapping)
de_info <- sim_info$de_info
de_genes <- sim_info$de_genes

kal_dirs <- file.path(base_dir, "exp_1", 1:n, "kallisto")
sample_to_condition <- get_sample_to_condition(n_a, n_b, kal_dirs)

###
# run the isoform analysis
###
sir <- run_sleuth(sample_to_condition, lift_genes = FALSE)

isoform_cds <- get_filtered_isoform_cds(sir$so, sample_to_condition)

gene_methods <- list(
  DESeq2 = runDESeq2,
  edgeR = runEdgeR,
  limmaVoom = runVoom
  # edgerRobust = runEdgeRRobust,
  # EBSeq = runEBSeq
  )
isoform_results <- lapply(gene_methods,
  function(f) {
    f(isoform_cds, FALSE)
  })

all_results <- c(Filter(is.data.frame, sir) , isoform_results)

de_bench <- new_de_benchmark(all_results, names(all_results), de_info)

fdr_nde_plot(de_bench, FALSE) +
  xlim(0, 3000) +
  ylim(0, 0.10) +
  theme(legend.position = c(0.1, 0.80)) +
  xlab('number of genes called DE') +
  ylab('FDR') +
  theme_cowplot(25)

fdr_tpr_plot(de_bench) +
  xlim(0, 0.25) +
  theme_cowplot(25)

fdr_efdr_plot(de_bench) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0, 0.25) +
  ylim(0, 0.25) +
  theme_cowplot(25)

# gene lifting
sgr <- run_sleuth(sample_to_condition, min_reads = 5, lift_genes = TRUE)

gene_counts <- get_gene_counts(sim_name)
gene_filter <- apply(gene_counts, 1, basic_filter)
gene_counts_filtered <- gene_counts[gene_filter, ]
gene_cds <- make_count_data_set(gene_counts_filtered, sample_to_condition)

gene_methods_count <- list(
  DESeq2 = runDESeq2,
  edgeR = runEdgeR,
  limmaVoom = runVoom
  # edgerRobust = runEdgeRRobust,
  # EBSeq = runEBSeq
  )
gene_results <- lapply(gene_methods,
  function(f) {
    f(gene_cds, FALSE)
  })

all_gene_results <- c(Filter(is.data.frame, sgr), gene_results)

de_gene_bench <- new_de_benchmark(all_gene_results, names(all_gene_results),
  de_genes)

fdr_nde_plot(de_gene_bench, FALSE) +
  xlim(0, 4200) +
  ylim(0, 0.10) +
  theme(legend.position = c(0.1, 0.80)) +
  xlab('number of genes called DE') +
  ylab('FDR') +
  theme_cowplot(25)

fdr_tpr_plot(de_gene_bench) +
  xlim(0, 0.25) +
  theme_cowplot(25)

fdr_efdr_plot(de_gene_bench) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0, 0.25) +
  ylim(0, 0.25) +
  theme_cowplot(25)
