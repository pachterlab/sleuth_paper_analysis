source("benchmark_methods.R")

sim_name <- 'gene_3_3_1_1_1_1'
source("gene_common.R")

kal_dirs <- file.path(base_dir, "exp_1", 1:6, "kallisto")

sample_to_condition <- data.frame(sample = paste0("sample_", 1:6),
  condition = rep(c("A", "B"), each = 3), stringsAsFactors = FALSE)
sample_to_condition <- dplyr::mutate(sample_to_condition, path = kal_dirs)

min_reads <- 7

###
# read in the count data
###

counts <- read.table(file.path(base_dir, "exp_1", 'counts.tsv'), header = TRUE,
  stringsAsFactors = FALSE, row.names = 1)
counts <- as.matrix(counts)
colnames(counts) <- sub('X', '', colnames(counts))

genes_filter <- apply(counts, 1, basic_filter, min_reads = min_reads, min_prop = 0.5)
counts_filtered <- counts[genes_filter, ]
gene_cds <- make_count_data_set(counts_filtered, sample_to_condition)


###
# run the methods
###

sr <- run_sleuth(sample_to_condition, min_reads = min_reads)
sgr <- run_sleuth(sample_to_condition, min_reads = min_reads, lift_genes = TRUE)

gene_methods <- list(
  DESeq2 = runDESeq2,
  edgeR = runEdgeR,
  edgerRobust = runEdgeRRobust,
  EBSeq = runEBSeq
  )
gene_results <- lapply(gene_methods,
  function(f) {
    f(gene_cds, FALSE)
  })

gene_results <- c(sgr, gene_results)

de_bench <- new_de_benchmark(gene_results, names(gene_results), de_genes)

fdr_nde_plot(de_bench, FALSE) +
  # xlim(0, 2000) +
  xlim(0, 3000) +
  ylim(0, 0.10) +
  theme(legend.position = c(0.1, 0.80))+
  xlab('number of genes called DE') +
  ylab('FDR')
