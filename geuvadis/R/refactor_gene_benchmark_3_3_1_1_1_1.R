source("benchmark_methods.R")

# sim_name <- 'gene_3_3_1_1_1_1'
sim_name <- 'gene_diff_fc_3_3_1_1_1_1'
source("gene_common.R")

true_counts <- read.table(paste0('../results/', sim_name, '/exp_1/counts.tsv'))

kal_dirs <- file.path(base_dir, "exp_1", 1:6, "kallisto")

sample_to_condition <- data.frame(sample = paste0("sample_", 1:6),
  condition = rep(c("A", "B"), each = 3), stringsAsFactors = FALSE)
sample_to_condition <- dplyr::mutate(sample_to_condition, path = kal_dirs)

min_reads <- 5

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

sum(genes_filter)

###
# run the methods
###

# sr <- run_sleuth(sample_to_condition, min_reads = min_reads)
# sgr <- run_sleuth(sample_to_condition, min_reads = 5, lift_genes = TRUE)

so <- sleuth_prep(sample_to_condition, ~condition, min_reads = 5,
  gene_mode = 'ens_gene', target_mapping = transcript_gene_mapping)

bs_summary <- so$bs_summary
gene_names_filter <- names(genes_filter[which(genes_filter)])
gene_names_filter <- intersect(rownames(bs_summary$obs_counts),
  gene_names_filter)

so$bs_summary$obs_counts <- bs_summary$obs_counts[gene_names_filter, ]
so$bs_summary$sigma_q_sq <- bs_summary$sigma_q_sq[gene_names_filter]

so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')

so <- sleuth_wt(so, which_beta = 'conditionB')
so <- sleuth_lrt(so, 'reduced', 'full')

lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt', show_all = FALSE)[,
  c('target_id', 'pval', 'qval')]
wt <- sleuth_results(so, 'conditionB', show_all = FALSE)[,
  c('target_id', 'pval', 'qval')]
lrt <- lrt[complete.cases(lrt), ]
wt <- wt[complete.cases(wt), ]
sgr <- list(sleuth.lrt = lrt, sleuth.wt = wt)


gene_methods <- list(
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


library('cowplot')
all_results <- c(sgr, gene_results)

de_bench <- new_de_benchmark(all_results, names(all_results), de_genes)

fdr_nde_plot(de_bench, FALSE) +
  xlim(0, 3000) +
  ylim(0, 0.10) +
  theme(legend.position = c(0.1, 0.80)) +
  xlab('number of genes called DE') +
  ylab('FDR') +
  theme_bw(25)


###

tmp <- semi_join(de_info, data.table(ens_gene = gene_names_filter), by = "ens_gene")

pl_bench <- new_de_benchmark(all_results, names(all_results), pl_genes)

fdr_nde_plot(pl_bench, FALSE) +
  xlim(0, 25) +
  ylim(0, 0.35) +
  theme(legend.position = c(0.1, 0.80)) +
  xlab('number of genes called DE') +
  ylab('FDR') +
  theme_bw(25)
