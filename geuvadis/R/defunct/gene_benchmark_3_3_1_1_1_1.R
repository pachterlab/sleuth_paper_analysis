library("cowplot")
library("sleuth")
library("mamabear")


sim_name <- 'gene_3_3_1_1_1_1'
sims_dir <- file.path('../sims', sim_name)
base_dir <- file.path('../results', sim_name)

# de_info contains the 'true' fold change as well as which transcripts are DE
de_info <- read.table(gzfile(file.path(sims_dir, "de_info.tsv.gz")),
  header = TRUE, stringsAsFactors = FALSE)

de_genes <- select(de_info, ens_gene, is_de) %>%
  group_by(ens_gene) %>%
  summarize(is_de = any(is_de))
de_genes <- rename(de_genes, target_id = ens_gene)
de_genes <- mutate(de_genes, log_fc = NA)

kal_dirs <- file.path(base_dir, "exp_1", 1:6, "kallisto")

s2c <- data.frame(sample = paste0("sample_", 1:6),
  condition = rep(c("A", "B"), each = 3), stringsAsFactors = FALSE)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host="www.ensembl.org")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

min_reads <- 7
# so <- sleuth_prep(s2c, ~ condition, max_bootstrap = 30,
#   target_mapping = t2g, min_reads = min_reads, min_prop = 0.5)
so <- sleuth_prep(s2c, ~ condition, norm_fun_counts = all_ones, max_bootstrap = 30,
  target_mapping = t2g, min_reads = min_reads, min_prop = 0.5)
# so <- sleuth_prep(s2c, ~ condition, norm_fun_counts = all_ones, min_prop = 0.5, min_reads = 10)

so <- sleuth_fit(so)

models(so)

so <- sleuth_wt(so, 'conditionB')
models(so)

# sleuth_live(so)

sres <- sleuth_results(so, 'conditionB') %>%
  dplyr::select(target_id, pval, qval)

so <- sleuth_fit(so, ~1, "reduced")

so <- sleuth_lrt(so, "reduced", "full")
s_ratio <- sleuth_results(so, 'reduced:full', test_type = 'lrt')

# save(s_ratio, t2g, de_info, file = 'gene_data.RData')

gene_table <- sleuth_gene_table(so, 'reduced:full')
# s_ratio_genes <- left_join(s_ratio, t2g, by = "target_id")
s_ratio_genes <- gene_table
gene_lift <- s_ratio_genes %>%
  group_by(ens_gene) %>%
  do({
    min_index <- which.min(.$qval)
    result <- .[min_index, ]
    result <- select(result, -target_id)
    result <- rename(result, target_id = ens_gene)
    })

gene_table_wt <- sleuth_gene_table(so, 'conditionB', test_type = 'wt')
gene_lift_wt <- gene_table_wt %>%
  group_by(ens_gene) %>%
  do({
    min_index <- which.min(.$qval)
    result <- .[min_index, ]
    result <- select(result, -target_id)
    result <- rename(result, target_id = ens_gene)
    })


################################################################################
# counts analysis

counts <- read.table(file.path(base_dir, "exp_1", 'counts.tsv'), header = TRUE,
  stringsAsFactors = FALSE, row.names = 1)
counts <- as.matrix(counts)
colnames(counts) <- sub('X', '', colnames(counts))

genes_filter <- apply(counts, 1, basic_filter, min_reads = min_reads, min_prop = 0.5)
counts_filtered <- counts[genes_filter, ]

library("DESeq2")
debugonce(runDESeq2)

tmp <- runDESeq2(counts_filtered)

################################################################################
# TEMPORARY: use the true counts for the gene counting programs
#
# sim_info <- readRDS(file.path(sims_dir, 'sim.rds'))
# counts_true <- sim_info$counts[[1]]
# counts_true <- as.data.frame(counts_true, stringsAsFactors = FALSE)
# counts_true <- mutate(counts_true, target_id = rownames(counts_true))
# counts_true <- melt(counts_true, 'target_id', variable.name = 'sample')
# counts_true <- left_join(counts_true, t2g, by = 'target_id')
# counts_true <- group_by(counts_true, ens_gene, sample) %>%
#   summarize(value = sum(value))
# counts_matrix <- dcast(counts_true, ens_gene ~ sample)
# counts_matrix <- counts_matrix[complete.cases(counts_matrix), ]
# rownames(counts_matrix) <- counts_matrix$ens_gene
# counts_matrix$ens_gene <- NULL
#
# genes_filter <- apply(counts_matrix, 1, basic_filter, min_reads = min_reads,
#   min_prop = 0.5)
# counts_filtered <- counts_matrix[genes_filter, ]

################################################################################
# DESeq
################################################################################

dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_filtered,
  colData = select(s2c, condition),
  design = ~ condition)

# force to not estimate size factors
DESeq2::sizeFactors(dds) <- rep(1, 6)

#dds <- DESeq2::DESeq(dds, fitType='local')
dds <- DESeq2::DESeq(dds)
deseq_res <- DESeq2::results(dds)

deseq_res <- as.data.frame(deseq_res) %>%
  mutate(target_id = rownames(deseq_res)) %>%
  dplyr::rename(qval = padj, pval = pvalue)

################################################################################
# limma
################################################################################

dge <- edgeR::DGEList(counts = counts_filtered)

# not normalizing because something odd is happening
dge <- edgeR::calcNormFactors(dge)

limma_design <- model.matrix(~ condition, s2c)
v <- limma::voom(dge, limma_design, plot = TRUE)

limma_fit <- limma::lmFit(v, limma_design)
limma_fit <- limma::eBayes(limma_fit)
limma_res <- limma::topTreat(limma_fit, coef = ncol(limma_design), number = nrow(counts_filtered))
limma_res <- limma_res %>%
  mutate(target_id = rownames(limma_res)) %>%
  dplyr::select(target_id, pval = P.Value, qval = adj.P.Val)

################################################################################
# edgeR
################################################################################

dge <- edgeR::DGEList(counts = counts_filtered, group = s2c$condition)
dge <- edgeR::calcNormFactors(dge)
y <- edgeR::estimateCommonDisp(dge)
y <- edgeR::estimateTagwiseDisp(y)
et <- edgeR::exactTest(y)
edgeR_res <- edgeR::topTags(et, n = nrow(counts_filtered)) %>%
  as.data.frame() %>%
  mutate(., target_id = rownames(.)) %>%
  dplyr::select(target_id, pval = PValue, qval = FDR)

################################################################################
# benchmark
################################################################################
de_bench <- new_de_benchmark(
  list(
    gene_lift,
    gene_lift_wt,
    deseq_res,
    limma_res,
    edgeR_res
    ),
  c(
    "sleuth lrt",
    "sleuth wt",
    "DESeq2",
    "limma",
    "edgeR (tagwise)"
    ), de_genes)

fdr_nde_plot(de_bench, FALSE) +
  # xlim(0, 2000) +
  xlim(0, 3000) +
  ylim(0, 0.10) +
  theme(legend.position = c(0.1, 0.80))+
  xlab('number of genes called DE') +
  ylab('FDR')
