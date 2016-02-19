devtools::document("~/dev/sleuth")
devtools::install("~/dev/sleuth")

library("cowplot")
library("sleuth")

sim_name <- '3_3_1_1_2'
sims_dir <- file.path('../sims', sim_name)
base_dir <- file.path('../results', sim_name)

# de_info contains the 'true' fold change as well as which transcripts are DE
de_info <- read.table(gzfile(file.path(sims_dir, "de_info.tsv.gz")),
  header = TRUE, stringsAsFactors = FALSE)

sim <- readRDS(file.path(sims_dir, 'sim.rds'))
kal_dirs <- file.path(base_dir, "exp_1", 1:6, "kallisto")

s2c <- data.frame(sample = paste0("sample_", 1:6),
  condition = rep(c("A", "B"), each = 3), stringsAsFactors = FALSE)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
so <- sleuth_prep(s2c, ~ condition)


so <- sleuth_fit(so)

models(so)

so <- sleuth_wt(so, 'conditionB')
models(so)

sleuth_interact(so)

sres <- sleuth_results(so, 'conditionB') %>%
  dplyr::select(target_id, pval, qval)

so <- sleuth_fit(so, ~1, "reduced")

so <- sleuth_lrt(so, "reduced", "full")
s_ratio <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
s_ratio <- dplyr::select(s_ratio, target_id, pval, qval)

################################################################################
# DESeq 2
#library("DESeq2")
pass_filt_names <- so$filter_df[['target_id']]

obs_raw <- sleuth:::spread_abundance_by(so$obs_raw, "est_counts")
obs_raw_filt <- obs_raw[pass_filt_names,]

dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(obs_raw_filt),
  colData = select(s2c, condition),
  design = ~ condition)

#dds <- DESeq2::DESeq(dds, fitType='local')
dds <- DESeq2::DESeq(dds)
deseq_res <- DESeq2::results(dds)

deseq_res <- as.data.frame(deseq_res) %>%
  mutate(target_id = rownames(deseq_res)) %>%
  dplyr::rename(qval = padj, pval = pvalue)

################################################################################
# limma
################################################################################

dge <- edgeR::DGEList(counts = obs_raw_filt)

# not normalizing because something odd is happening
# ddge <- edgeR::calcNormFactors(dge)

limma_design <- model.matrix(~ condition, s2c)
v <- limma::voom(dge, limma_design, plot = TRUE)

limma_fit <- limma::lmFit(v, limma_design)
limma_fit <- limma::eBayes(limma_fit)
limma_res <- limma::topTreat(limma_fit, coef = ncol(limma_design), number = nrow(obs_raw_filt))
limma_res <- limma_res %>%
  mutate(target_id = rownames(limma_res)) %>%
  dplyr::select(target_id, pval = P.Value, qval = adj.P.Val)

################################################################################
# edgeR
################################################################################

dge <- edgeR::DGEList(counts = obs_raw_filt, group = s2c$condition)
y <- edgeR::calcNormFactors(dge)
# y@.Data[[2]]$norm.factors <- 1 / so$est_counts_sf
y <- edgeR::estimateCommonDisp(y)
y <- edgeR::estimateTagwiseDisp(y)
et <- edgeR::exactTest(y)
edgeR_res <- edgeR::topTags(et, n = nrow(obs_raw_filt)) %>%
  as.data.frame() %>%
  mutate(., target_id = rownames(.)) %>%
  dplyr::select(target_id, pval = PValue, qval = FDR)

########################################################################
# benchmarks
########################################################################
library('mamabear')

de_bench <- new_de_benchmark(
  list(
    deseq_res,
    limma_res,
    edgeR_res,
    sres,
    s_ratio
    ),
  c(
    "DESeq2",
    "voom",
    "edgeR (tagwise)",
    "sleuth_wt",
    "sleuth_lrt"
    ), de_info)

fdr_nde_plot(de_bench, FALSE) +
  xlim(0, 1500) +
  ylim(0, 0.10) +
  xlab('number of transcripts called DE') +
  ylab('FDR') +
  theme(legend.position = c(0.1, 0.80))
# Saving 12.7 x 5.76 in image
ggsave('../img/3_3_1_1_2.pdf')

ggsave('../img/3_3_1_1_2.png')

fdr_nde_plot(de_bench, FALSE) +
  #theme_bw(25) +
  xlim(0, 5000) +
  #xlim(5000, 7500) +
  #xlim(0, 7500) +
  #ylim(0, 0.2) +
  ylim(0, 0.5) +
  theme(legend.position = c(0.1, 0.85))

fdr_nde_plot(de_bench, TRUE) +
  #theme_bw(25) +
  xlim(0, 2000) +
  #xlim(5000, 7500) +
  #xlim(0, 7500) +
  #ylim(0, 0.2) +
  ylim(0, 0.15) +
  theme(legend.position = c(0.1, 0.85))
