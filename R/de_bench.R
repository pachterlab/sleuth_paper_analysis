document("~/dev/sleuth")
install("~/dev/sleuth")

################################################################################
# DESeq 2
library("DESeq2")


obs_raw <- spread_abundance_by(s_o$obs_raw, "est_counts")

dds <- DESeqDataSetFromMatrix(countData = round(obs_raw),
  colData = select(id_to_condition, condition),
  design = ~ condition)

# force to not estimate size factors
sizeFactors(dds) <- rep(1, 8)

dds <- DESeq(dds, fitType='local')
deseq_res <- results(dds)

deseq_res <- as.data.frame(deseq_res) %>%
  mutate(target_id = rownames(deseq_res)) %>%
  dplyr::rename(qval = padj, pval = pvalue)

################################################################################
# TODO: limma
################################################################################


################################################################################

# de_bench <- new_de_benchmark(
#   list(deseq_res, pvals, valid_obs_dat),
#   c("deseq2", "lmm", "glm"), de_truth)

de_bench <- new_de_benchmark(
  list(deseq_res, pvals, pvals_shrink),
  c("deseq2", "lmm", "lmm_naive_shrink"), de_truth)

fdr_tpr_plot(de_bench)

fdr_nde_plot(de_bench) +
  theme_bw(25) +
  xlim(0,12000) +
  ylim(0,0.250)

ggsave("../img/fdr_naive_shrink.png")

ggsave("../img/fdr_nde.png")


pval_distribution(de_bench)

