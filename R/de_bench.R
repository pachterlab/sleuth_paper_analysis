document("~/dev/sleuth")
install("~/dev/sleuth")

################################################################################
# DESeq 2
library("DESeq2")


obs_raw <- spread_abundance_by(s_o$obs_raw, "est_counts")

dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(obs_raw),
  colData = select(id_to_condition, condition),
  design = ~ condition)

# force to not estimate size factors
DESeq2::sizeFactors(dds) <- rep(1, 8)

dds <- DESeq2::DESeq(dds, fitType='local')
deseq_res <- DESeq2::results(dds)

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

pvals <- pvals %>%
  mutate(qval = p.adjust(pval, method = "BH"))
summary_lmm_ml_df <- summary_lmm_ml_df %>%
  mutate(qval = p.adjust(pval, method = "BH"))

de_bench <- new_de_benchmark(
  list(
    deseq_res,
    # pvals,
    #naive_shrink_pvals,
    summary_lmm_ml_df,
    # raw_pvals,
    group_shrink_pval,
    group_shrink_obs_pval,
    convex_pval
    ),
  c(
    "deseq2",
    # "lmm_reml",
    #"naive_shrink",
    "lmm_ml",
    # "hp_raw",
    "group_shrink",
    "group_shrink_obs",
    "convex_pval"
    ), de_truth)

fdr_tpr_plot(de_bench) +
  xlim(0, 0.2) +
  ylim(0.75, 1.0)

fdr_nde_plot(de_bench) +
  theme_bw(25) +
  #xlim(0, 12000) +
  xlim(0, 7500) +
  ylim(0, 0.250)

ggsave("../img/fdr_naive_shrink_upd.png")

ggsave("../img/fdr_nde.png")


pval_distribution(de_bench)

