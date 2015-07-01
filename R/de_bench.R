document("~/dev/sleuth")
install("~/dev/sleuth")

pass_filt <- s_o$obs_norm %>%
  group_by(target_id) %>%
  summarise( count_filt = sum(est_counts > 5) > 0.8 * n())

pass_filt_names <- filter(pass_filt, count_filt)$target_id

################################################################################
# DESeq 2
#library("DESeq2")


obs_raw <- spread_abundance_by(s_o$obs_raw, "est_counts")
obs_raw_filt <- obs_raw[pass_filt_names,]

dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(obs_raw_filt),
  colData = select(id_to_condition, condition),
  design = ~ condition)

# force to not estimate size factors
DESeq2::sizeFactors(dds) <- rep(1, 8)

#dds <- DESeq2::DESeq(dds, fitType='local')
dds <- DESeq2::DESeq(dds)
deseq_res <- DESeq2::results(dds)

deseq_res <- as.data.frame(deseq_res) %>%
  mutate(target_id = rownames(deseq_res)) %>%
  dplyr::rename(qval = padj, pval = pvalue)

################################################################################
# TODO: limma
################################################################################

dge <- edgeR::DGEList(counts = obs_raw_filt)

# not normalizing because something odd is happening
# ddge <- edgeR::calcNormFactors(dge)

limma_design <- model.matrix(~ condition, id_to_condition)
v <- limma::voom(dge, limma_design, plot = TRUE)

limma_fit <- limma::lmFit(v, limma_design)
limma_fit <- limma::eBayes(limma_fit)
limma_res <- limma::topTreat(limma_fit, coef = ncol(limma_design), number = nrow(obs_raw_filt))

limma_res <- limma_res %>%
  mutate(target_id = rownames(limma_res)) %>%
  select(target_id, pval = P.Value, qval = adj.P.Val)

################################################################################
# edgeR
################################################################################

dge <- edgeR::DGEList(counts = obs_raw_filt, group = id_to_condition$condition)
y <- edgeR::estimateCommonDisp(dge)
y <- edgeR::estimateTagwiseDisp(y)
et <- edgeR::exactTest(y)
edgeR_res <- edgeR::topTags(et, n = nrow(obs_raw_filt)) %>%
  as.data.frame() %>%
  mutate(., target_id = rownames(.)) %>%
  select(target_id, pval = PValue, qval = FDR)

z <- edgeR::estimateCommonDisp(dge)
z <- edgeR::estimateTrendedDisp(z)
et_z <- edgeR::exactTest(z)
edgeR_trended_res <- edgeR::topTags(et_z, n = nrow(obs_raw_filt)) %>%
  as.data.frame() %>%
  mutate(., target_id = rownames(.)) %>%
  select(target_id, pval = PValue, qval = FDR)

################################################################################

# de_bench <- new_de_benchmark(
#   list(deseq_res, pvals, valid_obs_dat),
#   c("deseq2", "lmm", "glm"), de_truth)

# pvals <- pvals %>%
#   mutate(qval = p.adjust(pval, method = "BH"))
# summary_lmm_ml_df <- summary_lmm_ml_df %>%
#   mutate(qval = p.adjust(pval, method = "BH"))

de_bench <- new_de_benchmark(
  list(
    deseq_res,
    # pvals,
    #naive_shrink_pvals,
    #summary_lmm_ml_df,
    # raw_pvals,
    #group_shrink_pval,
    #group_shrink_obs_pval,
    #ah,
    #ah_adj,
    ah2,
    #limma_res,
    edgeR_res
    #edgeR_trended_res
    #convex_pval
    ),
  c(
    "DESeq2",
    # "lmm_reml",
    #"naive_shrink",
    #"lmm_ml",
    # "hp_raw",
    #"group_shrink",
    #"group_shrink_obs",
    #"me",
    #"me_adj",
    "sleuth (ME)",
    #"voom",
    "edgeR (tagwise)"
    #"edgeR (trended)"
    #"convex_pval"
    ), de_truth)

fdr_tpr_plot(de_bench) +
  xlim(0, 0.2) +
  ylim(0.75, 1.0)

fdr_nde_plot(de_bench) +
  xlim(0, 12000) +
  #xlim(5000, 7500) +
  #xlim(0, 7500) +
  ylim(0, 0.250)

debugonce(fdr_nde_plot)

fdr_nde_plot(de_bench, FALSE) +
  #theme_bw(25) +
  xlim(1000, 11000) +
  #xlim(5000, 7500) +
  #xlim(0, 7500) +
  ylim(0, 0.20) +
  theme(legend.position = c(0.1, 0.85))

ggsave("../img/fdr_me_edger_deseq.pdf")

fdr_nde_plot(de_bench, TRUE) +
  #theme_bw(25) +
  xlim(1000, 11000) +
  #xlim(5000, 7500) +
  #xlim(0, 7500) +
  ylim(0, 0.20) +
  theme(legend.position = c(0.1, 0.85))
ggsave("../img/efdr_me_edger_deseq.pdf")

ggsave("../img/fdr_nde.png")


pval_distribution(de_bench)

