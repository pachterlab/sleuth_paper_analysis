library("cowplot")
library("sleuth")

base_dir <- "../results/norm/3_3_1_1"

# de_info contains the 'true' fold change as well as which transcripts are DE
de_info <- read.table(gzfile(file.path(base_dir, "de_info.tsv.gz")),
  header = TRUE, stringsAsFactors = FALSE)

kal_fnames <- file.path(base_dir, "sample_1", 1:6, "kallisto", "abundance.h5" )

kal <- lapply(kal_fnames, read_kallisto_h5, read_bootstrap = TRUE)

s2c <- data.frame(sample = paste0("sample_", 1:6),
  condition = rep(c("A", "B"), each = 3), stringsAsFactors = FALSE)

sobj <- new_sleuth(kal, s2c, ~ condition)
rm(kal)
gc()

sobj$obs_norm <- sobj$obs_raw


tmp_summary <- sleuth:::sleuth_summarize_bootstrap_col(sobj, "est_counts",
  function(x) log(x + 0.5))

bs_summary <- bs_sigma_summary(sobj, function(x) log(x + 0.5))

mes <- me_model_by_row(sobj, sobj$design_matrix, bs_summary)

# TODO:
pass_filt <- sobj$obs_norm %>%
  group_by(target_id) %>%
  summarise( count_filt = sum(est_counts > 5) > 0.8 * n())

mes_df <- rbindlist(mes) %>% as.data.frame
mes_df$target_id <- names(mes)
mes_df <- semi_join(mes_df, filter(pass_filt, count_filt), by = "target_id")

mes_df <- mes_df %>%
  mutate(sd_rate = sqrt(var_obs)/sqrt(sigma_q_sq)) %>%
  mutate(sd_rate_ecdf = ecdf(sd_rate)(sd_rate)) %>%
  mutate(sd_rate_group = cut(sd_rate_ecdf, 4)) %>%
  mutate(sigma_sq_pmax = pmax(sigma_sq, 0))

hi <- sliding_window_grouping(mes_df, "mean_obs", "sigma_sq_pmax", ignore_zeroes = TRUE)


loess_smooth <- shrink_df(hi, sqrt(sqrt(sigma_sq_pmax)) ~ mean_obs, "iqr")
loess_smooth <- loess_smooth^4
hi <- hi %>%
  mutate(smooth_sigma_sq = loess_smooth)
hi <- hi %>%
  mutate(smooth_max_sigma_sq = pmax(smooth_sigma_sq, sigma_sq))

# plotting mean log expression vs sigma
ggplot(hi, aes(mean_obs, sqrt(sqrt(sigma_sq_pmax)))) +
  geom_point(aes(colour = iqr), alpha = 0.2) +
  geom_line(aes(mean_obs, smooth_sigma_sq ^ (1/4))) +
  scale_colour_manual(values = c("black", "dodgerblue")) +
  xlim(1.5, 12) +
  ylim(0, 1.5) +
  theme(legend.position = "none") +
  xlab("mean( log( counts + 0.5 ) )") +
  ylab("sqrt( sigma )")

x_val <- sobj$design_mat[,"conditionB"]
Sxx <- sum( ( x_val - mean(x_val) ) ^ 2)
#debugonce(compute_t_me)
hi <- as.data.frame(hi)
ah <- compute_t_me(hi, "smooth_sigma_sq", Sxx, 6)
ah <- ah %>%
  select(target_id, pval) %>%
  mutate(qval = p.adjust(pval, method = "BH"))
hist(ah$pval)

ah_adj <- compute_t_me(hi, "smooth_sigma_sq", Sxx, 8, 1)
ah_adj <- ah_adj %>%
  select(target_id, pval) %>%
  mutate(qval = p.adjust(pval, method = "BH"))

ah2 <- compute_t_me(hi, "smooth_max_sigma_sq", Sxx, 8, 0.1) %>%
  select(target_id, pval) %>%
  mutate(qval = p.adjust(pval, method = "BH"))

########################################################################

pass_filt_names <- filter(pass_filt, count_filt)$target_id

################################################################################
# DESeq 2
#library("DESeq2")


obs_raw <- spread_abundance_by(sobj$obs_raw, "est_counts")
obs_raw_filt <- obs_raw[pass_filt_names,]

dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(obs_raw_filt),
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
# TODO: limma
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
  select(target_id, pval = P.Value, qval = adj.P.Val)

################################################################################
# edgeR
################################################################################

dge <- edgeR::DGEList(counts = obs_raw_filt, group = s2c$condition)
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


de_bench <- new_de_benchmark(
  list(
    deseq_res,
    ah2,
    limma_res,
    edgeR_res,
    edgeR_trended_res
    ),
  c(
    "DESeq2",
    "sleuth (ME)",
    "voom",
    "edgeR (tagwise)",
    "edgeR (trended)"
    ), de_info)


fdr_nde_plot(de_bench, FALSE) +
  #theme_bw(25) +
  xlim(1000, 11000) +
  #xlim(5000, 7500) +
  #xlim(0, 7500) +
  ylim(0, 0.20) +
  theme(legend.position = c(0.1, 0.85))

fdr_nde_plot(de_bench) +
  xlim(0, 12000) +
  #xlim(5000, 7500) +
  #xlim(0, 7500) +
  ylim(0, 0.250)



debugonce(de_rank_scatter)
gah <- de_rank_scatter(de_bench, 5000)

ggplot(gah, aes(relative_rank, est_rank)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~method)

gah %>%
  do(function(x)
    {
      x$rank_diff

    }(.))

gah <- gah %>%
  mutate(rank_diff = true_rank - est_rank,
    rel_rank_diff = relative_rank - est_rank)

ggplot(gah, aes(true_rank, est_rank)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~method)

ggplot(gah, aes(rank_diff)) +
  geom_histogram() +
  facet_wrap(~ method)


gah %>%
  filter(grepl("DESeq2", method) | grepl("sleuth", method)) %>%
  ggplot(aes(rank_diff)) +
    geom_density(aes(group = method, fill = method), alpha = 0.2)

gah %>%
  filter(grepl("DESeq2", method) | grepl("sleuth", method)) %>%
  ggplot(aes(rel_rank_diff)) +
    geom_density(aes(group = method, fill = method), alpha = 0.2)

document("~/dev/sleuth")
install("~/dev/sleuth")
