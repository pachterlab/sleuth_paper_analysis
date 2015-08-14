document("~/dev/sleuth")
install("~/dev/sleuth")

roxygenise("~/dev/sleuth")


library("cowplot")
library("sleuth")

base_dir <- "../results/3_3_1_1"

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

pass_filt <- sobj$obs_norm %>%
  group_by(target_id) %>%
  summarise( count_filt = sum(est_counts > 5) > 0.8 * n())

# pass_filt_2 <- sobj$obs_norm %>%
#   group_by(target_id) %>%
#   summarise( count_filt = mean(est_counts) > 3)

########################################################################
# pooled variance and technical variance estimator
########################################################################

deg_free <- 4

debugonce(me_equal_var)
mev <- me_equal_var(sobj, pass_filt, function(x) log(x + 0.5))

X <- sobj$design_matrix
A <- solve( t(X) %*% X )

ols_beta_covars <- lapply(mev$summary$rss,
  function(i) {
    (i / deg_free) * A
  })

all_se <- data.frame(
  target_id = names(mev$beta_covars),
  ols_var = sapply(ols_beta_covars, function(x) x[2,2]),
  me_var = sapply(mev$beta_covars, function(x) x[2,2])
  )

m_ols <- with(all_se, median(ols_var, na.rm = TRUE))
m_me <- with(all_se, median(me_var, na.rm = TRUE))
m_diff <- m_me/m_ols
all_se <- mutate(all_se, me_var_adj = ifelse(me_var == ols_var, me_var,
# all_se <- mutate(all_se, me_var_adj = ifelse(FALSE, me_var,
    me_var / m_diff))

b <- sapply(mev$mes,
  function(x) {
    x$ols_fit$coefficients[2]
  })
names(b) <- names(mev$mes)
b <- data.frame(target_id = names(b), b = b)

b <- inner_join(all_se, b, by = "target_id")

b_adj <- mutate(b, t_value = b / sqrt(me_var_adj))
b_adj <- mutate(b_adj, pval = 2 * pt(abs(t_value), deg_free, lower.tail = FALSE))
b_adj <- mutate(b_adj, qval = p.adjust(pval, method = "BH"))

ggplot(all_se, aes(log2(ols_se), log2(me_se))) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, colour = "blue")

# looking at only shrinking things with sigma > 0

mev$summary <- mutate(mev$summary, sigma_sq_gt0_shrink =
  ifelse(sigma_sq > 0, pmax(sigma_sq, smooth_sigma_sq), 0))

sigma_gt0_covars <- with(mev$summary, lapply(sigma_q_sq + sigma_sq_gt0_shrink,
    function(x) {
      x * A
    }))
all_se$gt0_var <- sapply(sigma_gt0_covars, function(x) x[2,2])

b_gt0 <- b
b_gt0$gt0_var <- all_se$gt0_var

b_gt0 <- mutate(b_gt0, t_value = b / sqrt(gt0_var))
b_gt0 <- mutate(b_gt0, pval = 2 * pt(abs(t_value), deg_free, lower.tail = FALSE))
b_gt0 <- mutate(b_gt0, qval = p.adjust(pval, method = "BH"))

# basic equal variance model
mev_b <- sapply(mev$beta_covars, function(x) sqrt(x[2,2]))
mev_b <- data.frame(target_id = names(mev_b), se_b = mev_b,
  stringsAsFactors = FALSE)
rownames(mev_b) <- NULL
mev_b$b <- mev$summary$b1

mev_b <- inner_join(mev_b, b, by = "target_id")

mev_b <- mutate(mev_b, t_value = b / se_b)
mev_b <- mutate(mev_b, pval = 2 * pt(abs(t_value), deg_free,
    lower.tail = FALSE))
mev_b <- mutate(mev_b, qval = p.adjust(pval, method = "BH"))

bs_summary <- bs_sigma_summary(sobj, function(x) log(x + 0.5))

mes <- me_model_by_row(sobj, sobj$design_matrix, bs_summary)

mes_df <- rbindlist(mes) %>% as.data.frame
mes_df$target_id <- names(mes)
mes_df <- semi_join(mes_df, filter(pass_filt, count_filt), by = "target_id")

mes_df <- mes_df %>%
  mutate(sd_rate = sqrt(var_obs)/sqrt(sigma_q_sq)) %>%
  mutate(sd_rate_ecdf = ecdf(sd_rate)(sd_rate)) %>%
  mutate(sd_rate_group = cut(sd_rate_ecdf, 4)) %>%
  mutate(sigma_sq_pmax = pmax(sigma_sq, 0))

hi <- sliding_window_grouping(mes_df, "mean_obs", "sigma_sq_pmax", ignore_zeroes = TRUE)

l_smooth <- shrink_df(hi, sqrt(sqrt(sigma_sq_pmax)) ~ mean_obs, "iqr")
l_smooth <- mutate(l_smooth, smooth_sigma_sq = shrink ^ 4) %>%
  select(-shrink)
hi <- l_smooth %>%
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
ah <- compute_t_me(hi, "smooth_sigma_sq", Sxx, 8)



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
# compare the two diff versions of computing covar(b)
########################################################################

plot(ah2$se_b, b_sd$sd_b)
cor(ah2$se_b, b_sd$sd_b, use = "complete.obs")
cor(ah2$se_b, b_sd$sd_b, use = "complete.obs", method = "spearman")

########################################################################

pass_filt_names <- filter(pass_filt, count_filt)$target_id

#pass_filt_names <- filter(pass_filt_2, count_filt)$target_id


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

saveRDS(deseq_res, "../tmp/DESeq2.rds")

deseq_res <- readRDS("../tmp/DESeq2.rds")

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

saveRDS(limma_res, "../tmp/limma.rds")

limma_res <- readRDS("../tmp/limma.rds")

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
  dplyr::select(target_id, pval = PValue, qval = FDR)
saveRDS(edgeR_res, "../tmp/edgeR.rds")

edgeR_res <- readRDS("../tmp/edgeR.rds")

# z <- edgeR::estimateCommonDisp(dge)
# z <- edgeR::estimateTrendedDisp(z)
# et_z <- edgeR::exactTest(z)
# edgeR_trended_res <- edgeR::topTags(et_z, n = nrow(obs_raw_filt)) %>%
#   as.data.frame() %>%
#   mutate(., target_id = rownames(.)) %>%
#   dplyr::select(target_id, pval = PValue, qval = FDR)


# de_bench <- new_de_benchmark(
#   list(
#     deseq_res,
#     ah2,
#     limma_res,
#     edgeR_res,
#     #edgeR_trended_res,
#     betas
#     ),
#   c(
#     "DESeq2",
#     "sleuth (ME)",
#     "voom",
#     "edgeR (tagwise)",
#     #"edgeR (trended)",
#     "sleuth (White)"
#     ), de_info)

de_bench <- new_de_benchmark(
  list(
    deseq_res,
    #mev_b,
    b_adj,
    limma_res,
    edgeR_res,
    b_gt0
    #edgeR_trended_res,
    ),
  c(
    "DESeq2",
    "sleuth (MEV)",
    "voom",
    "edgeR (tagwise)",
    "sleuth (MEV gt0)"
    #"edgeR (trended)",
    #"sleuth (White)"
    ), de_info)

fdr_nde_plot(de_bench, TRUE) +
  #theme_bw(25) +
  xlim(2000, 7500) +
  #xlim(5000, 7500) +
  #xlim(0, 7500) +
  #ylim(0, 0.2) +
  ylim(0, 1) +
  theme(legend.position = c(0.1, 0.85))

fdr_nde_plot(de_bench, FALSE) +
  #theme_bw(25) +
  xlim(2000, 7500) +
  #xlim(5000, 7500) +
  #xlim(0, 7500) +
  ylim(0, 0.2) +
  theme(legend.position = c(0.1, 0.85))

# Saving 12.8 x 7.52 in image
ggsave("../img/de_bench.pdf")

########################################################################
# unequal variances
########################################################################

sample_bs_summary <- sleuth:::sleuth_summarize_bootstrap_col(sobj, "est_counts",
  function(x) log(x + 0.5))

# taking a look at the variability of the variability
samp_bs_cv <- sample_bs_summary %>%
  group_by(target_id) %>%
  summarise(
    var_mean = mean(bs_var_est_counts),
    var_se = sd(bs_var_est_counts)) %>%
  mutate(cv = var_se / var_mean)

samp_vs_cv <- inner_join(samp_bs_cv,
  dplyr::select(de_info, target_id, is_de, true_mean = baseMean),
  by = "target_id")

samp_vs_cv <- inner_join(samp_vs_cv, pass_filt, by = "target_id")

ggplot(samp_vs_cv, aes(log(true_mean + 0.5), cv)) +
  geom_point(aes(colour = count_filt), alpha = 0.2) +
  facet_wrap(~ is_de)


# cnts <- bs_summary$obs_counts['ENST00000000233',]
# bs_var <- sample_bs_summary %>% filter(target_id == "ENST00000000233") %>% .$bs_var_est_counts
# X <- sobj$design_matrix
# A <- solve(t(X) %*% X) %*% t(X)
# res <- sleuth:::me_white_model(sobj$design_matrix, cnts, bs_var, A)

debugonce(me_heteroscedastic_by_row)
tmp <- me_heteroscedastic_by_row(sobj, sobj$design_matrix, sample_bs_summary, bs_summary$obs_counts)

sig_summary <- lapply(tmp, function(x) x$df) %>% bind_rows

sig_summary <- inner_join(sig_summary, pass_filt, by = "target_id") %>%
  inner_join(dplyr::select(de_info, target_id, is_de), by = "target_id")

sig_summary %>%
  filter(count_filt) %>%
  ggplot(aes(mean_obs, pmax(0, sigma_sq) ^ (1/4))) +
  geom_point(alpha = 0.05) +
  facet_wrap(~ sample)


save.image(file = "stuff.RData")

load("stuff.RData", verbose = TRUE)

sig_summary <- sig_summary %>%
  mutate(sigma_sq_pmax = pmax(sigma_sq, 0))
sig_summary_filt <- semi_join(sig_summary, filter(pass_filt, count_filt), by = "target_id")

swg <- sig_summary_filt %>%
  group_by(sample) %>%
  do({
    sliding_window_grouping(., "mean_obs", "sigma_sq_pmax", ignore_zeroes = TRUE)
  })


system.time(loess_smooth <- swg %>%
  group_by(sample) %>%
  do(shrink_df(., sqrt(sqrt(sigma_sq_pmax)) ~ mean_obs, "iqr")))

loess_smooth <- loess_smooth %>%
  mutate(sigma_sq_shrink = shrink ^ 4) %>%
  select(-shrink)

# plotting mean log expression vs sigma
loess_smooth %>%
  ggplot(aes(mean_obs, sqrt(sqrt(sigma_sq_pmax)))) +
  geom_point(aes(colour = iqr), alpha = 0.2) +
  geom_line(aes(mean_obs, sigma_sq_shrink ^ (1/4)), colour = "red") +
  scale_colour_manual(values = c("black", "dodgerblue")) +
  theme(legend.position = "none") +
  xlab("mean( log( counts + 0.5 ) )") +
  ylab("sqrt( sigma )") +
  facet_wrap(~ sample)

loess_smooth <- loess_smooth %>%
  mutate(sigma_sq_smooth = pmax(sigma_sq_pmax, sigma_sq_shrink, na.rm = TRUE))

junk1 <- filter(loess_smooth, target_id == 'ENST00000226413')

X <- sobj$design_matrix
A <- solve(t(X) %*% X)

debugonce(sleuth:::me_white_test)
j1 <- sleuth:::me_white_covar(ungroup(junk1), 'sigma_sq_smooth', 'sigma_q_sq', X, A)

debugonce(sleuth:::me_white_var)

beta_var <- do(group_by(ungroup(loess_smooth), target_id),
  sleuth:::me_white_var(., 'sigma_sq_smooth', 'sigma_q_sq', X, A))

betas <- lapply(tmp, function(x) as.data.frame(t(x$ols$coefficients))) %>%
  rbind_all
betas$target_id <- names(tmp)

beta_var <- beta_var %>%
  mutate(sd_conditionB = sqrt(conditionB)) %>%
  select(target_id, sd_conditionB)
beta_var <- ungroup(beta_var)

betas <- inner_join(betas, beta_var, by = "target_id")
betas <- betas %>%
  mutate(t_stat = conditionB / sd_conditionB)

deg_free <- 4
betas <- mutate(betas,
  pval = 2 * pt(abs(t_stat), 4, lower.tail = FALSE))
betas <- mutate(betas, qval = p.adjust(pval, method = "BH"))
betas <- mutate(betas, qval = pval)

########################################################################
# equal bio variance, but diff tech var
########################################################################

eq_sig_neq_tech <- swg %>%
  group_by(target_id) %>%
  summarise(sigma_sq_eq = {function(r, q){ mean(pmax(r - q, 0)) }}(r_sq, sigma_q_sq))

tmp_df <- sig_summary_filt %>%
  select(mean_obs, count_filt, target_id) %>%
  distinct

eq_sig_neq_tech <- inner_join(eq_sig_neq_tech, tmp_df, by = "target_id")

l_esnt <- shrink_df(
  sliding_window_grouping(eq_sig_neq_tech, "mean_obs", "sigma_sq_eq"),
  sqrt(sqrt(sigma_sq_eq)) ~ mean_obs,
  "iqr")
l_esnt <- select(mutate(l_esnt, sigma_sq_shrink = shrink ^ 4), -shrink)

l_esnt %>%
  ggplot(aes(mean_obs, sqrt(sqrt(sigma_sq_eq)))) +
  geom_point(aes(colour = iqr), alpha = 0.2) +
  geom_line(aes(mean_obs, sigma_sq_shrink ^ (1/4)), colour = "red") +
  scale_colour_manual(values = c("black", "dodgerblue")) +
  theme(legend.position = "none") +
  xlab("mean( log( counts + 0.5 ) )") +
  ylab("sqrt( sigma )")

l_esnt <- left_join(l_esnt, select(swg, target_id, sample, sigma_q_sq), by = c("target_id"))
l_esnt <- mutate(l_esnt, sigma_sq_shrink_pmax = pmax(sigma_sq_eq, sigma_sq_shrink))

beta_sd_esnt <- do(group_by(l_esnt, target_id),
  sleuth:::me_white_var(., 'sigma_sq_shrink_pmax', 'sigma_q_sq', X, A)) %>%
    mutate(sd_conditionB = sqrt(conditionB)) %>%
    select(target_id, sd_conditionB)

betas <- lapply(tmp, function(x) as.data.frame(t(x$ols$coefficients))) %>%
  bind_rows
betas$target_id <- names(tmp)

deg_free <- 4

beta_sd_esnt <- inner_join(betas, beta_sd_esnt, by = "target_id")
beta_sd_esnt <- mutate(beta_sd_esnt,
  t_stat = conditionB / sd_conditionB,
  pval = 2 * pt(abs(t_stat), 4, lower.tail = FALSE))
beta_sd_esnt <- mutate(beta_sd_esnt, qval = p.adjust(pval, method = "BH")

########################################################################
# bitseq
########################################################################

# load bitseq results
bitseq_res <- readRDS("../results/3_3_1_1/sample_1/bitseq_de.rds")
bitseq_df <- as.data.frame(bitseq_res$pplr)
bitseq_df$target_id <- rownames(bitseq_df)
bitseq_df <- select(bitseq_df, target_id, pplr.1.2, log2fc = log2FC.1.2)
bitseq_df$pde <- abs(0.5-bitseq_df$pplr)
bitseq_df <- bitseq_df %>%
  arrange(desc(pde), desc(abs(log2fc))) %>%
  mutate(rank = 1:nrow(.))
bitseq_df <- select(bitseq_df, target_id, pval = rank)
bitseq_df <- mutate(bitseq_df, qval = NA_real_)

########################################################################
# EBSeq with kallisto results
########################################################################

library("biomaRt")

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

gene_names <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)

rownames(gene_names) <- gene_names$ensembl_transcript_id
gene_names <- gene_names[rownames(obs_raw),]

library("EBSeq")

all.equal(rownames(gene_names), rownames(obs_raw))

ngs <- GetNg(rownames(gene_names), gene_names$ensembl_gene_id)
iso_ng_trunc <- ngs$IsoformNgTrun
iso_sizes <- rep(1, 6)
iso_eb_out_kal <- EBTest(Data = obs_raw,
  NgVector = iso_ng_trunc,
  Conditions = s2c$condition, sizeFactors = iso_sizes,
  maxround = 15)

iso_eb_kal <- as.data.frame(iso_eb_out_kal$PPMat)
iso_eb_kal$target_id <- rownames(iso_eb_kal)
iso_eb_kal <- iso_eb_kal %>%
  mutate(pval = 1 - PPDE, qval = pval)
iso_eb_kal <- dplyr::select(iso_eb_kal, target_id, pval, qval)

########################################################################
# EBSeq with RSEM results
########################################################################

rsem_fnames <- paste0("../results/3_3_1_1/sample_1/",1:6, "/rsem/out.isoforms.results")

rsem_res <- lapply(rsem_fnames, read.table, header = TRUE, stringsAsFactors = FALSE)

rsem_mat <- matrix(NA_real_, ncol = 6, nrow = nrow(rsem_res[[1]]))
for (i in 1:6) {
  rsem_mat[,i] <- rsem_res[[i]]$expected_count
}
rownames(rsem_mat) <- rsem_res[[1]]$transcript_id
colnames(rsem_mat) <- s2c$sample

ngs <- GetNg(rownames(rsem_mat), gene_names[rownames(rsem_mat),]$ensembl_gene_id)

iso_eb_out_rsem <- EBTest(Data = rsem_mat,
  NgVector = ngs$IsoformNgTrun,
  Conditions = s2c$condition, sizeFactors = iso_sizes,
  maxround = 15)

iso_eb_rsem <- as.data.frame(iso_eb_out_rsem$PPMat)
iso_eb_rsem$target_id <- rownames(iso_eb_rsem)
iso_eb_rsem <- iso_eb_rsem %>%
  mutate(pval = 1 - PPDE, qval = pval)
iso_eb_rsem <- dplyr::select(iso_eb_rsem, target_id, pval, qval)

########################################################################
# DESeq2 with RSEM
########################################################################

rsem_filt <- rsem_mat[pass_filt_names,]

run_deseq2 <- function(counts, samp2cond, size_factors = rep(1, ncol(counts))) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = dplyr::select(samp2cond, condition),
    design = ~ condition)

  DESeq2::sizeFactors(dds) <- size_factors

  dds <- DESeq2::DESeq(dds)
  deseq_res <- DESeq2::results(dds)

  deseq_res <- as.data.frame(deseq_res) %>%
    mutate(target_id = rownames(deseq_res)) %>%
    dplyr::rename(qval = padj, pval = pvalue)

  deseq_res
}


debugonce(run_deseq2)
deseq_rsem <- run_deseq2(rsem_filt, s2c)

########################################################################
# benchmarks
########################################################################

de_bench2 <- new_de_benchmark(
  list(
    deseq_res,
    ah2,
    betas,
    bitseq_df,
    iso_eb_kal,
    iso_eb_rsem,
    deseq_rsem
    ),
  c(
    "DESeq2 (kallisto)",
    "sleuth (ME)",
    "sleuth (white)",
    "BitSeq",
    "EBSeq (kallisto)",
    "EBSeq (RSEM)",
    "DESeq2 (RSEM)"
    ), de_info)


debugonce(fdr_nde_plot)

fdr_nde_plot(de_bench2, FALSE) +
  #theme_bw(25) +
  #xlim(1000, 11000) +
  #xlim(5000, 7500) +
  xlim(0, 7500) +
  ylim(0, 0.50) +
  theme(legend.position = c(0.1, 0.85))

fdr_nde_plot(de_bench2, TRUE) +
  #theme_bw(25) +
  #xlim(1000, 11000) +
  #xlim(5000, 7500) +
  xlim(0, 7500) +
  #ylim(0, 0.20) +
  ylim(0, 0.50) +
  theme(legend.position = c(0.1, 0.85))



########################################################################
# good candidates for heteroscedastic plots
########################################################################

samp_vs_cv %>% filter(count_filt) %>% arrange(desc(cv), desc(true_mean), desc(var_mean))

# ENST00000359635
# ENST00000356239
# ENST00000390619
# ENST00000492709

interesting <- c(
  "ENST00000359635",
  "ENST00000356239",
  "ENST00000390619",
  'ENST00000492709')


tmp <- loess_smooth %>%
  mutate(id = "none")

tmp <- tmp %>%
  mutate(id = ifelse(target_id %in% interesting, target_id, id))

# plotting mean log expression vs sigma
tmp %>%
  filter(id == 'none') %>%
  ggplot(aes(mean_obs, sqrt(sqrt(sigma_sq_pmax)))) +
  geom_point(colour = "black", alpha = 0.02) +
  geom_point(aes(colour = factor(id)), data = filter(tmp, id != 'none'), size = 4) +
  geom_line(aes(mean_obs, sigma_sq_shrink ^ (1/4)), colour = "red") +
  #scale_colour_manual(values = c("black", "dodgerblue")) +
  theme(legend.position = "none") +
  xlab("mean( log( counts + 0.5 ) )") +
  ylab("sqrt( sigma )") +
  facet_wrap(~ sample)
ggsave("../img/shrink_persample.png")
