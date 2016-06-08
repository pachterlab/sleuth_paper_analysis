# temporary for debugging
n_cpu <- 20
sim_name <- 'isoform_3_3_20_1_1'
# sim_name <- 'gfr_3_3_20_42_2'

n_cpu <- args[1]
sim_name <- args[2]


source("benchmark_methods.R")
source("gene_common.R")

library('cowplot')
library("parallel")
library("mamabear")
options(mc.cores = n_cpu)

transcript_gene_mapping <- get_human_gene_names()

###
# run each method on their own filter
###

sim <- parse_simulation(sim_name)

sample_info <- lapply(1:1,
  function(i) {
    n <- sim$a + sim$b
    kal_dirs <- file.path('..', 'sims', sim_name, paste0("exp_", i),
      1:n, "kallisto")
    get_sample_to_condition(sim$a, sim$b, kal_dirs)
  })

sleuth_results <- lapply(1:1,
  function(i) {
    si <- sample_info[[i]]
    res <- run_sleuth(si)
    res$counts <- sleuth:::spread_abundance_by(res$so$obs_raw, "est_counts")
    res
  })

all_counts <- lapply(sleuth_results, '[[', 'counts')

dummy_count_filter <- rep(TRUE, nrow(all_counts[[1]]))
names(dummy_count_filter) <- rownames(all_counts[[1]])

each_filter <- list()

# The code below is a slightly modified version of the code from `DESeq2paper`
# http://www-huber.embl.de/DESeq2paper/
runDESeq2_all <- function(e, as_gene = TRUE, compute_filter = FALSE) {
  dds <- DESeqDataSetFromMatrix(exprs(e), DataFrame(pData(e)), ~ condition)
  if (compute_filter) {
    # Section 1.3.6 in DESeq2 vignette
    # https://www.bioconductor.org/packages/3.3/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
    dds <- dds[ rowSums(counts(dds)) > 1, ]
  }
  dds <- DESeq(dds,quiet=TRUE)
  res <- results(dds)
  res <- as.data.frame(res)
  # beta <- res$log2FoldChange
  # pvals <- res$pvalue
  # padj <- res$padj
  # pvals[is.na(pvals)] <- 1
  # # pvals[rowSums(exprs(e)) == 0] <- NA
  # padj[is.na(padj)] <- 1

  # rename_target_id(
  #   data.frame(target_id = rownames(res),
  #     pval = pvals, qval = padj, beta = beta,
  #     stringsAsFactors = FALSE),
  #   as_gene = as_gene)
  res
}

DESeq2_filter_and_run_intersect_all <- function(count_matrix, stc, sleuth_filter) { # nolint
  count_matrix <- round(count_matrix)
  mode(count_matrix) <- 'integer'
  which_targets <- DESeq2_filter(count_matrix)
  sleuth_filter <- sleuth_filter & which_targets
  cds <- make_count_data_set(count_matrix[sleuth_filter, ], stc)
  res <- runDESeq2_all(cds, FALSE, FALSE)

  sleuth_filter <- names(which(sleuth_filter))

  list(result = res, filter = sleuth_filter)
}

each_filter$DESeq2 <- mclapply(1:1,
  function(i) {
    si <- sample_info[[i]]
    counts <- all_counts[[i]]
    DESeq2_filter_and_run_intersect_all(counts, si, dummy_count_filter)$result
  })

i <- 1
si <- sample_info[[i]]
counts <- all_counts[[i]]
cds <- make_count_data_set(round(counts),si)
dds <- DESeqDataSetFromMatrix(exprs(cds), DataFrame(pData(cds)), ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1,]
dds <- DESeq(dds)

library(sleuth)
plotDispEsts(dds)
plot_mean_var(sleuth_results[[1]]$so)

res <- results(dds)
res <- as.data.frame(res)


mu <- mcols(dds)$baseMean
sel <-(mu>0)
mu <- mu[sel]
phi <- mcols(dds)$dispGeneEst[sel]
phi <- dispersions(dds)[sel]

cv_deseq2 <- phi + 1/mu

sr <- sleuth::sleuth_results(sleuth_results[[1]]$so, 'reduced:full', 'lrt')
cv_sleuth <- sleuth_results$final_sigma_sq/sr$mean_obs^2

# sigma2 <- mu + mu^2 * phi
#
log_mu <- function(mu,sigma2) {log(mu) - sigma2/(2 * mu^2)}
log_sigma2 <- function(mu,sigma2) {sigma2/(mu^2) - 0.25 * sigma2^2/(mu^4)}
#
# lmu <- log_mu(mu,sigma2)
# ls2 <- log_sigma2(mu,sigma2)


save(sample_info, all_counts, dds, sleuth_results,file='junk.RData')
load('junk.RData')
sims <- readRDS('../sims/isoform_3_3_20_1_1/sims.rds')
sim1 <- sims[[1]]

ds_names <- dds@rowRanges@partitioning@NAMES
ds_df <- data.frame(target_id=ds_names, mu=mcols(dds)$baseMean, phi=dispersions(dds))

ds_df <- inner_join(ds_df, dplyr::select(sim1$info,target_id,baseMean, phi=dispGeneEst))

sl_df <- inner_join(sr, dplyr::select(sim1$info,target_id,baseMean,phi=dispGeneEst))
sl_df <- dplyr::mutate(sl_df, cv = (phi + 1 / baseMean),
  cv_est = (final_sigma_sq + tech_var) / mean_obs ^ 2)
sl_df <- dplyr::mutate(sl_df, cv_inferential = tech_var / mean_obs ^ 2,
  cv_expected = 1 / mean_obs)

ggplot(sl_df, aes(mean_obs, cv_inferential)) +
  geom_point(alpha = 0.2) +
  geom_line(aes(mean_obs, cv_expected), color = 'red', size = 1.2) +
  theme_cowplot(25) +
  coord_trans(y = 'log') +
  xlab('mean( log(x + 0.5) )') +
  ylab('squared coefficient of variation')
