library("ggplot2")
library("cowplot")

# library("DESeq2")

load("../results/fin_females.RData", verbose = TRUE)

se_fin_females <- GenomicRanges::SummarizedExperiment(round(fin_females))
dds_fin_females <- DESeq2::DESeqDataSet(se_fin_females, ~ 1)
dds_fin_females <- DESeq2::estimateSizeFactors(dds_fin_females)
dds_fin_females <- DESeq2::estimateDispersionsGeneEst(dds_fin_females)
dds_fin_females <- DESeq2::estimateDispersionsFit(dds_fin_females)

#' Prepare a DESeqDataSet for simulation
#'
#' A basic routine that filters and adjusts dispersion estimates so that we can
#' simulate from it
#'

#' @param dds a DESeqDataSet which has \code{dispGeneEst}
#' @param min_mean the minumum mean counts to consider for differential
#' expression
#' @param min_dispersion the minimum dispersion to use. Default pulled from
#' DESeq2paper package
#' @return a data.frame that with can simulate from. Adds column sim_filt to
#' denote whether or target should be considered for DE.
#' @export
prep_dds_sim <- function(dds, min_mean = 5, min_dispersion = 1e-6) {
  fit <- as.data.frame(GenomicRanges::mcols(dds))
  fit <- fit %>%
    mutate(target_id = names(GenomicRanges::rowData(dds)))

  # we want to keep dispersion estimates within a reasonable range.
  # range was extracted from DESeq2paper package
  fit <- fit %>%
    mutate(disp_filt = dispGeneEst > min_dispersion) %>%
    mutate(disp_filt = ifelse(is.na(disp_filt), FALSE, disp_filt))

  # for dispersion estimates that are hard to estimate (due to low counts),
  # simply use a fixed estimate based on those which pass the filter
  med_dispersion <- fit %>%
    filter(disp_filt) %>%
    .$dispGeneEst %>%
    median

  # we now have dispersion estimates for every single transcript
  fit <- fit %>%
    mutate(disp_final = ifelse(disp_filt, dispGeneEst, med_disp))

  # we don't want to simulate DE from transcripts that are unlikely expressed
  # things, so let's omit transcripts that are zero in every condition
  fit <- fit %>%
    mutate(sim_filt = !allZero & baseMean > min_mean)

  fit
}

# the data.frame to simulate from
prep_fin <- prep_dds_sim(dds_fin_females)

save(prep_fin, file = "../results/prep_fin.RData")

########################################################################
# exploratory analysis
########################################################################

fin_df_sim %>%
  filter(sim_filt) %>%
  ggplot(aes(baseMean, disp_final)) +
  geom_point(alpha = 0.4) +
  coord_trans(xtrans = "log", ytrans = "log")

summary(fin_df_sim)


filter(fin_females_fit, dispGeneEst > 1e-6) %>%
  ggplot(aes(baseMean, dispGeneEst)) +
  geom_point(alpha = 0.2) +
  coord_trans(xtrans = "log2", ytrans = "log2") +
  xlim(1e-2, 1e4) +
  ylim(1e-8, 1e1)

plotDispEsts(dds_fin_females, ylim = c(1e-8, 50))

# get dispersion estimates for
fin_females_filt <- fin_females_fit %>%
  filter(dispGeneEst > 1e-6)

unfilt <- anti_join(fin_females_fit, fin_females_filt,
  by = "target_id")

unfilt %>%
  filter(!allZero) %>%
  ggplot(aes(baseMean)) +
  geom_histogram(binwidth=0.1)

fin_females_filt %>%
  ggplot(aes(dispGeneEst)) +
  geom_histogram(binwidth=0.1) +
  xlim(0, 15)

fin_females_filt %>%
  ggplot(aes(1,dispGeneEst)) +
  geom_boxplot()
