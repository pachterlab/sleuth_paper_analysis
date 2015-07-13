# library("ggplot2")
# library("cowplot")
source("../../simulation_core/R/simulate_de.R")

# load DF prep_fin which is ready to simulate from
load("../results/fin_females.RData", verbose = TRUE)

se_fin_females <- GenomicRanges::SummarizedExperiment(round(fin_females))
dds_fin_females <- DESeq2::DESeqDataSet(se_fin_females, ~ 1)
dds_fin_females <- DESeq2::estimateSizeFactors(dds_fin_females)
cat("Cox-Reid estimate\n")
dds_fin_females <- DESeq2::estimateDispersionsGeneEst(dds_fin_females)
cat("Smoothed estimate\n")
dds_fin_females <- DESeq2::estimateDispersionsFit(dds_fin_females)

# the data.frame to simulate from
cat("preparing the data.frame for DE simulations\n")
prep_fin <- prep_dds_sim(dds_fin_females)

cat("saving data.frame\n")
save(prep_fin, file = "../results/prep_fin.RData")

# sim_3_3 <- simulate_counts(prep_fin,
#   n_sim = 10L,
#   n_a = 3L,
#   n_b = 3L,
#   prop_de = 0.20,
#   seed = 42)

########################################################################
# exploratory analysis
########################################################################

# fin_df_sim %>%
#   filter(sim_filt) %>%
#   ggplot(aes(baseMean, disp_final)) +
#   geom_point(alpha = 0.4) +
#   coord_trans(xtrans = "log", ytrans = "log")
#
# summary(fin_df_sim)
#
#
# filter(fin_females_fit, dispGeneEst > 1e-6) %>%
#   ggplot(aes(baseMean, dispGeneEst)) +
#   geom_point(alpha = 0.2) +
#   coord_trans(xtrans = "log2", ytrans = "log2") +
#   xlim(1e-2, 1e4) +
#   ylim(1e-8, 1e1)
#
# plotDispEsts(dds_fin_females, ylim = c(1e-8, 50))
#
# # get dispersion estimates for
# fin_females_filt <- fin_females_fit %>%
#   filter(dispGeneEst > 1e-6)
#
# unfilt <- anti_join(fin_females_fit, fin_females_filt,
#   by = "target_id")
#
# unfilt %>%
#   filter(!allZero) %>%
#   ggplot(aes(baseMean)) +
#   geom_histogram(binwidth=0.1)
#
# fin_females_filt %>%
#   ggplot(aes(dispGeneEst)) +
#   geom_histogram(binwidth=0.1) +
#   xlim(0, 15)
#
# fin_females_filt %>%
#   ggplot(aes(1,dispGeneEst)) +
#   geom_boxplot()
