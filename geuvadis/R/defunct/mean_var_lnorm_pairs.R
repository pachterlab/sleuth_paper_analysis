# load DF prep_fin which is ready to simulate from
load("../results/fin_females.RData", verbose = TRUE)

# load the negative binomial analog to this so we can keep the same filters
load("../results/prep_fin.RData", verbose = TRUE)

size_factors <- DESeq2::estimateSizeFactorsForMatrix(fin_females)
fin_females_norm <- t(t(fin_females) / size_factors)

prep_fin_lnorm <- prep_fin

not_all_zero <- prep_fin$target_id[!prep_fin$allZero]

# small epsilon to avoid -Inf in calculations for low counts
eps <- 0.01

cat("Computing log-normal MLEs\n")

means <- apply(
  fin_females_norm[not_all_zero,],
  1,
  function(x) mean(log(x + eps)))

vars <- apply(
  fin_females_norm[not_all_zero,],
  1,
  function(x) var(log(x + eps)))
n <- ncol(fin_females_norm)
vars <- vars * (n - 1) / (n)

rownames(prep_fin_lnorm) <- prep_fin_lnorm$target_id
prep_fin_lnorm[not_all_zero,'baseMean'] <- means
prep_fin_lnorm[not_all_zero,'baseVar'] <- vars

prep_fin_lnorm[prep_fin_lnorm$allZero,'baseMean'] <- -Inf

prep_fin_lnorm <- dplyr::select(prep_fin_lnorm,
  baseMean, baseVar, allZero, target_id, disp_filt, sim_filt)
prep_fin_lnorm <- as.data.frame(prep_fin_lnorm)
rownames(prep_fin_lnorm) <- rownames(prep_fin)

cat("saving data.frame\n")
save(prep_fin_lnorm, file = "../results/prep_fin_lnorm.RData")
