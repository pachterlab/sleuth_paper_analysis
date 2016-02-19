args <- commandArgs(trailingOnly = TRUE)

if( length(args) != 5) {
  stop("Usage: gen_sim.R BASE_OUT N_SIMULATIONS N_A N_B")
}

cat("Args: ", args, "\n")

sim_name <- args[1]
n_sim <- as.integer(args[2])
n_a <- as.integer(args[3])
n_b <- as.integer(args[4])
seed <- as.integer(args[5])

sim_base <- '../sims'
sim_out <- file.path(sim_base, sim_name)

source("../../simulation_core/R/simulate_de.R")
source("../../simulation_core/R/sim_to_rsem.R")

load("../results/prep_fin.RData", verbose = TRUE)

make_sim_norm <- function(sim_df, X, size_factors = c(1, nrow(X))) {
  n <- nrow(sim_df)
  p <- nrow(X)

  beta_full <- cbind(log(sim_df$baseMean), sim_df$log_fc)
  mu <- exp(beta_full %*% t(X))

  counts <- matrix(rnorm(n * p, mean = mu, sd = sqrt(sim_df$baseVar)), ncol = p)
  counts <- round(counts)
  counts[sim_df$allZero,] <- 0
  counts[counts < 0] <- 0

  counts
}


cat("Generating counts\n")

sim <- simulate_counts(prep_fin,
  n_sim = n_sim,
  n_a = n_a,
  n_b = n_b,
  prop_de = 0.20,
  seed = seed,
  log_fc_sd = 1,
  sim_function = make_sim_norm)

# debugonce(make_sim_lnorm)
# sim <- simulate_counts(prep_fin_lnorm,
#   n_sim = n_sim,
#   n_a = n_a,
#   n_b = n_b,
#   prop_de = 0.20,
#   seed = seed,
#   log_fc_sd = 1,
#   sim_function = make_sim_lnorm)
#
# sim_means <- rowMeans(sim$counts[[1]])
# emp_means <- prep_fin$baseMean
# w_gt <- which(sim_means > 7 * emp_means)
#
# plot(sim_means[w_gt], emp_means[w_gt])
#
# fin_gt <- prep_fin[w_gt,] %>% arrange(desc(baseVar))
#
# i <- 1
# tid <- fin_gt$target_id[i]
#
# norm_glm <- glm(y ~ 1,
#   family = gaussian(),
#   data = data.frame(y = log(fin_females_norm[tid,] + eps)))
#
# gam_glm <- glm(y ~ 1,
#   family = Gamma(link = "log"),
#   data = data.frame(y = log(fin_females_norm[tid,] + eps)))
#
# hist(fin_females_norm[tid,], breaks = ncol(fin_females_norm) / 2,
#   freq = TRUE)
# hist(log(fin_females_norm[tid,]), breaks = ncol(fin_females_norm) / 2,
#   freq = TRUE)
#

rsem_fhandle <- gzfile("../results/rsem/HG00365_7/out.isoforms.results.gz", open = "r")
rsem_res <- read.table(rsem_fhandle, header = TRUE, stringsAsFactors = FALSE)

dir.create(sim_out, showWarnings = FALSE, recursive = TRUE)
cat("sim_out: ", sim_out, "\n")

cat("Writing out DE information\n")
de_handle <- gzfile(file.path(sim_out, "de_info.tsv.gz"), open = "w")
write.table(sim$info,
  file = de_handle,
  sep = '\t',
  row.names = FALSE,
  quote = FALSE,
  col.names = TRUE)
close(de_handle)

invisible(lapply(seq_along(sim$counts),
  function(i) {
    cat("Writing simulation ", i, "\n")
    invisible(counts_to_simulation(sim$counts[[i]], sim$info, rsem_res,
      file.path(sim_out, paste0("sample_", i))))
    NULL
  }))
