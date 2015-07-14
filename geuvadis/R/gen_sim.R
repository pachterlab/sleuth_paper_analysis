args <- commandArgs(trailingOnly = TRUE)

if( length(args) != 3) {
  stop("Usage: gen_sim.R N_SIMULATIONS N_A N_B")
}

cat("Args: ", args, "\n")

n_sim <- args[1]
n_a <- args[2]
n_b <- args[3]

source("../../simulation_core/R/simulate_de.R")
source("../../simulation_core/R/sim_to_rsem.R")

load("../results/prep_fin.RData", verbose = TRUE)

cat("Generating counts\n")

sim <- simulate_counts(prep_fin,
  n_sim = n_sim,
  n_a = n_a,
  n_b = n_b,
  prop_de = 0.20,
  seed = 37L,
  log_fc_sd = 1)

rsem_fhandle <- gzfile("../results/rsem/HG00365_7/out.isoforms.results.gz")
rsem_res <- read.table(rsem_fhandle, header = TRUE, stringsAsFactors = FALSE)
#close(rsem_fhandle)

debugonce(counts_to_simulation)
ah <- counts_to_simulation(sim$counts[[1]], sim$info, rsem_res, "tmp")
