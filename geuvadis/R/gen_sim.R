args <- commandArgs(trailingOnly = TRUE)

if( length(args) != 6) {
  stop("Usage: gen_sim.R BASE_OUT N_SIMULATIONS N_A N_B SEED SIZE_FACTOR_TYPE")
}

cat("Args: ", args, "\n")

sim_name <- args[1]
n_sim <- as.integer(args[2])
n_a <- as.integer(args[3])
n_b <- as.integer(args[4])
seed <- as.integer(args[5])
sf_type <- as.integer(args[6])

# tests go here
# sim_name <- 'mySimulation'
# n_sim <- 6
# n_a <- 3
# n_b <- 3
# seed <- 42
# sf_type <- 1

n_total <- n_a + n_b

set.seed(seed)

sim_base <- '../sims'
sim_out <- file.path(sim_base, sim_name)

source("../../simulation_core/R/simulate_de.R")
source("../../simulation_core/R/sim_to_rsem.R")

sf_mode <- switch(sf_type,
  sf_mode_1,
  sf_mode_2
  )

load("../results/prep_fin.RData", verbose = TRUE)

cat("Generating counts\n")
# debugonce(simulate_counts)
#
# debugonce(make_sim)

sim <- simulate_counts(prep_fin,
  n_sim = n_sim,
  n_a = n_a,
  n_b = n_b,
  prop_de = 0.20,
  seed = seed,
  log_fc_sd = 1,
  size_factor_function = sf_mode,
  constant_fc = NULL
  )

# get the RSEM 'template'
rsem_fhandle <- gzfile("../results/rsem/HG00365_7/out.isoforms.results.gz",
  open = "r")
rsem_res <- read.table(rsem_fhandle, header = TRUE, stringsAsFactors = FALSE)
close(rsem_fhandle)

dir.create(sim_out, showWarnings = FALSE, recursive = TRUE)
cat("sim_out: ", sim_out, "\n")

saveRDS(sim, file.path(sim_out, 'sims.rds'))

maximum_seed <- as.integer(2 ^ 16) - 2
seeds <- sample.int(maximum_seed, length(sim))

invisible(lapply(seq_along(sim),
  function(i) {
    cat("Writing simulation ", i, "\n")
    cur_dir <- file.path(sim_out, paste0('exp_', i))
    invisible(counts_to_simulation(
      sim[[i]]$counts,
      sim[[i]]$info,
      rsem_res,
      seeds[i],
      cur_dir
      ))
    NULL
  }))
