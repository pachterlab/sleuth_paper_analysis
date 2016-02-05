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
# n_sim <- 1
# n_a <- 3
# n_b <- 3
# seed <- 42
# sf_type <- 1

n_total <- n_a + n_b

sf_mode_1 <- function(n) {
  rep.int(1, n)
}

sf_mode_2 <- function(n) {
  nsamp <- 0
  res <- numeric(n)

  while (nsamp < n) {
    if ( (n - nsamp) == 1 ) {
      nsamp <- nsamp + 1
      res[nsamp] <- 1
    } else {
      eq <- sample(c(TRUE, FALSE), 1)
      if (eq) {
        nsamp <- nsamp + 1
        res[nsamp] <- 1
      } else {
        nsamp <- nsamp + 1
        res[nsamp] <- 1 / 3
        nsamp <- nsamp + 1
        res[nsamp] <- 3
      }
    }
  }

  res[sample.int(length(res))]
}

set.seed(seed)
sfs <- switch(sf_type,
  sf_mode_1(n_total),
  sf_mode_2(n_total)
  )


sim_base <- '../sims'
sim_out <- file.path(sim_base, sim_name)

source("../../simulation_core/R/simulate_de.R")
source("../../simulation_core/R/sim_to_rsem.R")

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
  size_factors = sfs,
  constant_fc = NULL
  )
message(paste0('Size factors: ', sfs))

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

saveRDS(sim, file.path(sim_out, 'sim.rds'))

invisible(lapply(seq_along(sim$counts),
  function(i) {
    cat("Writing simulation ", i, "\n")
    invisible(counts_to_simulation(sim$counts[[i]], sim$info, rsem_res,
      file.path(sim_out, paste0("exp_", i))))
    NULL
  }))
