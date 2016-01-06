args <- commandArgs(trailingOnly = TRUE)

if( length(args) != 7) {
  stop("Usage: gen_sim.R BASE_OUT N_SIMULATIONS N_A N_B SEED SIZE_FACTOR_TYPE")
}

cat("Args: ", args, "\n")

sim_name <- args[1]
n_sim <- as.integer(args[2])
n_a <- as.integer(args[3])
n_b <- as.integer(args[4])
seed <- as.integer(args[5])
sf_type <- as.integer(args[6])
log_fc <- as.integer(args[7])

# temporary
sim_name <- 'tmp'
n_sim <- 1L
n_a <- 3L
n_b <- 3L
seed <- 1
sf_type <- 2
log_fc <- 1
# end temporary

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
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host="www.ensembl.org")
ttg <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)
# ttg <- readRDS('~/ttg.rds')

# rename the columns: ensure we have 'target_id' and also make the others a bit
# shorter
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# debugonce(simulate_gene_counts)
sim <- simulate_gene_counts(prep_fin,
  target_mapping = ttg,
  gene_label = "ens_gene",
  n_sim = n_sim,
  n_a = n_a,
  n_b = n_b,
  prop_de = 0.20,
  seed = seed,
  log_fc = log_fc,
  # log_fc_sd = 1,
  size_factors = sfs)
message(paste0('Size factors: ', sfs))

# # verify that `gene_fold_change` correctly assigns the fold change
# filtered <- sim$info %>%
#   filter(ens_gene == 'ENSG00000000419' | ens_gene == 'ENSG00000006210')
#
# undebug(gene_fold_change_same_direction)
# arrange(gene_fold_change(filtered), ens_gene)
#
# # ensure simulation is valid
# filtered <- sim$info %>%
#   filter(is_de) %>%
#   group_by(ens_gene) %>%
#   summarize(all_same_sign = all(sign(log_fc)))
# all(filtered$all_same_sign)


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
