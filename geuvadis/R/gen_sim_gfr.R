args <- commandArgs(trailingOnly = TRUE)

if( length(args) != 6) {
  stop("Usage: gen_sim_gfr.R BASE_OUT N_SIMULATIONS N_A N_B SEED SIZE_FACTOR_TYPE")
}

cat("Args: ", args, "\n")

# temporary
# sim_name <- 'tmp'
# n_sim <- 2L
# n_a <- 3L
# n_b <- 3L
# seed <- 1
# sf_type <- 2
# log_fc <- 1
# end temporary

sim_name <- args[1]
n_sim <- as.integer(args[2])
n_a <- as.integer(args[3])
n_b <- as.integer(args[4])
seed <- as.integer(args[5])
sf_type <- as.integer(args[6])
# log_fc <- as.integer(args[7])

n_total <- n_a + n_b

sim_base <- '../sims'
sim_out <- file.path(sim_base, sim_name)

source("../../simulation_core/R/simulate_de.R")
source("../../simulation_core/R/sim_to_rsem.R")

set.seed(seed)
sf_function <- switch(sf_type,
  sf_mode_1,
  sf_mode_2
  )

load("../results/prep_fin.RData", verbose = TRUE)

cat("Generating counts\n")

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = "may2015.archive.ensembl.org")
  # host="www.ensembl.org")
ttg <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)

# rename the columns: ensure we have 'target_id' and also make the others a bit
# shorter
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

fc <- readRDS('../../cuffdiff2_analysis/results/fc.rds')

gfcr <- function(df) {
  # remove the things that are invalid
  possible <- unique(fc$n_transcripts)
  invalid <- setdiff(unique(df$n_transcripts), possible)
  invalid_rows <- df$n_transcripts %in% invalid
  df$is_de[invalid_rows] <- FALSE
  warning('Invalidated ', sum(invalid_rows),
    ' because no match in the reference')

  original_order <- df$target_id

  res <- gene_fold_change_reference(df, fc)

  res[match(original_order, res$target_id), ]
}

# currently a dummy
log_fc <- 1

sim <- simulate_gene_counts(prep_fin,
  target_mapping = ttg,
  gene_label = "ens_gene",
  n_sim = n_sim,
  n_a = n_a,
  n_b = n_b,
  prop_de = 0.20,
  seed = seed,
  # log_fc = log_fc,
  generate_fold_changes = gfcr,
  # log_fc_sd = 1,
  size_factor_function = sf_function)

rsem_fhandle <- gzfile("../results/rsem/HG00365_7/out.isoforms.results.gz", open = "r")
rsem_res <- read.table(rsem_fhandle, header = TRUE, stringsAsFactors = FALSE)

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
