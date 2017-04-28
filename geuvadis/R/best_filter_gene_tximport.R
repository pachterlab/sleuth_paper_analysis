args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop('Usage: Rscript bench_gene.R N_CPU SIM_NAME')
}

# temporary for debugging
n_cpu <- 20
sim_name <- 'gfr_3_3_20_42_2'

n_cpu <- args[1]
sim_name <- args[2]


source("benchmark_methods.R")
source("gene_common.R")

library("parallel")
library("mamabear")
library("tximport")
options(mc.cores = n_cpu)

transcript_gene_mapping <- get_human_gene_names()

###
# run each method with their respective filters
###

sim <- parse_simulation(sim_name)

sample_info <- lapply(1:20,
  function(i) {
    n <- sim$a + sim$b
    kal_dirs <- file.path('..', 'sims', sim_name, paste0("exp_", i),
      1:n, "kallisto")
    get_sample_to_condition(sim$a, sim$b, kal_dirs)
  })

tx_files <- lapply(sample_info,
  function(x) {
    new_path <- file.path(x$path, 'abundance.tsv')
    names(new_path) <- x$sample
    new_path
  })

tx_mapping <- transcript_gene_mapping
tx_mapping <- dplyr::select(tx_mapping, TXNAME = target_id, GENEID = ens_gene)

tx_counts <- mclapply(tx_files,
  function(filenames) {
    tximport(filenames, type = 'kallisto', tx2gene = tx_mapping)
  })

each_filter_txi <- list()

dummy_filter <- TRUE

each_filter_txi$DESeq2 <- mclapply(1:sim$n,
  function(i) {
    si <- sample_info[[i]]
    txi <- tx_counts[[i]]
    DESeq2_filter_and_run_intersect(txi, si, dummy_filter, FALSE)$result
  })

each_filter_txi$edgeR <- mclapply(1:sim$n,
  function(i) {
    si <- sample_info[[i]]
    txi <- tx_counts[[i]]
    edgeR_filter_and_run(txi, si, dummy_filter, FALSE)$result
  })

# only used with limma
tx_scaled_counts <- mclapply(tx_files,
  function(filenames) {
    tximport(filenames, type = 'kallisto', tx2gene = tx_mapping,
      countsFromAbundance = "lengthScaledTPM")
  })

each_filter_txi$voom <- mclapply(1:sim$n,
  function(i) {
    si <- sample_info[[i]]
    txi <- tx_scaled_counts[[i]]
    current_filter <- edgeR_filter(txi$counts)
    limma_filter_and_run(txi$counts, si, current_filter)$result
  })

saveRDS(each_filter_txi, file = paste0('../results/', sim_name, '/txi.rds'))
