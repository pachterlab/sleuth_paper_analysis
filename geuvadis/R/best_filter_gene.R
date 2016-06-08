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

valid_genes <- unique(transcript_gene_mapping$ens_gene)
valid_genes <- valid_genes[!is.na(valid_genes)]

all_counts <- lapply(1:sim$n,
  function(i) {
    load_union_counts(sim, i)
  })
dummy_count_filter <- rep(FALSE, nrow(all_counts[[1]]))
names(dummy_count_filter) <- rownames(all_counts[[1]])
dummy_count_filter[intersect(valid_genes, names(dummy_count_filter))] <- TRUE
dummy_count_filter <- dummy_count_filter[rownames(all_counts[[1]])]

each_filter <- list()


each_filter$sleuth <- mclapply(1:sim$n,
# each_filter$sleuth <- lapply(1:sim$n,
  function(i) {
    si <- sample_info[[i]]
    # debugonce(sleuth_prep)
    res <- run_sleuth(si, gene_mode = 'aggregate', gene_column = 'ens_gene')
    res$so <- NULL
    res
  })
#
# i <- 1
# si <- sample_info[[i]]
# debugonce(sleuth_prep)
# res <- run_sleuth(si, gene_mode = 'aggregate', gene_column = 'ens_gene')

each_filter$DESeq <- mclapply(1:sim$n,
  function(i) {
    si <- sample_info[[i]]
    counts <- all_counts[[i]]
    DESeq_filter_and_run(counts, si, dummy_count_filter)$result
  })

each_filter$DESeq2 <- mclapply(1:sim$n,
  function(i) {
    si <- sample_info[[i]]
    counts <- all_counts[[i]]
    DESeq2_filter_and_run_intersect(counts, si, dummy_count_filter)$result
  })

each_filter$edgeR <- mclapply(1:sim$n,
  function(i) {
    si <- sample_info[[i]]
    counts <- all_counts[[i]]
    edgeR_filter_and_run(counts, si, dummy_count_filter)$result
  })

each_filter$limmaVoom <- mclapply(1:sim$n,
  function(i) {
    si <- sample_info[[i]]
    counts <- all_counts[[i]]
    current_filter <- edgeR_filter(counts)
    limma_filter_and_run(counts, si, current_filter)$result
  })

each_filter$EBSeq <- mclapply(1:sim$n,
  function(i) {
    si <- sample_info[[i]]
    counts <- all_counts[[i]]
    current_filter <- edgeR_filter(counts)
    EBSeq_gene_filter_and_run(counts, si, current_filter)$result
  })

each_filter$Cuffdiff2 <- mclapply(1:sim$n,
  function(i) {
    si <- sample_info[[i]]
    cuffdiff_filter_and_run_gene(NULL, si, dummy_count_filter)$result
  })

each_filter$sleuth <- lapply(each_filter$sleuth, '[[', 'sleuth.lrt')

# each_filter_benchmark <- lapply(1:20,
each_filter_benchmark <- mclapply(1:sim$n,
  function(i) {
    current <- lapply(each_filter, '[[', i)
    current_names <- names(each_filter)
    sim_info <- get_de_info(sim_name, i, transcript_gene_mapping)
    new_de_benchmark(current, current_names, sim_info$de_genes,
      join_mode = 'union')
  })

saveRDS(each_filter_benchmark, file = paste0('../results/', sim_name,
  '/gene_benchmarks_filter.rds'))
