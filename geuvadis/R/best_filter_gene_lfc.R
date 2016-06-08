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

lfc_filter <- list()

set.seed(101)
lfc_filter$LFC <- lapply(1:sim$n,
  function(i) {
    si <- sample_info[[i]]
    counts <- all_counts[[i]]
    res <- lfc_filter_and_run (counts, si, dummy_count_filter)$result
  })

set.seed(101)
lfc_filter$GLFC <- lapply(1:sim$n,
  function(i) {
    si <- sample_info[[i]]
    counts <- all_counts[[i]]
    res <- geometric_lfc_filter_and_run (counts, si, dummy_count_filter)$result
  })

###
# end temporary log fold change
###

each_filter_benchmark <- readRDS(paste0('../results/', sim_name,
  '/gene_benchmarks_filter.rds'))

each_filter <- lapply(each_filter_benchmark[[1]]$labels,
  function(label) {
    lapply(each_filter_benchmark, function(bench) {
      bench$original_data[[label]]
    })
  })
names(each_filter) <- each_filter_benchmark[[1]]$labels

each_filter_lfc <- c(each_filter, lfc_filter)
each_filter_benchmark_lfc <- lapply(1:sim$n,
  function(i) {
    current <- lapply(each_filter_lfc, '[[', i)
    current_names <- names(each_filter_lfc)
    sim_info <- get_de_info(sim_name, i, transcript_gene_mapping)
    new_de_benchmark(current, current_names, sim_info$de_genes,
      join_mode = 'union')
  })

saveRDS(each_filter_benchmark_lfc, file = paste0('../results/', sim_name,
  '/gene_benchmarks_filter_lfc.rds'))


# temporary place for colors
