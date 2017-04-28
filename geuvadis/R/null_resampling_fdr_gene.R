args <- commandArgs(trailing = TRUE)

# geuvadis_base <- '../results/finn_samples'

if (length(args) != 1) {
  stop('Must supply only 1 argument (the number of cores).')
}

cores <- 20

cores <- args[1]

source('benchmark_methods.R')
source('gene_common.R')

library('dplyr')
library('mamabear')
library('parallel')

options(mc.cores = cores)

transcript_gene_mapping <- get_human_gene_names()

sample_info <- readRDS('../results/finn_subsamples.rds')

gene_names <- unique(transcript_gene_mapping$ens_gene)
gene_names <- gene_names[!is.na(gene_names)]
# restrict to only the genes that are in the ensembl database
gene_names <- unique(transcript_gene_mapping$ens_gene)

dummy_filter <- rep(TRUE, length(gene_names))
names(dummy_filter) <- gene_names


###
# loading featureCounts data
###

training_sets <- lapply(sample_info,
  function(df) {
    sample_path <- sapply(strsplit(df$path, '/'), '[[', 4)
    df$featureCounts <- file.path('..', 'results', 'finn_samples', sample_path,
      'featureCounts.txt')
    df
  })

training_counts <- lapply(training_sets,
  function(df) {
    obs <- load_union_counts_general(df$featureCounts, df$sample)
    current_filter <- intersect(rownames(obs), names(dummy_filter))
    obs[current_filter, ]
  })

dummy_filter <- rownames(training_counts[[1]])
names(dummy_filter) <- dummy_filter
dummy_filter[1:length(dummy_filter)] <- TRUE
storage.mode(dummy_filter) <- 'logical'

###
# run the gene methods
###

all_training <- list()

sleuth_training <- mclapply(seq_along(training_sets),
  function(i) {
    training <- training_sets[[i]]
    sir <- run_sleuth(training, gene_mode = 'aggregate', gene_column = 'ens_gene')
    list(sleuth.lrt = sir$sleuth.lrt, sleuth.wt = sir$sleuth.wt)
  })
all_training$sleuth <- lapply(sleuth_training, '[[', 'sleuth.lrt')

all_training$DESeq <- mclapply(seq_along(training_sets),
  function(i) {
    training <- training_sets[[i]]
    obs <- training_counts[[i]]
    DESeq_filter_and_run(obs, training, dummy_filter)$result
  })

all_training$DESeq2 <- mclapply(seq_along(training_sets),
  function(i) {
    training <- training_sets[[i]]
    obs <- training_counts[[i]]
    DESeq2_filter_and_run_intersect(obs, training, dummy_filter)$result
  })

edgeR_training <- mclapply(seq_along(training_sets),
  function(i) {
    training <- training_sets[[i]]
    obs <- training_counts[[i]]
    edgeR_filter_and_run(obs, training, dummy_filter)
  })
edgeR_filter_training <- lapply(edgeR_training, '[[', 'filter')
all_training$edgeR <- lapply(edgeR_training, '[[', 'result')

# use edgeR as the filter for limmaVoom and EBSeq
edgeR_filter_training <- lapply(edgeR_filter_training,
  function(x) {
    y <- dummy_filter
    y <- !y
    y[x] <- TRUE
    y
  })

all_training$voom <- mclapply(seq_along(training_sets),
  function(i) {
    training <- training_sets[[i]]
    obs <- training_counts[[i]]
    current_filter <- edgeR_filter_training[[i]]
    limma_filter_and_run(obs, training, current_filter)$result
  })

all_training$EBSeq <- mclapply(seq_along(training_sets),
  function(i) {
    training <- training_sets[[i]]
    obs <- training_counts[[i]]
    current_filter <- edgeR_filter_training[[i]]
    EBSeq_gene_filter_and_run(obs, training, current_filter)$result
  })

oracle <- data.frame(target_id = names(dummy_filter), is_de = FALSE,
  log_fc = NA, stringsAsFactors = FALSE)
oracle <- lapply(1:length(sleuth_training), function(i) oracle)


self_benchmark <- lapply(seq_along(all_training),
  function(i) {
    method <- names(all_training)[[i]]
    training <- all_training[[i]]
    Map(
      function(x, y) new_de_benchmark(list(x), method, y),
      training, oracle)
    })

saveRDS(self_benchmark, file = '../results/null_resampling/gene.rds')
