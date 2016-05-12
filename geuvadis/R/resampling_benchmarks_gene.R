args <- commandArgs(trailing = TRUE)

if (length(args) != 1) {
  stop('Must supply only 1 argument (the number of cores).')
}

cores <- 20

cores <- args[1]
sim_name <- 'gfr_10_11_1_43_2'

source('benchmark_methods.R')
source('gene_common.R')

library('dplyr')
library('mamabear')
library('parallel')

options(mc.cores = cores)

training_sets <- readRDS('../metadata/training_sets.rds')
validation_sets <- readRDS('../metadata/validation_sets.rds')

transcript_gene_mapping <- get_human_gene_names()

# # this is a temporary fix
# transcript_gene_mapping <- dplyr::mutate(transcript_gene_mapping,
#   target_id_backup = target_id)
# transcript_gene_mapping <- dplyr::mutate(transcript_gene_mapping,
#   target_id = ensembl_transcript_id)

###
# loading featureCounts data
###

training_sets <- lapply(training_sets,
  function(df) {
    sample_number <- sapply(strsplit(df$sample, '_'), '[[', 2)
    df$featureCounts <- file.path('..', 'sims', sim_name, 'exp_1', sample_number,
      'featureCounts.txt')
    df <- data.frame(df, stringsAsFactors = FALSE)
    rownames(df) <- df$sample
    df
  })

validation_sets <- lapply(validation_sets,
  function(df) {
    sample_number <- sapply(strsplit(df$sample, '_'), '[[', 2)
    df$featureCounts <- file.path('..', 'sims', sim_name, 'exp_1', sample_number,
      'featureCounts.txt')
    df <- data.frame(df, stringsAsFactors = FALSE)
    rownames(df) <- df$sample
    df
  })

# restrict to only the genes that are in the ensembl database
gene_names <- unique(transcript_gene_mapping$ens_gene)

dummy_filter <- rep(TRUE, length(gene_names))
names(dummy_filter) <- gene_names

training_counts <- lapply(training_sets,
  function(df) {
    obs <- load_union_counts_general(df$featureCounts, df$sample)
    current_filter <- intersect(rownames(obs), names(dummy_filter))
    obs[current_filter, ]
  })

validation_counts <- lapply(validation_sets,
  function(df) {
    obs <- load_union_counts_general(df$featureCounts, df$sample)
    current_filter <- intersect(rownames(obs), names(dummy_filter))
    obs[current_filter, ]
  })

# dummy_filter now contains the intersection between sleuth and featureCounts
dummy_filter <- rep(TRUE, nrow(validation_counts[[1]]))
names(dummy_filter) <- rownames(validation_counts[[1]])
dummy_filter_df <- data.frame(target_id = names(dummy_filter),
  stringsAsFactors = FALSE)

###
# run the gene methods
###

all_validation <- list()

sleuth_validation <- mclapply(seq_along(validation_sets),
  function(i) {
    validation <- validation_sets[[i]]
    sir <- run_sleuth(validation, gene_mode = 'aggregate', gene_column = 'ens_gene')
    list(sleuth.lrt = sir$sleuth.lrt, sleuth.wt = sir$sleuth.wt)
  })
all_validation$sleuth.lrt <- lapply(sleuth_validation, '[[', 'sleuth.lrt')
all_validation$sleuth.wt <- lapply(sleuth_validation, '[[', 'sleuth.wt')

all_validation$DESeq <- mclapply(seq_along(validation_sets),
  function(i) {
    validation <- validation_sets[[i]]
    obs <- validation_counts[[i]]
    DESeq_filter_and_run(obs, validation, dummy_filter)$result
  })

all_validation$DESeq2 <- mclapply(seq_along(validation_sets),
  function(i) {
    validation <- validation_sets[[i]]
    obs <- validation_counts[[i]]
    DESeq2_filter_and_run_intersect(obs, validation, dummy_filter)$result
  })

edgeR_results <- mclapply(seq_along(validation_sets),
  function(i) {
    validation <- validation_sets[[i]]
    obs <- validation_counts[[i]]
    edgeR_filter_and_run(obs, validation, dummy_filter)
  })
edgeR_filter_validation <- lapply(edgeR_results, '[[', 'filter')
all_validation$edgeR <- lapply(edgeR_results, '[[', 'result')

# use edgeR as the filter for limmaVoom and EBSeq
edgeR_filter_validation <- lapply(edgeR_filter_validation,
  function(x) {
    y <- dummy_filter
    y <- !y
    y[x] <- TRUE
    y
  })

all_validation$limmaVoom <- mclapply(seq_along(validation_sets),
  function(i) {
    validation <- validation_sets[[i]]
    obs <- validation_counts[[i]]
    current_filter <- edgeR_filter_validation[[i]]
    limma_filter_and_run(obs, validation, current_filter)$result
  })

all_validation$EBSeq <- mclapply(seq_along(validation_sets),
  function(i) {
    validation <- validation_sets[[i]]
    obs <- validation_counts[[i]]
    current_filter <- edgeR_filter_validation[[i]]
    EBSeq_gene_filter_and_run(obs, validation, current_filter)$result
  })

all_validation$Cuffdiff2 <- mclapply(seq_along(validation_sets),
  function(i) {
    path <- file.path('..', 'results', 'validation', sim_name, 'exp_1',
      i, 'results', 'cuffdiff')
    res <- get_cuffdiff(path)$gene
    res <- dplyr::filter(res, status == 'OK')
    # the next lines are relevant since cuffdiff will likely test on a superset
    # due to the larger gtf
    res <- dplyr::inner_join(res, dummy_filter_df, by = 'target_id')
    res <- dplyr::mutate(res, qval = p.adjust(pval, method = 'BH'))
    res
  })

all_validation <- lapply(all_validation,
  function(validation) {
    lapply(validation,
      function(x) {
        dplyr::mutate(x, is_de = ifelse(qval <= 0.05, TRUE, FALSE), log_fc = NA)
      })
  })

###
# run on the training set
###

all_training <- list()

sleuth_training <- mclapply(seq_along(training_sets),
  function(i) {
    training <- training_sets[[i]]
    sir <- run_sleuth(training, gene_mode = 'aggregate', gene_column = 'ens_gene')
    list(sleuth.lrt = sir$sleuth.lrt, sleuth.wt = sir$sleuth.wt)
  })
all_training$sleuth.lrt <- lapply(sleuth_training, '[[', 'sleuth.lrt')
all_training$sleuth.wt <- lapply(sleuth_training, '[[', 'sleuth.wt')

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

all_training$limmaVoom <- mclapply(seq_along(training_sets),
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

all_training$Cuffdiff2 <- mclapply(seq_along(training_sets),
  function(i) {
    path <- file.path('..', 'results', 'training', sim_name, 'exp_1',
      i, 'results', 'cuffdiff')
    res <- get_cuffdiff(path)$gene
    res <- dplyr::filter(res, status == 'OK')
    # the next lines are relevant since cuffdiff will likely test on a superset
    # due to the larger gtf
    res <- dplyr::inner_join(res, dummy_filter_df, by = 'target_id')
    res <- dplyr::mutate(res, qval = p.adjust(pval, method = 'BH'))
    res
  })

self_benchmark <- lapply(seq_along(all_training),
  function(i) {
    method <- names(all_training)[[i]]
    training <- all_training[[i]]
    validation <- all_validation[[method]]
    Map(
      function(x, y) new_de_benchmark(list(x), method, y),
      training, validation)
  })
names(self_benchmark) <- names(all_training)

saveRDS(self_benchmark, file.path('..', 'results', sim_name,
  'gene_self_benchmark.rds'))
