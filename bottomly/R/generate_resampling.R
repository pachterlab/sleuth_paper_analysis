library('dplyr')

# generate_samples <- function(df, grouping, n_group,
#   sort_column = 'sample') {
#   df <- dplyr::group_by_(df, grouping)
#   new_sample <- dplyr::sample_n(df, n_group)
#   new_sample <- dplyr::arrange_(new_sample, sort_column)
#   dplyr::ungroup(new_sample)
# }

metadata <- read.csv('../metadata/experiment.csv', stringsAsFactors = FALSE)

metadata <- dplyr::select(metadata,
  sample = run_accession, strain, library_name)

metadata <- dplyr::mutate(metadata,
  condition = factor(strain, levels = sort(unique(strain)),
    labels = c('A', 'B')))

###
# generate random permutations
###

set.seed(99)
n_samples <- 20

# training_sets <- lapply(1:n_samples,
#   function(discard) {
#     ret <- generate_samples(metadata, 'strain', 3)
#
#     dplyr::mutate(ret,
#       condition = factor(strain, levels = sort(unique(strain)),
#         labels = c('A', 'B')))
#   })
#
# ids <- lapply(training_sets,
#   function(x) {
#     paste(x$sample, collapse = '')
#   })

# # ensure that every single permutation is unique
# # if not, the user should manually come in and fix this
# # with the current seed and number of samples, it should be okay
# stopifnot(length(ids) == length(unique(ids)))
#
# validation_sets <- lapply(training_sets,
#   function(x) {
#     validation_sample <- setdiff(metadata$sample, x$sample)
#     validation_sample <- data.frame(sample = validation_sample,
#       stringsAsFactors = FALSE)
#     ret <- dplyr::inner_join(metadata, validation_sample, by = 'sample')
#     ret <- data.frame(ret, stringsAsFactors = FALSE)
#     rownames(ret) <- ret$sample
#
#     ret
#   })

# from the DESeq2paper package:
# http://www-huber.embl.de/DESeq2paper/DESeq2paper_1.3.tar.gz
random_subsets <- scan('../metadata/random_subsets.txt', what = character(),
  sep = '\n')

# now split the subsets
random_subsets <- lapply(random_subsets, function(x) strsplit(x, ' ')[[1]])

training <- lapply(random_subsets, function(x) x[1:6])
validation <- lapply(random_subsets, function(x) x[7:length(x)])

training <- training[1:20]
validation <- validation[1:20]

training_sets <- lapply(training,
  function(sample_names) {
    res <- dplyr::inner_join(metadata,
      data.frame(sample = sample_names, stringsAsFactors = FALSE),
      by = 'sample')
    res <- dplyr::arrange(res, condition)
    rownames(res) <- res$sample
    res
  })

validation_sets <- lapply(validation,
  function(sample_names) {
    res <- dplyr::inner_join(metadata,
      data.frame(sample = sample_names, stringsAsFactors = FALSE),
      by = 'sample')
    res <- dplyr::arrange(res, condition)
    rownames(res) <- res$sample
    res
  })


saveRDS(training_sets, file = '../metadata/training_sets.rds')
saveRDS(validation_sets, file = '../metadata/validation_sets.rds')

###
# dealing with plaintext so that we can read into cuffdiff2
###

samples <- list()

# the format in plain text is as follows:
# 1 file for each condition
# 1 line for a specific permutation
# each accession number separated by a space
# example:
# run1 run7 run3
# run7 run4 run2
# ...

permutations_to_string <- function(permutations, prefix) {
  ret <- list()
  ret$A <- sapply(permutations,
    function(df) {
      ret <- df$sample[df$condition == 'A']
      paste(ret, collapse = ' ')
    })

  ret$B <- sapply(permutations,
    function(df) {
      ret <- df$sample[df$condition == 'B']
      paste(ret, collapse = ' ')
    })

  cat(ret$A, file = paste0(prefix, '_a.txt'), sep = '\n')
  cat(ret$B, file = paste0(prefix, '_b.txt'), sep = '\n')

  ret
}

training_strings <- permutations_to_string(training_sets,
  '../metadata/training')
validation_strings <- permutations_to_string(validation_sets,
  '../metadata/validation')
