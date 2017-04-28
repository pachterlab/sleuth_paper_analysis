geuvadis_base <- '../finn_samples'

finn_females <- load('../results/fin_females.RData')

load('../results/fin_females.RData', verbose = TRUE)

load('../results/prep_fin.RData', verbose = TRUE)

load('../metadata/geu_meta.RData', verbose = TRUE)

sample_names <- data.frame(sample = colnames(fin_females),
  stringsAsFactors = FALSE)

finn_samples <- dplyr::inner_join(geu_meta, sample_names, by = 'sample')

finn_summary <- dplyr::group_by(finn_samples, lab)
finn_summary <- dplyr::summarize(finn_summary, n = n())
finn_summary <- dplyr::filter(finn_summary, n > 6)

finn_subset <- dplyr::inner_join(
  finn_samples,
  dplyr::select(finn_summary, lab),
  by = 'lab')

balanced_sample <- function(labels, n, n_samples) {
  proportions <- table(labels) / length(labels)
  indices <- 1:length(labels)
  lapply(1:n_samples,
    function(i) {
      which_group <- sample(names(proportions), 1, prob = proportions)
      sample(indices[labels == which_group], n, replace = FALSE)
    })
}

# test <- balanced_sample(finn_subset$lab, 6, 1000000)
# table(sapply(test, function(i) finn_subset[i, 'lab'][1])) / 1000000
# table(finn_subset$lab) / nrow(finn_subset)

set.seed(107)
finn_subsamples <- balanced_sample(finn_subset$lab, 6, 20)
finn_subsamples <- lapply(finn_subsamples, function(i) finn_subset[i, ])

# get the technical replicate name
all_geuvadis <- dir(geuvadis_base)

proper_names <- lapply(finn_subset$sample,
  function(x) {
    grep(x, all_geuvadis, value = TRUE)
  })

# ensure there are no weird duplicates (replicates)
tmp <- sapply(proper_names, length)
stopifnot(length(unique(tmp)) == 1)

write(unlist(proper_names), sep = '\n', file = '../finn_samples.txt')

# put the path in and randomly assign to groups
finn_subsamples <- lapply(finn_subsamples,
  function(df) {
    ret <- dplyr::mutate(df, path = file.path(geuvadis_base, true_name_mapping[sample],
      'abundance.h5'),
      condition = 'A')
    condition_b <- sample.int(6, 3, replace = FALSE)
    ret$condition[condition_b] <- 'B'
    ret
  })

finn_subsamples <- lapply(finn_subsamples,
  function(df) {
    df <- data.frame(df, stringsAsFactors = TRUE)
    rownames(df) <- df$sample
    df
  })

saveRDS(finn_subsamples, file = '../results/finn_subsamples.rds')

text_names <- sapply(finn_subsamples,
  function(x) {
    paste(c(true_name_mapping[x$sample[x$condition == 'A']],
      true_name_mapping[x$sample[x$condition == 'B']]),
      collapse = ' ', sep = '')
  })

write(text_names, sep = '\n', file = '../finn_null_experiment.txt')
