source('gene_common.R')
source('benchmark_methods.R')

ttg <- get_human_gene_names()

sim_name <- 'gfr_10_11_1_43_2'

###
# get the full design
###
which_sample <- 1
# which_sample <- as.integer(which_sample)
sim <- parse_simulation(sim_name)
n <- sim$a + sim$b
kal_dirs <- file.path('..', 'sims', sim_name, paste0("exp_", which_sample),
  1:n, "kallisto")
metadata <- get_sample_to_condition(sim$a, sim$b, kal_dirs)

set.seed(100)
n_samples <- 20

training_sets <- lapply(1:n_samples,
  function(discard) {
    generate_samples(sample_to_condition, 'condition', 3, sort_column = 'sample')
  })

ids <- lapply(training_sets,
  function(x) {
    paste(x$sample, collapse = '')
  })

# ensure that every single permutation is unique
# if not, the user should manually come in and fix this
# with the current seed and number of samples, it should be okay
stopifnot(length(ids) == length(unique(ids)))

validation_sets <- lapply(training_sets,
  function(x) {
    validation_sample <- setdiff(metadata$sample, x$sample)
    validation_sample <- data.frame(sample = validation_sample,
      stringsAsFactors = FALSE)
    ret <- dplyr::inner_join(metadata, validation_sample, by = 'sample')
    ret <- data.frame(ret, stringsAsFactors = FALSE)
    rownames(ret) <- ret$sample

    ret
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
      gsub('sample_', '', paste(ret, collapse = ' '))
    })

  ret$B <- sapply(permutations,
    function(df) {
      ret <- df$sample[df$condition == 'B']
      gsub('sample_', '', paste(ret, collapse = ' '))
    })

  cat(ret$A, file = paste0(prefix, '_a.txt'), sep = '\n')
  cat(ret$B, file = paste0(prefix, '_b.txt'), sep = '\n')

  ret
}

training_strings <- permutations_to_string(training_sets,
  '../metadata/training')
validation_strings <- permutations_to_string(validation_sets,
  '../metadata/validation')
