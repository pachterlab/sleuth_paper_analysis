# Merge counts from `featureCounts` output

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript merge_counts.R DIRECTORY OUTPUT.tsv")
}

suppressPackageStartupMessages(library('data.table'))
suppressPackageStartupMessages(library('dplyr'))

which_directory <- args[1]
output <- args[2]

cat('The base directory: ', which_directory, '\n')
cat('Counts output: ', output, '\n')

# this code assumes the following structure:
# which_directory/1/featureCounts.txt
# which_directory/2/featureCounts.txt
# ...
# which_directory/N/featureCounts.txt
directories <- Sys.glob(file.path(which_directory, '*', 'featureCounts.txt'))
samples <- sub(file.path(which_directory, ''), '', directories)
samples <- sub(file.path('', 'featureCounts.txt'), '', samples)
# to avoid issues due to lexicographic sorting
samples <- sort(as.integer(samples))
directories <- file.path(which_directory, samples, 'featureCounts.txt')

print(samples)

cat('The number of samples: ', length(samples), '\n')
cat('The featureCounts results to be read in:\n')
cat(paste('\n\t', directories))
cat('\n')

readFeatureCounts <- function(path) {
  result <- suppressWarnings( data.table::fread(path, data.table = FALSE) )

  # pull out the columns 'gene_id' and 'count'
  result <- result[, c(1, 7)]
  data.table::setnames(result, c('gene_id', 'count'))

  result
}

all_counts <- lapply(seq_along(directories),
  function(i) {
    cur_directory <- directories[i]
    cur_sample <- samples[i]
    counts <- readFeatureCounts(cur_directory)
    data.table::setnames(counts, c('count'), as.character(cur_sample))

    counts
  })

counts_table <- Reduce(
  function(x, y) inner_join(x, y, by = 'gene_id'),
  all_counts
  )

# a sanity check to see the total number of reads mapped
apply(counts_table[, 2:7], 2, sum)

write.table(counts_table, output, quote = FALSE, sep = "\t", row.names = FALSE)
