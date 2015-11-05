args <- commandArgs(trailingOnly = TRUE)

print(args)

if (length(args) != 2) {
  stop("Usage: Rscript merge_counts.R DIRECTORY OUTPUT.tsv")
}

suppressPackageStartupMessages(library('data.table'))
suppressPackageStartupMessages(library('dplyr'))

which_directory <- args[1]
output <- args[2]

directories <- Sys.glob(file.path(which_directory, '*', 'featureCounts.txt'))
samples <- sub(file.path(which_directory, ''), '', directories)
samples <- sub(file.path('', 'featureCounts.txt'), '', samples)

readFeatureCounts <- function(path) {
  result <- suppressWarnings( data.table::fread(path, data.table = FALSE) )
  result <- result[, c(1, 7)]
  setnames(result, c('gene_id', 'count'))

  result
}

all_counts <- lapply(seq_along(directories),
  function(i) {
    cur_directory <- directories[i]
    cur_sample <- samples[i]
    counts <- readFeatureCounts(cur_directory)
    setnames(counts, c('count'), cur_sample)

    counts
    })

counts_table <- Reduce(function(x, y) inner_join(x, y, by = 'gene_id'),
  all_counts)

apply(counts_table[, 2:7], 2, sum)

write.table(counts_table, output, quote = FALSE, sep = "\t", row.names = FALSE)
