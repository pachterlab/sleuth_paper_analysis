source('load_info.R')
source('../../geuvadis/R/benchmark_methods.R')
source('../../simulation_core/R/simulate_de.R')

library('DESeq2')

file_names <- file.path(dirname(info$path), 'featureCounts.txt')

counts_list <- lapply(seq_along(file_names),
  function(i) {
    f <- file_names[i]
    # df <- data.table::fread(f, skip = 1, data.table = FALSE)
    df <- read.table(f, sep = '\t', skip = 1, header = TRUE)
    df <- dplyr::select(df, 1, 7)
    data.table::setnames(df, c('ens_gene', info$sample[i]))
    df
  })

counts <- Reduce(function(x, y) dplyr::left_join(x, y, by = 'ens_gene'),
  counts_list)
rownames(counts) <- counts$ens_gene
counts$ens_gene <- NULL
counts <- as.matrix(counts)

pass_filter <- DESeq2_filter(counts)
counts <- counts[pass_filter, ]

rownames(info) <- info$sample
cds <- make_count_data_set(counts, info)

res <- runDESeq2(cds)
saveRDS(res, '../results/DESeq2.rds')

counts <- Reduce(function(x, y) dplyr::left_join(x, y, by = 'ens_gene'),
  counts_list)
rownames(counts) <- counts$ens_gene
counts$ens_gene <- NULL
counts <- as.matrix(counts)

all.equal(info$sample, colnames(counts))

dr <- res
dr <- dplyr::rename(dr, target_id = ens_gene)
lr <- limma_filter_and_run(counts, info, TRUE)
lr <- lr$result

dl_significant <- dplyr::inner_join(
  dplyr::rename(dplyr::filter(dr, qval < 0.10),
    qval_DESeq2 = qval, pval_DESeq2 = pval),
  dplyr::rename(dplyr::filter(lr, qval < 0.10),
    qval_voom = qval, pval_voom = pval),
  by = 'target_id')
saveRDS(dl_significant, '../results/dl_significant.rds')

normalized_counts <- t(t(counts / DESeq2::estimateSizeFactorsForMatrix(counts)))
