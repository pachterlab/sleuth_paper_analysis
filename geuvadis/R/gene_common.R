
get_de_info <- function(sim_name, ttg) {
  sims_dir <- file.path('../sims', sim_name)
  base_dir <- file.path('../results', sim_name)

  # de_info contains the 'true' fold change as well as which transcripts are DE
  de_info <- read.table(
    gzfile(file.path(sims_dir, "de_info.tsv.gz")),
    header = TRUE,
    stringsAsFactors = FALSE
    )

  tmp <- dplyr::inner_join(de_info, ttg, by = 'target_id')

  de_genes <- dplyr::select(tmp, ens_gene, is_de)
  de_genes <- dplyr::group_by(de_genes, ens_gene)
  de_genes <- dplyr::summarize(de_genes, is_de = any(is_de))

  de_genes <- dplyr::rename(de_genes, target_id = ens_gene)
  de_genes <- dplyr::mutate(de_genes, log_fc = NA)

  list(de_info = de_info, de_genes = de_genes)
}

get_gene_counts <- function(sim_name, sim_num = 1) {
  base_dir <- file.path('../results', sim_name)

  experiment_name <- paste0('exp_', sim_num)
  counts <- read.table(file.path(base_dir, experiment_name, 'counts.tsv'),
    header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  counts <- as.matrix(counts)
  colnames(counts) <- sub('X', 'sample_', colnames(counts))

  counts
}
