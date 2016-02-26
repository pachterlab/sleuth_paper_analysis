#' Load the differential expression truth
#'
#' This function loads the isoform level information and gene level information.
#' The gene level differential expression is defined by at least 1 isoform of the gene being differentially expressed.
#'
#' @param sim_name the name of the simulation (with respect to the results location)
#' @param ttg the 'target_to_gene_mapping' coming from ensembl
#' @return a list containing the gene and isoform differential expression information
get_de_info <- function(sim_name, sim_number, ttg) {
  sims_dir <- file.path('../sims', sim_name)

  # de_info contains the 'true' fold change as well as which transcripts are DE
  all_sim <- readRDS(file.path(sims_dir, 'sims.rds'))
  de_info <- all_sim[[sim_number]]$info

  tmp <- dplyr::inner_join(de_info, ttg, by = 'target_id')

  de_genes <- dplyr::select(tmp, ens_gene, is_de)
  de_genes <- dplyr::group_by(de_genes, ens_gene)
  de_genes <- dplyr::summarize(de_genes, is_de = any(is_de))

  de_genes <- dplyr::rename(de_genes, target_id = ens_gene)
  de_genes <- dplyr::mutate(de_genes, log_fc = NA)

  list(de_info = de_info, de_genes = de_genes)
}

#' Load the differential expression truth (DEPRECATED)
#'
#' This function loads the isoform level information and gene level information.
#' The gene level differential expression is defined by at least 1 isoform of the gene being differentially expressed.
#'
#' @param sim_name the name of the simulation (with respect to the results location)
#' @param ttg the 'target_to_gene_mapping' coming from ensembl
#' @return a list containing the gene and isoform differential expression information
get_de_info_old <- function(sim_name, ttg) {
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

#' Get gene counts from featureCounts
#'
#' Assumes \code{featureCounts} has been run and the results are in the appropriate place.
#'
#' @param sim_name the name of the simulation (with respect to the results location)
#' @param sim_num the integer number of the simulation (with respect to \code{sim_name})
#' @return a matrix containing gene counts
get_gene_counts <- function(sim_name, sim_num = 1) {
  base_dir <- file.path('../results', sim_name)

  experiment_name <- paste0('exp_', sim_num)
  counts <- read.table(file.path(base_dir, experiment_name, 'counts.tsv'),
    header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  counts <- as.matrix(counts)
  colnames(counts) <- sub('X', 'sample_', colnames(counts))

  counts
}
