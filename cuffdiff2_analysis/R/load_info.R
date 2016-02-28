get_human_gene_names <- function() {
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "hsapiens_gene_ensembl",
    # host = "may2015.archive.ensembl.org")
    host = "ensembl.org")
  ttg <- biomaRt::getBM(
    attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
    mart = mart)
  ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
    ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

  ttg
}

transcript_gene_mapping <- get_human_gene_names()
transcript_gene_mapping <- dplyr::select(transcript_gene_mapping,
  target_id, ens_gene)

info <- read.table('../metadata/hiseq_info.txt',
  header = TRUE, stringsAsFactors = FALSE)
info <- dplyr::select(info, sample = run_accession, condition)

info <- dplyr::mutate(info,
  path = file.path('../results/paired', info$sample, 'kallisto'))
