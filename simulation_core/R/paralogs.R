paralogs <- read.table('../data/dgd_Hsa_all_v71.tsv', stringsAsFactors = FALSE,
  sep = '\t', header = TRUE)

paralogs <- dplyr::select(paralogs, group_id, ens_gene = ENS_ID,
  external_gene = Name, num_genes = NB_Genes)



pl_genes <- inner_join(de_genes, paralogs, by = c('target_id' = 'ens_gene'))

temporary <- tmp %>%
  group_by(group_id) %>%
  filter(any(is_de))
