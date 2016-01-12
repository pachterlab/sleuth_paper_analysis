# XXX: be sure to set 'sim_name' appropriately before sourcing this file

library("sleuth")
library("mamabear")

sims_dir <- file.path('../sims', sim_name)
base_dir <- file.path('../results', sim_name)

# de_info contains the 'true' fold change as well as which transcripts are DE
de_info <- read.table(gzfile(file.path(sims_dir, "de_info.tsv.gz")),
  header = TRUE, stringsAsFactors = FALSE)

de_genes <- select(de_info, ens_gene, is_de) %>%
  group_by(ens_gene) %>%
  summarize(is_de = any(is_de))
de_genes <- dplyr::rename(de_genes, target_id = ens_gene)
de_genes <- dplyr::mutate(de_genes, log_fc = NA)
