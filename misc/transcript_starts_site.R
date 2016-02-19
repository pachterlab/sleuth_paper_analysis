
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
att <- biomaRt::listAttributes(mart)
tss <- biomaRt::getBM(
  attributes = c(
    "ensembl_transcript_id",
    "transcription_start_site",
    "ensembl_gene_id",
    "external_gene_name"), mart = mart)

tss_summary <- tss %>%
  group_by(ensembl_gene_id, transcription_start_site) %>%
  summarize(n = n())

ggplot(tss_summary, aes(n)) +
  geom_histogram(binwidth = 1)

ggplot(tss_summary, aes(n)) +
  stat_ecdf()
