# source("http://bioconductor.org/biocLite.R")
# biocLite("SRAdb")

library("SRAdb")
library("dplyr")

sqlfile <- "../db/SRAmetadb.sqlite"
if(!file.exists(sqlfile)) {
  sqlfile <<- getSRAdbFile()
}

sra_con <- dbConnect(SQLite(), sqlfile)

proj_id <- 'SRP012607'

all_samples <- dbGetQuery(sra_con, sprintf("SELECT
  e.experiment_ID, r.run_ID, r.run_accession, e.library_strategy,
  e.library_name, e.experiment_accession, s.spots, e.title
  FROM experiment as e, run as r, sra as s
  WHERE
  e.study_accession='%s' AND
  r.experiment_accession=e.experiment_accession AND
  s.experiment_accession=e.experiment_accession", proj_id))

extract_info <- function(df) {
  df$library_name %>%
    strsplit(" ") %>%
    lapply(function(x) x[2]) %>%
    unlist() %>%
    strsplit("_") %>%
    lapply(
      function(x) {
        data.frame(condition = x[2], sequencer = x[3],
          sample = sub("rep", "", x[4]),stringsAsFactors = FALSE)
      }) %>%
      rbind_all() %>%
    mutate(run_accession = df$run_accession)
}

sample_info <- all_samples %>%
  extract_info()

all_samples <- all_samples %>%
  select(run_accession, experiment_accession, spots)

all_samples <- all_samples %>%
  left_join(sample_info, by = c("run_accession"))

hiseq_samples <- all_samples %>%
  filter(sequencer == "hiseq")

write(hiseq_samples$run_accession, file = "../metadata/hiseq_accession.txt")
write.table(hiseq_samples, file = "../metadata/hiseq_info.txt",
  row.names = FALSE, quote = FALSE)
