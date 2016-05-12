
library("SRAdb")
library("dplyr")

sqlfile <- "../../db/SRAmetadb.sqlite"
if(!file.exists(sqlfile)) {
  new_file <- getSRAdbFile()
  file.rename(new_file, sqlfile)
}

sra_con <- dbConnect(SQLite(), sqlfile)

project_id <- 'SRP004777'

all_samples <- dbGetQuery(sra_con, sprintf("SELECT
  e.experiment_ID, r.run_ID, r.run_accession, e.library_strategy,
  e.library_name, e.experiment_accession, s.spots, e.title
  FROM experiment as e, run as r, sra as s
  WHERE
  e.study_accession='%s' AND
  r.experiment_accession=e.experiment_accession AND
  s.experiment_accession=e.experiment_accession", project_id))

all_samples <- dplyr::mutate(all_samples,
  strain = simplify2array(strsplit(library_name, '_'))[1, ])

###
# outputs
###
write(all_samples$run_accession, file = '../metadata/accession.txt')

write.csv(all_samples, file = '../metadata/experiment.csv', row.names = FALSE)
