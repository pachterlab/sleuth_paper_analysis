library("dplyr")

# This is information from the 1000 genomes website:
# http://www.1000genomes.org/data
#
# We use this to get sample => sex mapping
samp_info_1k <- read.csv("../metadata/20130606_sample_info.csv", header = TRUE,
  stringsAsFactors = FALSE)
samp_info_1k <- samp_info_1k %>%
  select(sample = Sample, population = Population,
    population_desc = Population.Description, sex = Gender)

# Let's now look at the samples that have RNA-Seq data
# According to: http://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-3/?keywords=E-GEUV*
# E-GEU-1 contains only the sequencing that passed quality filters
# E-GEU-3 contains all the sequencing

egeu1 <- read.table("../metadata/E-GEUV-1.sdrf.txt", sep = "\t", header = TRUE,
  stringsAsFactors = FALSE)

egeu1 <- egeu1 %>%
  select(
    sample = Source.Name,
    population = Characteristics.population.,
    performer = Performer,
    factor_population = Factor.Value.population.,
    lab = Factor.Value.laboratory.
    ) %>%
  distinct()

all.equal(egeu1$population, egeu1$factor_population)
# let's drop 'factor_population' because it's the same as 'population'
egeu1 <- select(egeu1, -factor_population)

select(egeu1, performer, lab) %>%
  distinct()

# let's also drop "lab" and rename performer "lab"
egeu1 <- egeu1 %>%
  select(-lab) %>%
  rename(lab = performer)

################################################################################

# We can now join the two metadata tables
geu_meta <- inner_join(samp_info_1k, egeu1, by = c("sample"))
all.equal(nrow(geu_meta), nrow(egeu1)) # looks good

all.equal(geu_meta$population.x, geu_meta$population.y)
# looks good, let's drop one column
geu_meta <- geu_meta %>%
  select(-population.y)

geu_meta <- geu_meta %>%
  rename(population = population.x)

# high-level stats
geu_meta %>%
  group_by(population, sex) %>%
  summarise(total = n())

# high-level stats
geu_meta %>%
  group_by(population) %>%
  summarise(total = n())

# high-level stats
geu_meta %>%
  group_by(population, sex, lab) %>%
  summarise(total = n()) %>%
  head(nrow(.))

save(geu_meta, file = "../metadata/geu_meta.RData")
