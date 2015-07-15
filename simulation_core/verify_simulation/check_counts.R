truth <- read.table("../geuvadis/sims/3_3_1_1/sample_1/1/1.isoforms.results", header = TRUE,
  stringsAsFactors = FALSE)

sims <- read.table("../geuvadis/sims/3_3_1_1/sample_1/1/sim_1_1.counts", header = FALSE)
colnames(sims) <- c("sid", "count")

sims2 <- read.table("../geuvadis/sims/3_3_1_1/sample_1/1/sim_1_2.counts", header = FALSE)
colnames(sims2) <- c("sid", "count")

truth <- truth %>%
  mutate(sid = 1:nrow(.)) %>%
  select(sid, count = expected_count) %>%
  filter(count > 0)

all.equal(truth, sims)
all.equal(truth, sims2)
