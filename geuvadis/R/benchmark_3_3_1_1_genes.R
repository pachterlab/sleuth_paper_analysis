library("cowplot")
library("sleuth")
library("mamabear")

all_ones <- function(x) {
  p <- ncol(x)
  sf <- rep.int(1, p)
  names(sf) <- colnames(x)

  sf
}

sim_name <- '3_3_1_1_1'
sims_dir <- file.path('../sims', sim_name)
base_dir <- file.path('../results', sim_name)

# de_info contains the 'true' fold change as well as which transcripts are DE
de_info <- read.table(gzfile(file.path(sims_dir, "de_info.tsv.gz")),
  header = TRUE, stringsAsFactors = FALSE)

kal_dirs <- file.path(base_dir, "exp_1", 1:6, "kallisto")

s2c <- data.frame(sample = paste0("sample_", 1:6),
  condition = rep(c("A", "B"), each = 3), stringsAsFactors = FALSE)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

so <- sleuth_prep(s2c, ~ condition, norm_fun_counts = all_ones, max_bootstrap = 30)
# so <- sleuth_prep(s2c, ~ condition, norm_fun_counts = all_ones, min_prop = 0.5, min_reads = 10)

so <- sleuth_fit(so)

models(so)

so <- sleuth_wt(so, 'conditionB')
models(so)

# sleuth_live(so)

sres <- sleuth_results(so, 'conditionB') %>%
  dplyr::select(target_id, pval, qval)

so <- sleuth_fit(so, ~1, "reduced")

so <- sleuth_lrt(so, "reduced", "full")
s_ratio <- sleuth_results(so, 'reduced:full', test_type = 'lrt')

save(s_ratio, t2g, de_info, file = 'gene_data.RData')

################################################################################
load('gene_data.RData')

de_genes <- inner_join(de_info, t2g, by = 'target_id')
de_genes <- de_genes %>%
  group_by(ens_gene) %>%
  summarize(is_de = any(is_de))
de_genes <- rename(de_genes, target_id = ens_gene)
de_genes <- mutate(de_genes, log_fc = NA)


s_ratio_genes <- left_join(s_ratio, t2g, by = "target_id")
gene_lift <- s_ratio_genes %>%
  group_by(ens_gene) %>%
  do({
    min_index <- which.min(.$qval)
    result <- .[min_index, ]
    result <- select(result, -target_id)
    result <- rename(result, target_id = ens_gene)
    })

de_bench <- new_de_benchmark(
  list(
    gene_lift
    ),
  c(
    "sleuth lrt"
    ), de_genes)

fdr_nde_plot(de_bench, FALSE) +
  # xlim(0, 2000) +
  xlim(0, 4000) +
  ylim(0, 0.10) +
  theme(legend.position = c(0.1, 0.80))+
  xlab('number of genes called DE') +
  ylab('FDR')
