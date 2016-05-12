metadata <- read.csv('../metadata/experiment.csv', stringsAsFactors = FALSE)

extract_metadata <- function(library_name) {
  ret <- lapply(strsplit(library_name, '_'),
    function(x) {
      data.frame(strain = x[1], experiment = x[2], lane = x[3],
        stringsAsFactors = FALSE)
    })
  dplyr::bind_rows(ret)
}

metadata <- dplyr::select(metadata, -strain)
metadata <- dplyr::bind_cols(metadata, extract_metadata(metadata$library_name))
metadata <- dplyr::select(metadata, run_accession, library_name, strain,
  experiment, lane)

metadata <- dplyr::mutate(metadata,
  path = file.path('..', 'results', 'single', run_accession, 'kallisto'))

mouse_genes <- get_mouse_gene_names()

library('sleuth')

metadata <- dplyr::rename(metadata, sample = run_accession)

so <- sleuth_prep(metadata, ~strain, target_mapping = mouse_genes)
so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt')
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
sleuth_null <- dplyr::filter(sleuth_table, qval > 0.05)
sleuth_live(so)

###
# run DESeq2
###

obs_raw <- sleuth:::spread_abundance_by(so$obs_raw, "est_counts")
obs_raw <- round(obs_raw)

metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
rownames(metadata) <- metadata$sample
dds <- DESeqDataSetFromMatrix(obs_raw, metadata, ~strain)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds,quiet=TRUE)
res <- results(dds)
DESeq2_result <- as.data.frame(res, stringsAsFactors = FALSE)
DESeq2_result <- dplyr::mutate(DESeq2_result, target_id = rownames(DESeq2_result))
DESeq2_result <- dplyr::rename(DESeq2_result, qval = padj)
DESeq2_significant <- dplyr::filter(DESeq2_result, qval <= 0.05)

###
# run voom
###

design <- model.matrix(~ strain, metadata)
dgel <- DGEList(obs_raw)

# Section 2.6 in edgeR vignette
# https://www.bioconductor.org/packages/3.3/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
keep <- rowSums(cpm(dgel) > 1) >= 2
dgel <- dgel[keep, , keep.lib.sizes=FALSE]

dgel <- calcNormFactors(dgel)
v <- voom(dgel,design,plot=FALSE)
fit <- lmFit(v,design)
fit <- eBayes(fit)
tt <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")
voom_result <- data.frame(tt, target_id = rownames(tt))
voom_result <- dplyr::rename(voom_result, qval = adj.P.Val)
voom_significant <- dplyr::filter(voom_result, qval <= 0.05)

###
# intersection between DESeq2 and voom
###

DESeq2_voom <- intersect(DESeq2_significant$target_id,
  voom_significant$target_id)

target_set <- data.frame(target_id = setdiff(sleuth_null$target_id, DESeq2_voom))

target_set <- dplyr::inner_join(sleuth_table, target_set, by = 'target_id')

target_set <- dplyr::arrange(target_set, desc(qval), desc(tech_var))

# Weird things
# ENSMUST00000027800.14
# ENSMUST00000084497.11
# ENSMUST00000031414.14

# ENSMUST00000113388.2

# ENSMUST00000043866.7


kt <- kallisto_table(so, use_filtered = TRUE)

kt <- dplyr::inner_join(kt, head(dplyr::select(target_set, target_id), 20),
  by = 'target_id')

ggplot(kt, aes(sample, est_counts)) +
  geom_point(aes(color = strain)) +
  facet_wrap(~target_id)
