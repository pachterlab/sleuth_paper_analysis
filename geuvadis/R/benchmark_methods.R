# This file contains functions to run several different differential expression
# methods, usually simply by providing the "count matrix"

library("Biobase")
library("DESeq2")
library("EBSeq")
library("edgeR")
library("limma")
library("sleuth")

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

#' parse the simulation name
#'
#' Get back a named list with the simulation info
#'
#' @param sim_name a simulation name in the form 'isoform_3_3_20_1_1'
parse_simulation <- function(sim_name) {
  sim_name_split <- strsplit(sim_name, '_')[[1]]

  res <- list()
  res[['type']] <- sim_name_split[[1]]
  res[['a']] <- as.integer(sim_name_split[[2]])
  res[['b']] <- as.integer(sim_name_split[[3]])
  res[['n']] <- as.integer(sim_name_split[[4]])
  res[['seed']] <- as.integer(sim_name_split[[5]])
  res[['sf']] <- as.integer(sim_name_split[[6]])

  res
}

sleuth_filter <- function(mat, ...) {
  apply(mat, 1, sleuth::basic_filter, ...)
}

edgeR_filter <- function(mat, ...) {
  rowSums(cpm(mat) > 1) >= 2
}

DESeq2_filter <- function(mat, ...) {
  rowSums(mat) > 1
}

#' create benchmark objects against a reference
#'
#' Given several differential expression tables, compare them against a
#' reference in order to make benchmark objects using \code{new_de_benchmark}
#' @param results_list a named list with a set of results (differential expression tables)
#' @param de_info the "truth" that was used to simulate differential expression
#' @param reference a character string vector with a set of references
compare_reference <- function(results_list, de_info,
  reference = c('sleuth.wt', 'sleuth.lrt')) {

  other_methods <- names(results_list)[!(names(results_list) %in% reference)]
  if (length(other_methods) == 0) {
    stop('No other methods (or missing names) in results_list')
  }

  res <- lapply(other_methods,
    function(m) {
      ms <- c(reference, m)
      new_de_benchmark(results_list[ms], ms, de_info)
    })
  names(res) <- other_methods

  res
}

#' @param sim_name a simulation name such as 'isoform_3_3_20_1_1'
#' @param which_sample which sample (replication) to load (an integer from 1 to N)
#' @param method_filtering if \code{TRUE}, use the methods own filtering.
#' Otherwise, use the filtering provided from sleuth.
#' @param ... additional arguments passed to \code{run_sleuth}
#' NOTE: cuffdiff uses a filter method internally and it cannot be changed
load_isoform_results <- function(
  sim_name,
  which_sample,
  method_filtering = FALSE,
  ...) {
  sim <- parse_simulation(sim_name)

  if (which_sample > sim$n) {
    stop('which_sample must be less than the total number of replications: ', sim$n)
  }

  which_sample <- as.integer(which_sample)
  n <- sim$a + sim$b

  sim_info <- get_de_info(sim_name, which_sample, transcript_gene_mapping)
  de_info <- sim_info$de_info
  de_genes <- sim_info$de_genes

  kal_dirs <- file.path('..', 'sims', sim_name, paste0("exp_", which_sample),
    1:n, "kallisto")
  sample_to_condition <- get_sample_to_condition(sim$a, sim$b, kal_dirs)

  sir <- run_sleuth(sample_to_condition, gene_mode = NULL, ...)

  isoform_cds <- get_filtered_isoform_cds(sir$so, sample_to_condition,
    method_filtering)
  # sir$so <- NULL

  gene_methods <- list(
    # DESeq2 = runDESeq2,
    # edgeR = runEdgeR,
    limmaVoom = runVoom
    # edgerRobust = runEdgeRRobust,
    # EBSeq = runEBSeq
    )
  isoform_results <- lapply(gene_methods,
    function(f) {
      f(isoform_cds, FALSE, method_filtering)
    })
  all_results <- c(Filter(is.data.frame, sir) , isoform_results)

  # TODO: adjust the fdr based off of the filtering scheme
  # e.g. take the intersection of the tests and recompute the fdr
  # This ensures that the calibration tests can be comparable

  # cr <- get_cuffdiff(
  #   file.path('..', 'sims', sim_name, paste0("exp_", which_sample),
  #     'results', 'cuffdiff')
  #   )[['isoform']]
  # all_results <- c(all_results, list(Cuffdiff2 = cr))

  all_results
}

#' @param sim_name a simulation name such as 'isoform_3_3_20_1_1'
#' @param which_sample which sample (replication) to load (an integer from 1 to N)
#' @param method_filtering if \code{TRUE}, use the methods own filtering.
#' Otherwise, use the filtering provided from sleuth.
#' @param ... additional arguments passed to \code{run_sleuth}
#' NOTE: cuffdiff uses a filter method internally and it cannot be changed
load_isoform_results_intersect <- function(
  sim_name,
  which_sample,
  method_label,
  method_fit_function,
  ...) {
  sim <- parse_simulation(sim_name)

  if (which_sample > sim$n) {
    stop('which_sample must be less than the total number of replications: ', sim$n)
  }

  which_sample <- as.integer(which_sample)
  n <- sim$a + sim$b

  sim_info <- get_de_info(sim_name, which_sample, transcript_gene_mapping)
  de_info <- sim_info$de_info
  de_genes <- sim_info$de_genes

  kal_dirs <- file.path('..', 'sims', sim_name, paste0("exp_", which_sample),
    1:n, "kallisto")
  sample_to_condition <- get_sample_to_condition(sim$a, sim$b, kal_dirs)

  # simply do this so we can read in the data
  so_data <- sleuth_prep(sample_to_condition, ~1, max_bootstrap = 3)
  obs_raw <- sleuth:::spread_abundance_by(so_data$obs_raw, "est_counts")
  rm(so_data)

  s_which_filter <- sleuth_filter(obs_raw, ...)
  # method_filters <- list()
  # method_filters$edgeR <- edgeR_filter(obs_raw)
  # method_filters$edgeR <- method_filters$edgeR & s_which_filter
  #
  # method_filters$DESeq2 <- DESeq2_filter(obs_raw)
  # method_filters$DESeq2 <- method_filters$DESeq2 & s_which_filter

  # limma

  # sir_limma <- run_sleuth(sample_to_condition, gene_mode = NULL,
  #   filter_target_id = names(which(s_which_filter)), ...)
  #   # filter_target_id = names(which(method_filters$edgeR)), ...)
  # sir_limma <- Filter(is.data.frame, sir_limma)
  # # cds <- make_count_data_set(obs_raw[method_filters$edgeR, ], sample_to_condition)
  # cds <- make_count_data_set(obs_raw[s_which_filter, ], sample_to_condition)
  # limma_results <- runVoom(cds, FALSE, FALSE)

  method_result <- method_fit_function(obs_raw, sample_to_condition,
    s_which_filter)
  sir <- run_sleuth(sample_to_condition, gene_mode = NULL,
    filter_target_id = method_result$filter)

  # isoform_cds <- get_filtered_isoform_cds(sir$so, sample_to_condition,
  #   method_filtering)
  # # sir$so <- NULL
  #
  # gene_methods <- list(
  #   # DESeq2 = runDESeq2,
  #   # edgeR = runEdgeR,
  #   limmaVoom = runVoom
  #   # edgerRobust = runEdgeRRobust,
  #   # EBSeq = runEBSeq
  #   )
  # isoform_results <- lapply(gene_methods,
  #   function(f) {
  #     f(isoform_cds, FALSE, method_filtering)
  #   })
  # all_results <- c(Filter(is.data.frame, sir) , isoform_results)

  # TODO: adjust the fdr based off of the filtering scheme
  # e.g. take the intersection of the tests and recompute the fdr
  # This ensures that the calibration tests can be comparable

  # cr <- get_cuffdiff(
  #   file.path('..', 'sims', sim_name, paste0("exp_", which_sample),
  #     'results', 'cuffdiff')
  #   )[['isoform']]
  # all_results <- c(all_results, list(Cuffdiff2 = cr))

  all_results <- Filter(is.data.frame, sir)
  all_results[[method_label]] <- method_result$result

  all_results
}

limma_filter_and_run <- function(count_matrix, stc, sleuth_filter) {
  cds <- make_count_data_set(count_matrix[sleuth_filter, ], stc)
  res <- runVoom(cds, FALSE, FALSE)
  sleuth_filter <- names(which(sleuth_filter))
  list(result = res, filter = sleuth_filter)
}

# DEPRECATED
DESeq2_filter_and_run <- function(count_matrix, stc, sleuth_filter) {
  # we should check if taking the intersection of results does better
  count_matrix <- round(count_matrix)
  mode(count_matrix) <- 'integer'
  which_targets <- DESeq2_filter(count_matrix)
  cds <- make_count_data_set(count_matrix[which_targets, ], stc)
  res <- runDESeq2(cds, FALSE, FALSE)

  sleuth_filter <- sleuth_filter & which_targets
  sleuth_filter <- names(which(sleuth_filter))

  list(result = res, filter = sleuth_filter)
}

DESeq2_filter_and_run_intersect <- function(count_matrix, stc, sleuth_filter) { # nolint
  # we should check if taking the intersection of results does better
  count_matrix <- round(count_matrix)
  mode(count_matrix) <- 'integer'
  which_targets <- DESeq2_filter(count_matrix)
  sleuth_filter <- sleuth_filter & which_targets
  cds <- make_count_data_set(count_matrix[sleuth_filter, ], stc)
  res <- runDESeq2(cds, FALSE, FALSE)

  sleuth_filter <- names(which(sleuth_filter))

  list(result = res, filter = sleuth_filter)
}

edgeR_filter_and_run <- function(count_matrix, stc, sleuth_filter) {
  count_matrix <- round(count_matrix)
  mode(count_matrix) <- 'integer'
  which_targets <- edgeR_filter(count_matrix)
  sleuth_filter <- sleuth_filter & which_targets

  cds <- make_count_data_set(count_matrix[sleuth_filter, ], stc)
  res <- runEdgeR(cds, FALSE, FALSE)

  sleuth_filter <- names(which(sleuth_filter))

  list(result = res, filter = sleuth_filter)
}

cuffdiff_filter_and_run <- function(count_matrix, stc, sleuth_filter) {
  # ignore the count_matrix
  base_name <- dirname(dirname(stc$path[1]))
  base_name <- file.path(base_name, 'results', 'cuffdiff')
  res <- get_cuffdiff(base_name)$isoform
  res <- dplyr::filter(res, status == 'OK')

  sleuth_filter <- names(which(sleuth_filter))

  sleuth_filter <- res$target_id

  list(result = res, filter = sleuth_filter)
}

EBSeq_isoform_filter_and_run <- function(count_matrix, stc, sleuth_filter) {
  # a subset of things don't have gene names in biomaRt
  sleuth_filter <- intersect(names(which(sleuth_filter)),
    names(NG_LIST$IsoformNgTrun))
  cds <- make_count_data_set(count_matrix[sleuth_filter, ], stc)
  res <- EBSeq_isoform(cds, FALSE, FALSE)

  list(result = res, filter = sleuth_filter)
}

# transcript_gene_mapping <- get_human_gene_names()

#' generate `sample_to_condition` that sleuth is expecting
#'
#' @param n_a the number of samples in condition A
#' @param n_a the number of samples in condition B
#' @param kal_dirs if not NULL, then add the appropriate `path` column
#' @return a \code{data.frame} in the proper sleuth form
get_sample_to_condition <- function(n_a, n_b, kal_dirs = NULL) {
  n <- n_a + n_b

  sample_to_condition <- data.frame(
    sample = paste0("sample_", 1:n),
    condition = factor(c(rep("A", n_a), rep("B", n_b))),
    stringsAsFactors = FALSE)

  if (!is.null(kal_dirs)) {
    sample_to_condition <- dplyr::mutate(sample_to_condition, path = kal_dirs)
  }
  rownames(sample_to_condition) <- sample_to_condition$sample

  sample_to_condition
}

#' Generate equal size factors
#'
#' Generate size factors all equal to 1. Helpful for simulations.
#'
#' @param x the count matrix
#' @return a vector of all ones
all_ones <- function(x) {
  p <- ncol(x)
  sf <- rep.int(1, p)
  names(sf) <- colnames(x)

  sf
}

#' @param counts \code{matrix} of counts with transcripts on the rows and
#' samples on the columns
#' @param sample_info \code{data.frame} of sample information with at least a
#' column called \code{condition}
make_count_data_set <- function(counts, sample_info) {
  ExpressionSet(counts, AnnotatedDataFrame(sample_info))
}

run_sleuth_prep <- function(sample_info, max_bootstrap = 30,
  filter_target_id = NULL, ...) {
  so <- sleuth_prep(sample_info, ~ condition, max_bootstrap = max_bootstrap,
    target_mapping = transcript_gene_mapping,
    filter_target_id = filter_target_id,
    ...)
  so <- sleuth_fit(so)

  so
}

#' @param gene_mode if NULL, do isoform mode, if 'lift' do gene lifting, if 'aggregate', do gene aggregation
run_sleuth <- function(sample_info,
  max_bootstrap = 30,
  gene_mode = NULL,
  filter_target_id = NULL,
  ...) {
  so <- run_sleuth_prep(sample_info, max_bootstrap = max_bootstrap,
    filter_target_id = filter_target_id, ...)
  so <- sleuth_wt(so, 'conditionB')
  so <- sleuth_fit(so, ~ 1, 'reduced')
  so <- sleuth_lrt(so, 'reduced', 'full')

  res <- NULL
  if (is.null(gene_mode)) {
    lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt',
      show_all = FALSE)[, c('target_id', 'pval', 'qval')]
    wt <- sleuth_results(so, 'conditionB',
      show_all = FALSE)[, c('target_id', 'pval', 'qval')]
    res <- list(sleuth.lrt = lrt, sleuth.wt = wt)
  } else if (gene_mode == 'lift') {
    # test every isoform
    lrt <- get_gene_lift(so, 'reduced:full', test_type = 'lrt')
    wt <- get_gene_lift(so, 'conditionB', test_type = 'wt')
    res <- list(sleuth.lrt = lrt, sleuth.wt = wt)
  } else if (gene_mode == 'aggregate') {
    stop('TODO: implement me!')
  } else {
    stop('Unrecognized mode for "run_sleuth"')
  }

  res$so <- so

  res
}

# if method_filtering is false, then use the sleuth filter
get_filtered_isoform_cds <- function(so, stc, method_filtering = FALSE) {
  pass_filt_names <- so$filter_df[['target_id']]

  obs_raw <- sleuth:::spread_abundance_by(so$obs_raw, "est_counts")
  if (!method_filtering) {
    obs_raw <- obs_raw[pass_filt_names,]
  }
  isoform_cds <- make_count_data_set(round(obs_raw), stc)

  isoform_cds
}

#' @param obj a sleuth object
#' @param ... additional arguments to \code{sleuth_gene_table}
#' @example \code{get_gene_lift(so, 'reduced:full', test_type = 'lrt')}
get_gene_lift <- function(obj, ...) {
  sgt <- sleuth_gene_table(obj, ...)
  sgt <- group_by(sgt, ens_gene)
  do(sgt, {
    min_index <- which.min(.$qval)
    result <- .[min_index, ]
    result <- dplyr::select(result, -target_id)
    result <- dplyr::rename(result, target_id = ens_gene)
    # dplyr::select(result, target_id, pval, qval)
    result
  })
}

rename_target_id <- function(df, as_gene = FALSE) {
  if (as_gene) {
    dplyr::rename(df, ens_gene = target_id)
  } else {
    df
  }
}

get_cuffdiff <- function(results_path) {
  isoform_de <- suppressWarnings(
    data.table::fread(file.path(results_path, 'isoform_exp.diff'),
      header = TRUE, stringsAsFactors = FALSE, sep = '\t', data.table = FALSE))
  isoform_de <- dplyr::select(isoform_de,
    target_id = test_id,
    ens_gene = gene_id,
    status,
    pval = p_value,
    qval = q_value,
    beta = one_of("log2(fold_change)")
    )

  gene_de <- suppressWarnings(
    data.table::fread(file.path(results_path, 'gene_exp.diff'),
      header = TRUE, stringsAsFactors = FALSE, sep = '\t', data.table = FALSE))
  gene_de <- dplyr::select(gene_de,
    target_id = test_id,
    ens_gene = gene_id,
    status,
    pval = p_value,
    qval = q_value,
    beta = one_of("log2(fold_change)")
    )

  list(isoform = isoform_de, gene = gene_de)
}


EBSeq_isoform <- function(e, as_gene = TRUE, method_filtering = FALSE) {
  # this assumes that we have a global variable `NG_LIST`
  # that has been generated using EBSeq::GetNg(transcript_gene_mapping...)
  isoforms_per_gene <- NG_LIST$IsoformNgTrun
  isoforms_per_gene <- isoforms_per_gene[rownames(exprs(e))]

  sizes <- MedianNorm(exprs(e))
  out <- capture.output({
    suppressMessages({
      res <- EBTest(Data = exprs(e),
        NgVector = isoforms_per_gene,
        Conditions = pData(e)$condition,
        sizeFactors = sizes,
        maxround = 15)
    })
  })
  print("Alpha")
  print(res$Alpha)
  print("Beta")
  print(res$Beta)
  n_it <- nrow(res$Alpha)
  convergence_limit <- 1e-3
  cat("Alpha converged: ",
    abs(res$Alpha[n_it, 1] - res$Alpha[n_it - 1, 1]) < convergence_limit, "\n")
  cat("Beta converged: ",
    abs(res$Beta[n_it, ] - res$Beta[n_it - 1, ]) < convergence_limit, "\n")
  padj <- rep(1, nrow(exprs(e)))
  # we use 1 - PPDE for the FDR cutoff as this is recommended in the EBSeq vignette
  padj[match(rownames(res$PPMat), rownames(e))] <- res$PPMat[,"PPEE"]
  beta <- rep(0, nrow(exprs(e)))

  rename_target_id(
    data.frame(target_id = rownames(res$PPMat),
      pval = padj,
      qval = padj,
      beta = beta
      ),
    as_gene = as_gene)
}

# The code below is a slightly modified version of the code from `DESeq2paper`
# http://www-huber.embl.de/DESeq2paper/
runDESeq2 <- function(e, as_gene = TRUE, compute_filter = FALSE) {
  dds <- DESeqDataSetFromMatrix(exprs(e), DataFrame(pData(e)), ~ condition)
  if (compute_filter) {
    # Section 1.3.6 in DESeq2 vignette
    # https://www.bioconductor.org/packages/3.3/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
    dds <- dds[ rowSums(counts(dds)) > 1, ]
  }
  dds <- DESeq(dds,quiet=TRUE)
  res <- results(dds)
  beta <- res$log2FoldChange
  pvals <- res$pvalue
  padj <- res$padj
  pvals[is.na(pvals)] <- 1
  # pvals[rowSums(exprs(e)) == 0] <- NA
  padj[is.na(padj)] <- 1

  rename_target_id(
    data.frame(target_id = rownames(res),
      pval = pvals, qval = padj, beta = beta,
      stringsAsFactors = FALSE),
    as_gene = as_gene)
}

runEdgeR <- function(e, as_gene = TRUE, compute_filter = FALSE) {
  design <- model.matrix(~ pData(e)$condition)
  dgel <- DGEList(exprs(e))
  if (compute_filter) {
    # Section 2.6 in edgeR vignette
    # https://www.bioconductor.org/packages/3.3/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
    keep <- rowSums(cpm(dgel) > 1) >= 2
    dgel <- dgel[keep, , keep.lib.sizes=FALSE]
  }
  dgel <- calcNormFactors(dgel)
  dgel <- estimateGLMCommonDisp(dgel, design)
  dgel <- estimateGLMTrendedDisp(dgel, design)
  dgel <- estimateGLMTagwiseDisp(dgel, design)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  # predbeta <- predFC(exprs(e), design, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  # predbeta10 <- predFC(exprs(e), design, prior.count=10, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  predbeta <- predFC(dgel$counts, design, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  predbeta10 <- predFC(dgel$counts, design, prior.count=10, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  pvals <- edger.lrt$table$PValue
  # pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1

  rename_target_id(
    data.frame(
      target_id = rownames(edger.lrt$table),
      pval = pvals,
      qval = padj,
      beta = log2(exp(1)) * edger.fit$coefficients[,"pData(e)$conditionB"],
      predbeta = predbeta[,"pData(e)$conditionB"],
      predbeta10 = predbeta10[,"pData(e)$conditionB"],
      stringsAsFactors = FALSE),
    as_gene = as_gene)
}

runEdgeRRobust <- function(e, as_gene = TRUE) {
  design <- model.matrix(~ pData(e)$condition)
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  # settings for robust from robinson_lab/edgeR_robust/robust_simulation.R
  dgel <- estimateGLMRobustDisp(dgel, design, maxit=6)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  predbeta <- predFC(exprs(e), design, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  pvals <- edger.lrt$table$PValue
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  # list(pvals=pvals, padj=padj, beta=log2(exp(1)) * edger.fit$coefficients[,"pData(e)$conditionB"],
  #      predbeta=predbeta[,"pData(e)$conditionB"])
  rename_target_id(
    data.frame(
      target_id = rownames(edger.lrt$table),
      pval = pvals,
      qval = padj,
      beta=log2(exp(1)) * edger.fit$coefficients[,"pData(e)$conditionB"],
      predbeta=predbeta[,"pData(e)$conditionB"],
      stringsAsFactors = FALSE),
    as_gene = as_gene)
}

runVoom <- function(e, as_gene = TRUE, compute_filter = FALSE) {
  design <- model.matrix(~ condition, pData(e))
  dgel <- DGEList(exprs(e))
  if (compute_filter) {
    # Section 2.6 in edgeR vignette
    # https://www.bioconductor.org/packages/3.3/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
    keep <- rowSums(cpm(dgel) > 1) >= 2
    dgel <- dgel[keep, , keep.lib.sizes=FALSE]
  }
  dgel <- calcNormFactors(dgel)
  v <- voom(dgel,design,plot=FALSE)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  tt <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")
  pvals <- tt$P.Value
  # pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1

  rename_target_id(data.frame(target_id = rownames(tt),
    pval = pvals,
    qval = padj,
    beta = tt$logFC,
    stringsAsFactors = FALSE),
  as_gene = as_gene)
}

runEBSeq <- function(e, as_gene = TRUE, method_filtering = FALSE) {
  sizes <- MedianNorm(exprs(e))
  out <- capture.output({
    suppressMessages({
      res <- EBTest(Data = exprs(e),
                    Conditions = pData(e)$condition,
                    sizeFactors = sizes,
                    maxround = 5)
    })
  })
  padj <- rep(1, nrow(exprs(e)))
  # we use 1 - PPDE for the FDR cutoff as this is recommended in the EBSeq vignette
  padj[match(rownames(res$PPMat), rownames(e))] <- res$PPMat[,"PPEE"]
  beta <- rep(0, nrow(exprs(e)))

  rename_target_id(
    data.frame(target_id = rownames(res$PPMat),
      pval = padj,
      qval = padj,
      beta = beta
      ),
    as_gene = as_gene)
}

# these methods used right now and need to be updated

runDSS <- function(e) {
  X <- as.matrix(exprs(e))
  colnames(X) <- NULL
  designs <- as.character(pData(e)$condition)
  seqData <- newSeqCountSet(X, designs)
  seqData <- estNormFactors(seqData)
  seqData <- estDispersion(seqData)
  result <- waldTest(seqData, "B", "A")
  result <- result[match(rownames(seqData),rownames(result)),]
  pvals <- result$pval
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj, beta=( log2(exp(1)) * result$lfc ))
}

runDESeq2Outliers <- function(e, retDDS=FALSE) {
  dds <- DESeqDataSetFromMatrix(exprs(e), DataFrame(pData(e)), ~ condition)
  ddsDefault <- DESeq(dds, quiet=TRUE)
  ddsNoRepl <- ddsDefault
  if (ncol(e) >= 14) {
    # insert original maximum Cook's distances
    # so the rows with replacement will be filtered
    # this avoid re-running with minReplicateForReplace=Inf
    mcols(ddsNoRepl)$maxCooks <- apply(assays(ddsNoRepl)[["cooks"]], 1, max)
  }
  resDefault <- results(ddsDefault)
  resNoFilt <- results(ddsDefault, cooksCutoff=FALSE)
  resNoRepl <- results(ddsNoRepl)
  resList <- list("DESeq2"=resDefault, "DESeq2-noFilt"=resNoFilt, "DESeq2-noRepl"=resNoRepl)
  resOut <- lapply(resList, function(res) {
    pvals <- res$pvalue
    padj <- res$padj
    pvals[is.na(pvals)] <- 1
    pvals[rowSums(exprs(e)) == 0] <- NA
    padj <- p.adjust(pvals,method="BH")
    padj[is.na(padj)] <- 1
    list(pvals=pvals, padj=padj)
  })
  return(resOut)
}
