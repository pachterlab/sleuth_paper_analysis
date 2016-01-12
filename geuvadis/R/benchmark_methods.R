library("Biobase")
library("DESeq2")
library("EBSeq")
library("edgeR")
library("limma")

get_human_gene_names <- function() {
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",
    host="www.ensembl.org")
  ttg <- biomaRt::getBM(
    attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
    mart = mart)
  ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
    ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

  ttg
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

transcript_gene_mapping <- get_human_gene_names()

run_sleuth_prep <- function(sample_info, max_bootstrap = 30, ...) {
  so <- sleuth_prep(sample_info, ~ condition, max_bootstrap = max_bootstrap,
    target_mapping = transcript_gene_mapping, ...)
  so <- sleuth_fit(so)

  so
}

# run_sleuth_wt <- function(obj) {
#   obj <- sleuth_wt(obj, 'conditionB')
#
#   sleuth_results(obj, 'conditionB')[, c('target_id', 'pval', 'qval')]
# }
#
# run_sleuth_lrt <- function(obj) {
#   obj <- sleuth_fit(obj, ~ 1, 'reduced')
#   obj <- sleuth_lrt(obj, 'reduced', 'full')
#
#   sleuth_results(obj, 'reduced:full', test_type = 'lrt')[,
#     c('target_id', 'pval', 'qval')]
# }
#
run_sleuth <- function(sample_info, max_bootstrap = 30, lift_genes = FALSE, ...) {
  so <- run_sleuth_prep(sample_info, max_bootstrap = max_bootstrap, ...)
  so <- sleuth_wt(so, 'conditionB')
  so <- sleuth_fit(so, ~ 1, 'reduced')
  so <- sleuth_lrt(so, 'reduced', 'full')

  res <- NULL
  if (lift_genes) {
    lrt <- get_gene_lift(so, 'reduced:full', test_type = 'lrt')
    wt <- get_gene_lift(so, 'conditionB', test_type = 'wt')
    res <- list(sleuth.lrt = lrt, sleuth.wt = wt)
  } else {
    lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt')[,
      c('target_id', 'pval', 'qval')]
    wt <- sleuth_results(so, 'conditionB')[, c('target_id', 'pval', 'qval')]
    res <- list(sleuth.lrt = lrt, sleuth.wt = wt)
  }

  res
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
    dplyr::select(result, target_id, pval, qval)
  })
}

rename_target_id <- function(df, as_gene = FALSE) {
  if (as_gene) {
    dplyr::rename(df, ens_gene = target_id)
  } else {
    df
  }
}

# The code below is a slightly modified version of the code from `DESeq2paper`
# http://www-huber.embl.de/DESeq2paper/
runDESeq2 <- function(e, as_gene = TRUE) {
  dds <- DESeqDataSetFromMatrix(exprs(e), DataFrame(pData(e)), ~ condition)
  dds <- DESeq(dds,quiet=TRUE)
  res <- results(dds)
  beta <- res$log2FoldChange
  pvals <- res$pvalue
  padj <- res$padj
  pvals[is.na(pvals)] <- 1
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj[is.na(padj)] <- 1

  rename_target_id(
    data.frame(target_id = rownames(res),
      pval = pvals, qval = padj, beta = beta,
      stringsAsFactors = FALSE),
    as_gene = as_gene)
}

runEdgeR <- function(e, as_gene = TRUE) {
  design <- model.matrix(~ pData(e)$condition)
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  dgel <- estimateGLMCommonDisp(dgel, design)
  dgel <- estimateGLMTrendedDisp(dgel, design)
  dgel <- estimateGLMTagwiseDisp(dgel, design)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  predbeta <- predFC(exprs(e), design, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  predbeta10 <- predFC(exprs(e), design, prior.count=10, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  pvals <- edger.lrt$table$PValue
  pvals[rowSums(exprs(e)) == 0] <- NA
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

runVoom <- function(e, as_gene = TRUE) {
  design <- model.matrix(~ condition, pData(e))
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  v <- voom(dgel,design,plot=FALSE)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  tt <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")
  pvals <- tt$P.Value
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1

  rename_target_id(data.frame(target_id = rownames(tt),
    pvals = pvals,
    qval = padj,
    beta = tt$logFC,
    stringsAsFactors = FALSE),
  as_gene = as_gene)
}

runEBSeq <- function(e, as_gene = TRUE) {
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
