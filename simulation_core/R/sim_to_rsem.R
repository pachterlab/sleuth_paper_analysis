#' Generate fold changes
#'
#' Generate fold changes, similar to the method used in "Jun Li and Robert
#' Tibshirani (2011) Finding consistent patterns: a nonparametric approach for
#' identifying differential expression in RNA-Seq data."
#'
#' @param valid a bool vector where TRUE represents if the feature passes some
#' filter, FALSE otherwise.
#' @param p_de the proportion of differentially expressed features
#' @param p_up the proportion of differentially expressed features whose
#' expression should go up
#' @return a named list with \code{fold_changes}, \code{n_de}, \code{n_up} and
#' \code{n_down}
gen_foldchanges <- function(valid, p_de = 0.3, p_up = 0.8) {
  which_valid <- which(valid)
  n_de <- ceiling(p_de * length(which_valid))
  n_up <- ceiling(p_up * n_de)
  n_down <- n_de - n_up

  fc <- rep(1.0, length(valid))

  # which are DE?
  w_de <- sample(which_valid, n_de, replace = FALSE)
  w_up <- sample(w_de, n_up)
  w_down <- setdiff(w_de, w_up)

  fc_up <- exp( abs(rnorm(n_up)) )
  fc_down <- exp( -abs(rnorm(n_down)) )

  fc[w_up] <- fc_up
  fc[w_down] <- fc_down

  is_de <- rep(FALSE, length(valid))
  is_de[w_de] <- TRUE

  list(fold_changes = data.frame(valid = valid, fold_change = fc, is_de = is_de),
    n_de = n_de,
    n_up = n_up,
    n_down = n_down)
}

sim_nb <- function(base_mean, size) {
  valid_mean <- which( base_mean > .Machine$double.eps )

  counts <- rep.int(0, length(base_mean))
  counts[valid_mean] <- rnbinom(
    n = length(valid_mean),
    mu = base_mean[valid_mean],
    size = size[valid_mean])

  counts
}

to_rsem <- function(target_id, counts, len, eff_len) {
  data.frame(transcript_id = target_id, gene_id = target_id,
    length = len, effective_length = eff_len,
    expected_count = counts,
    TPM = counts_to_tpm(counts, eff_len),
    FPKM = counts_to_fpkm(counts, eff_len),
    IsoPct = ifelse(counts > 0, 100.00, 0.0),
    stringsAsFactors = FALSE)
}

#' Generate DE counts
#'
#' Generate differential expression counts from transcripts
generate_counts <- function(base_df, fold_changes, out_dir = "gen_sim",
  n_sim = c(4, 4), seed = 42) {
  set.seed(seed)
  size <- base_df$counts / 3
  exp_means <- ceiling( base_df$counts * fold_changes )
  exp_means <- exp_means[, c(rep.int(1, n_sim[1]), rep.int(2, n_sim[2]))]

  counts <- lapply(1:ncol(exp_means),
    function(col)
    {
      sim_nb(exp_means[,col], size)
    })

  rsem_tbl <- lapply(counts,
    function(cnt)
    {
      to_rsem(base_df$target_id, cnt, base_df$length, base_df$eff_len)
    })

  counts <- data.frame(counts)
  colnames(counts) <- 1:ncol(exp_means)
  counts$target_id <- base_df$target_id

  dir.create(out_dir)
  write.table(base_df,
    file = file.path(out_dir, "de.txt"),
    sep = "\t",
    row.names = FALSE, quote = FALSE, col.names = TRUE)
  write.table(counts,
    file = file.path(out_dir, "counts.txt"),
    sep = "\t",
    row.names = FALSE, quote = FALSE, col.names = TRUE)

  invisible(lapply(seq_along(rsem_tbl),
    function(i)
    {
      dir.create(file.path(out_dir, i))
      out_fname <- file.path(out_dir, i, paste0(i, ".isoforms.results"))
      write.table(rsem_tbl[[i]],
        file = out_fname,
        sep = "\t", eol = "\n",
        row.names = FALSE, col.names = TRUE,
        quote = FALSE)
      out_str <- paste("#!/bin/bash",
        "",
        'if [ "$#" -ne 3 ]; then',
        ' echo "Requires two arguments RSEM_REFERENCE RSEM_MODEL_FILE"',
        ' exit 1',
        'fi',
        '',
        "RSEM_REFERENCE=$1",
        "MODEL_FILE=$2",
        "OUT_DIR=$3",
        sep = "\n")
      theta <- 0.0
      total_reads <- sum(rsem_tbl[[i]]$expected_count)
      out_str <- c(out_str, "",
        paste("rsem-simulate-reads", "${RSEM_REFERENCE}", "${MODEL_FILE}",
          paste0("${OUT_DIR}/", i, ".isoforms.results"),
          theta, total_reads, paste0("${OUT_DIR}/sim_", i), "--seed", i, "\n"))
      out_str <- c(out_str,
        "",
        "gzip *.fq", "\n")

      out_str <- paste(out_str, collapse = "\n")

      script_fname <- file.path(out_dir, i, paste0("sim_", i, ".sh"))
      script_con <- file(script_fname)
      cat(out_str, file = script_con)
      close(script_con)
      Sys.chmod(script_fname, mode = "0755")

      NULL
    }))


  list(counts = counts, size = size, rsem = rsem_tbl)
}

#' Counts to rsem-simulate-expression
#'
#'
#' Take a set of counts and make several scripts to simulate counts for each
#' sample using \code{rsem-simulate-expression} with the --deterministic flag
#' @param counts a matrix with counts from \code{simulate_counts}
#' @param info a \code{data.frame} output from \code{simulate_counts}
#' @param rsem_example an example rsem \code{data.frame} read which is used for
#' length and effective length
#' @param out_dir the base directory to output simulation scripts
generate_counts <- function(counts, info, rsem_example, out_dir = "gen_sim") {
  stopifnot( is(counts, "matrix") )
  stopifnot( is(info, "data.frame") )

  rsem_tbl <- lapply(1:ncol(counts),
    function(cnt)
    {
      to_rsem(info$target_id, cnt, info$length, base_df$eff_len)
    })

  counts <- data.frame(counts)
  colnames(counts) <- 1:ncol(exp_means)
  counts$target_id <- base_df$target_id

  dir.create(out_dir)
  write.table(base_df,
    file = file.path(out_dir, "de.txt"),
    sep = "\t",
    row.names = FALSE, quote = FALSE, col.names = TRUE)
  write.table(counts,
    file = file.path(out_dir, "counts.txt"),
    sep = "\t",
    row.names = FALSE, quote = FALSE, col.names = TRUE)

  invisible(lapply(seq_along(rsem_tbl),
    function(i)
    {
      dir.create(file.path(out_dir, i))
      out_fname <- file.path(out_dir, i, paste0(i, ".isoforms.results"))
      write.table(rsem_tbl[[i]],
        file = out_fname,
        sep = "\t", eol = "\n",
        row.names = FALSE, col.names = TRUE,
        quote = FALSE)
      out_str <- paste("#!/bin/bash",
        "",
        'if [ "$#" -ne 3 ]; then',
        ' echo "Requires two arguments RSEM_REFERENCE RSEM_MODEL_FILE"',
        ' exit 1',
        'fi',
        '',
        "RSEM_REFERENCE=$1",
        "MODEL_FILE=$2",
        "OUT_DIR=$3",
        sep = "\n")
      theta <- 0.0
      total_reads <- sum(rsem_tbl[[i]]$expected_count)
      out_str <- c(out_str, "",
        paste("rsem-simulate-reads", "${RSEM_REFERENCE}", "${MODEL_FILE}",
          paste0("${OUT_DIR}/", i, ".isoforms.results"),
          theta, total_reads, paste0("${OUT_DIR}/sim_", i), "--seed", i, "\n"))
      out_str <- c(out_str,
        "",
        "gzip *.fq", "\n")

      out_str <- paste(out_str, collapse = "\n")

      script_fname <- file.path(out_dir, i, paste0("sim_", i, ".sh"))
      script_con <- file(script_fname)
      cat(out_str, file = script_con)
      close(script_con)
      Sys.chmod(script_fname, mode = "0755")

      NULL
    }))


  list(counts = counts, size = size, rsem = rsem_tbl)
}


# Generate outline:
# - Generate fold changes as above
# - Output this directory
# - Call snakemake to make a bunch of reads
# - Quantify
# - Call DE
sim_res <- generate_counts(sim_fc, fc, out_dir = "de_4_4", n_sim = c(4, 4))
