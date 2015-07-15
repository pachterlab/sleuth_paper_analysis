library("dplyr")

#' Output a valid RSEM data.frame
#'
#'
to_rsem <- function(target_id, counts, len, eff_len) {
  count_df <- data.frame(transcript_id = target_id, gene_id = target_id,
    length = len, effective_length = eff_len,
    expected_count = counts,
    TPM = sleuth::counts_to_tpm(counts, eff_len),
    FPKM = sleuth::counts_to_fpkm(counts, eff_len),
    IsoPct = ifelse(counts > 0, 100.00, 0.0),
    stringsAsFactors = FALSE)

  discard_counts_eff_len <- count_df$effective_length < .Machine$double.eps
  cat("Need to discard ", sum(discard_counts_eff_len),
    " transcripts from simulation due to small counts\n")
  cat("Total counts to discard: ",
    sum(count_df$expected_count[discard_counts_eff_len]), "\n")

  count_df$expected_count[discard_counts_eff_len] <- 0
  count_df$IsoPct[discard_counts_eff_len] <- 0.0

  count_df
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
counts_to_simulation <- function(counts, info, rsem_example, out_dir = "gen_sim") {
  stopifnot( is(counts, "matrix") )
  stopifnot( is(info, "data.frame") )

  # TODO: check that rsem_example$target_id and info$target_id are the same and
  # in order.
  if ( nrow(rsem_example) != nrow(info) ) {
    stop("Inconsistent number of rows between 'info' and 'rsem_example'")
  }

  if ( !all(rsem_example$transcript_id == info$target_id) ) {
    stop("transcript_ids are in different order.\n")
  }

  # TODO: check and see if we assign DE to things that have effective length zero

  rsem_tbl <- lapply(1:ncol(counts),
    function(cnt) {
      to_rsem(info$target_id, counts[,cnt], rsem_example$length,
        rsem_example$effective_length)
    })

  which_lt_eps <- info$effective_length < .Machine$double.eps
  if ( any(info$is_de[which_lt_eps]) ) {
    cat("***Warning: ", sum(info$is_de[which_lt_eps]),
      " rows that are DE with effective length zero.\n")
    print(info[info$is_de & which_lt_eps,])
  }

  dir.create(out_dir, showWarnings = FALSE)

  counts_handle <- gzfile(file.path(out_dir, "counts.tsv.gz"), open = "w")
  write.table(counts,
    file = counts_handle,
    sep = "\t",
    row.names = TRUE, quote = FALSE, col.names = TRUE)
  close(counts_handle)

  invisible(lapply(seq_along(rsem_tbl),
    function(i) {
      dir.create(file.path(out_dir, i))
      out_fname <- file.path(out_dir, i, paste0(i, ".isoforms.results"))
      write.table(rsem_tbl[[i]],
        file = out_fname,
        sep = "\t", eol = "\n",
        row.names = FALSE, col.names = TRUE,
        quote = FALSE)
      out_str <- paste('#!/bin/bash',
        '',
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
        paste(
          "rsem-simulate-reads",
          "${RSEM_REFERENCE}",
          "${MODEL_FILE}",
          paste0("${OUT_DIR}/", i, ".isoforms.results"),
          theta,
          total_reads,
          paste0("${OUT_DIR}/sim_", i),
          "--seed", i,
          "--deterministic",
          "\n")
        )
      out_str <- c(out_str,
        "",
        "gzip ${OUT_DIR}/*.fq", "\n")

      out_str <- paste(out_str, collapse = "\n")

      script_fname <- file.path(out_dir, i, paste0("sim_", i, ".sh"))
      script_con <- file(script_fname)
      cat(out_str, file = script_con)
      close(script_con)
      Sys.chmod(script_fname, mode = "0755")

      NULL
    }))

  list(counts = counts, rsem = rsem_tbl)
}
