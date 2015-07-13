to_rsem <- function(target_id, counts, len, eff_len) {
  data.frame(transcript_id = target_id, gene_id = target_id,
    length = len, effective_length = eff_len,
    expected_count = counts,
    TPM = counts_to_tpm(counts, eff_len),
    FPKM = counts_to_fpkm(counts, eff_len),
    IsoPct = ifelse(counts > 0, 100.00, 0.0),
    stringsAsFactors = FALSE)
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
