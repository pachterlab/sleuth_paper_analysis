suppressWarnings(suppressMessages(library("dplyr")))

#' Prepare a DESeqDataSet for simulation
#'
#' A basic routine that filters and adjusts dispersion estimates so that we can
#' simulate from it
#'
#' @param dds a DESeqDataSet which has \code{dispGeneEst}
#' @param min_mean the minumum mean counts to consider for differential
#' expression
#' @param min_dispersion the minimum dispersion to use. Default pulled from
#' DESeq2paper package
#' @return a data.frame that with can simulate from. Adds column sim_filt to
#' denote whether or target should be considered for DE.
#' @export
prep_dds_sim <- function(dds, min_mean = 5, min_dispersion = 1e-6) {
  fit <- as.data.frame(GenomicRanges::mcols(dds))
  fit <- fit %>%
    mutate(target_id = names(GenomicRanges::rowData(dds)))

  # we want to keep dispersion estimates within a reasonable range.
  # range was extracted from DESeq2paper package
  fit <- fit %>%
    mutate(disp_filt = dispGeneEst > min_dispersion) %>%
    mutate(disp_filt = ifelse(is.na(disp_filt), FALSE, disp_filt))

  # for dispersion estimates that are hard to estimate (due to low counts),
  # simply use a fixed estimate based on those which pass the filter
  med_dispersion <- fit %>%
    dplyr::filter(disp_filt) %>%
    .$dispGeneEst %>%
    median

  # we now have dispersion estimates for every single transcript
  fit <- fit %>%
    mutate(disp_final = ifelse(disp_filt, dispGeneEst, med_dispersion))

  # we don't want to simulate DE from transcripts that are unlikely expressed
  # things, so let's omit transcripts that are zero in every condition
  fit <- fit %>%
    mutate(sim_filt = !allZero & baseMean > min_mean)

  # reorder so that matches RSEM targets
  fit[order(fit$target_id),]
}

#' @param sim_df a data.frame which contains columns baseMean, disp_final
#' @param X the design matrix
#' @param size_factors the scaling of each sample
#' @return a matrix of integer counts
make_sim <- function(sim_df, X, size_factors = rep(1, nrow(X))) {
  # number of targets
  n <- nrow(sim_df)
  p <- nrow(X)

  disp <- sim_df$disp_final
  beta_full <- cbind( log(sim_df$baseMean), sim_df$log_fc )
  mu <- exp( beta_full %*% t(X) )

  mu <- sweep(mu, 2, size_factors, `*`)

  matrix(rnbinom(n * p, mu = mu, size = 1 / disp), ncol = p)
}

# all size factors equal to 1
sf_mode_1 <- function(n) {
  rep.int(1, n)
}

# random size factors such that the geometric mean is 1
sf_mode_2 <- function(n) {
  res <- rep.int(1L, n)

  done <- FALSE

  while (all(res == 1L)) {
    nsamp <- 0
    while (nsamp < n) {
      if ( (n - nsamp) == 1 ) {
        nsamp <- nsamp + 1
        res[nsamp] <- 1
      } else {
        eq <- sample(c(TRUE, FALSE), 1)
        if (eq) {
          nsamp <- nsamp + 1
          res[nsamp] <- 1
        } else {
          nsamp <- nsamp + 1
          res[nsamp] <- 1 / 3
          nsamp <- nsamp + 1
          res[nsamp] <- 3
        }
      }
    }
  }


  res[sample.int(length(res))]
}

set.seed(3)

tmp <- sapply(1:100,
  function(x) {
    sf_mode_2(n = 6)
  }) %>% t
tmp[sapply(1:nrow(tmp), function(i) all(tmp[i, ] == 1))]

#' Simulate counts from a negative binomial
#'
#' Simulate counts from a negative binomial using a data.frame prepared from
#' prep_dds_sim. This model assumes that everything is independent.
#'
#' Counts are generated from a truncated normal such that the minimum fold change
#' is a reasonable amount
#'
#' @param prep_df result from running \code{prep_dds_sim}
#' @param n_sim the number of simulations to draw from this particular (random)
#' configuration
#' @param n_a number of samples from condition A
#' @param n_b number of samples from condition B
#' @param prop_de proportion of targets to be differentially expressed
#' @param seed an integer seed
#' @param log_fc_sd the standard deviation to use for the log fold change
#' @param constant_fc if not \code{NULL}, every differentially expressed
#' transcript will have this fold change
#' @param sim_function the functions used to generate the simulation
#' @param min_magnitude the minimum magnitude for the truncated normal
#' @return a named list with components (1) counts a list of matrices with
#' counts, (2) prep_df - a data.frame with information about the differential
#' expression, (3) condition a factor vector indicating the condition for each
#' sample
simulate_counts <- function(prep_df,
  n_sim = 1,
  n_a = 3L,
  n_b = 3L,
  prop_de = 0.2,
  seed = 37L,
  log_fc_sd = 1,
  constant_fc = NULL,
  sim_function = make_sim,
  size_factor_function = sf_mode_1,
  min_magnitude = log(1.5)) {
  stopifnot( is(prep_df, "data.frame") )

  set.seed(seed)
  # n_de <- as.integer(ceiling(nrow(prep_df) * prop_de))
  # n_null <- as.integer(nrow(prep_df) - n_de)
  n_pass_filt <- sum(prep_df$sim_filt)
  n_de <- as.integer(ceiling(n_pass_filt * prop_de))
  n_null <- as.integer(n_pass_filt - n_de)

  n_samples <- n_a + n_b

  condition <- factor(c(rep("A", n_a), rep("B", n_b)))
  X <- model.matrix(~ condition)

  sim <- lapply(1:n_sim,
    function(i) {
      # this needs to be re factored so that every different n_sim has different
      # isoforms differentially expressed
      row_idx_de <- sample(which(prep_df$sim_filt), n_de, replace = FALSE)

      # choose which will be DE
      prep_df <- dplyr::mutate(prep_df, is_de = FALSE)

      prep_df$is_de[row_idx_de] <- TRUE

      # assign log fold change using a normal centered at zero
      # log_fc <- rnorm(n_de, mean = 0, sd = log_fc_sd)
      log_fc <- NULL
      if (is.null(constant_fc)) {
        log_fc <- truncated_normal(n_de, min_magnitude = min_magnitude,
          sd = log_fc_sd)
      } else {
        log_fc <- rep.int(constant_fc, n_de)
        log_fc <- sample(c(-1, 1), n_de, replace = TRUE) * log_fc
      }
      prep_df <- dplyr::mutate(prep_df, log_fc = 0)
      prep_df$log_fc[row_idx_de] <- log_fc

      size_factors <- size_factor_function(n_samples)
      s <- sim_function(prep_df, X, size_factors = size_factors)
      colnames(s)
      colnames(s) <- paste0(condition, c(1:n_a, 1:n_b))
      rownames(s) <- prep_df$target_id

      list(counts = s, info = prep_df, condition = condition,
        size_factors = size_factors)
    })

  # list(counts = sim, info = prep_df, condition = condition,
  #   size_factors = size_factors)
  sim
}

simulate_gene_counts <- function(prep_df,
  target_mapping,
  gene_label,
  n_sim = 1,
  n_a = 3,
  n_b = 3,
  prop_de = 0.2,
  seed = 37L,
  log_fc = 1,
  log_fc_sd = NA,
  sim_function = make_sim,
  # size_factors = rep(1, n_a + n_b),
  size_factor_function = sf_mode_1,
  generate_fold_changes = gene_fold_change
  ) {
  stopifnot( is(prep_df, "data.frame") )

  set.seed(seed)

  # a gene passes the filter is at least 1 transcript passes the filter
  prep_df <- left_join(prep_df, target_mapping, by = "target_id")
  genes_pass_filter <- prep_df %>%
    group_by_(gene_label) %>%
    summarize(gene_filter = any(sim_filt))
  prep_df <- left_join(prep_df, genes_pass_filter, by = gene_label)

  # "gene_names" refers to the genes that pass the filter
  gene_names <- dplyr::filter(prep_df, gene_filter)[[gene_label]] %>% unique
  n_genes <- length(gene_names)
  n_de <- ceiling(prop_de * n_genes)

  n_samples <- n_a + n_b

  condition <- factor(c(rep("A", n_a), rep("B", n_b)))
  X <- model.matrix(~ condition)

  # prep_df_init is a template that should be copied at the beginning of every
  # iteration
  prep_df_init <- prep_df

  sim <- lapply(1:n_sim,
    function(i) {
      prep_df <- prep_df_init

      # select the genes to be differentially expressed and choose their fold change
      which_de <- sample(gene_names, n_de, replace = FALSE)

      # this will initialize everything to `log_fc`, but will get overwritten
      # `generate_fold_changes()` below
      which_de <- data.frame(gene_name = which_de, log_fc = log_fc, is_de = TRUE,
        stringsAsFactors = FALSE)
      data.table::setnames(which_de, "gene_name", gene_label)

      prep_df <- dplyr::left_join(prep_df, which_de, by = gene_label)
      # prep_df <- dplyr::mutate(prep_df,
      #   log_fc = ifelse(is.na(log_fc), 0, log_fc),
      #   is_de = ifelse(is.na(is_de), FALSE, is_de)
      # )
      prep_df <- dplyr::mutate(prep_df,
      is_de = ifelse(is.na(is_de), FALSE, is_de)
      )

      # instead of complicating our lives with another join, generate fold changes for everything
      prep_df <- generate_fold_changes(prep_df)

      size_factors <- size_factor_function(n_samples)
      message(paste0("size factors: ", paste(size_factors, collapse = ' ')))

      # clean out the old fold changes
      prep_df <- dplyr::mutate(prep_df,
        log_fc = ifelse(is_de, log_fc, 0))
      s <- sim_function(prep_df, X, size_factors = size_factors)
      colnames(s)
      colnames(s) <- paste0(condition, c(1:n_a, 1:n_b))
      rownames(s) <- prep_df$target_id
      list(counts = s, info = prep_df, condition = condition,
        size_factors = size_factors)
    })

  #
  # sim <- lapply(1:n_sim,
  #   function(i) {
  #
  #   })
  #
  # list(counts = sim, info = prep_df, condition = condition,
  #   size_factors = size_factors)
  sim
}
# debugonce(simulate_counts)
# debugonce(make_sim)

gene_fold_change <- function(df) {
  df <- dplyr::group_by(df, ens_gene)
  dplyr::mutate(df, log_fc = gfc_same_direction(length(ens_gene)))
}

#' Generate gene fold counts in the same direction
#'
#' This function takes a number of transcripts and generates log fold changes in the same direction
#' @param n_transcripts the number of transcripts at this locus
#' @param direction the direction/scaling to use. If `NULL`, then a random direction will be chosen (-1, 1)
#' @param min_magnitude the minimum magnitude to truncate the normal distribution
#' @return generate \code{n_transcripts} samples using \code{truncated_normal}
gfc_same_direction <- function(n_transcripts, direction = NULL,
  min_magnitude = log(1.5)) {

  if (is.null(direction)) {
    direction <- sample(c(-1, 1), 1)
  }

  s <- truncated_normal(n_transcripts, min_magnitude)
  s <- abs(s) * direction

  s
}

# gfc_single_transcript <- function(n_transcripts, direction = NULL,
#   min_magnitude = log(1.5)) {
#
#   if (is.null(direction)) {
#     direction <- sample(c(-1, 1), 1)
#   }
#
#   s <-
# }
#' simulate from a truncated normal
#'
#' simulate from a truncated normal which is missing the center part of the mass.
#' simulation is performed using rejection sampling.
#'
#' @param n the number of samples
#' @param min_magnitude the minimum magnitude of a sample to accept
#' @param ... additional parameters sent to \code{rnorm}
#' @return \code{n} samples that all have magnitude greater than \code{min_magnitude}
truncated_normal <- function(n, min_magnitude, ...) {
  samples <- numeric(n)
  i <- 1
  while (i <= n) {
    current_sample <- rnorm(1, ...)
    if (abs(current_sample) >= min_magnitude) {
      samples[i] <- current_sample
      i <- i + 1
    }
  }

  samples
}


# let's get the deepest sample from female Finns to use as our RSEM example for
# effective lengths and lengths.
#
# XXX: we need to quantify this using RSEM now...

# deepest_fin_female <- names( which.max( apply(fin_females, 2, sum) ) )
#
# sim1 <- simulate_counts(prep_fin, n_sim = 10, seed = 42)
#
# mean_a <- apply(sim1$counts[[1]][, sim1$condition %in% "A"], 1, mean)
# mean_b <- apply(sim1$counts[[1]][, sim1$condition %in% "B"], 1, mean)
#
# lfc <- log((mean_b + 0.5)/ (mean_a + 0.5))
# lfc[is.nan(lfc)] <- 0
#
# cor(lfc, sim1$prep_df$log_fc)
#
# cor(lfc[which(sim1$prep_df$is_de)], filter(sim1$prep_df, is_de)$log_fc)
#
# plot(lfc[which(sim1$prep_df$is_de)], filter(sim1$prep_df, is_de)$log_fc)
