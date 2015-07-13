# load DF prep_fin which is ready to simulate from
load("../results/prep_fin.RData", verbose = TRUE)

#' @param sim_df a data.frame which contains columns baseMean, disp_final#' @param X the design matrix
#' @param size_factors the scaling of each sample
#' @return a matrix of integer counts
make_sim <- function(sim_df, X, size_factors = c(1, nrow(X))) {
  # number of targets
  n <- nrow(sim_df)
  p <- nrow(X)

  disp <- sim_df$disp_final
  beta_full <- cbind( log(sim_df$baseMean), sim_df$log_fc )
  mu <- exp( beta_full %*% t(X) )

  matrix(rnbinom(n * p, mu = mu, size = 1 / disp), ncol = p)
}

#' Simulate counts from a negative binomial
#'
#' Simulate counts from a negative binomial using a data.frame prepared from
#' prep_dds_sim
#'
#' @param prep_df result from running \code{prep_dds_sim}
#' @param n_sim the number of simulations to draw from this particular (random)
#' configuration
#' @param n_a number of samples from condition A
#' @param n_b number of samples from condition B
#' @param prop_de proportion of targets to be differentially expressed
#' @param seed an integer seed
#' @param log_fc_sd the standard deviation to use for the log fold change
#' @return a named list with components (1) counts a list of matrices with
#' counts, (2) prep_df - a data.frame with information about the differential
#' expression, (3) condition a factor vector indicating the condition for each
#' sample
simulate_counts <- function(prep_df, n_sim = 1, n_a = 3L, n_b = 3L, prop_de = 0.2,
  seed = 37L, log_fc_sd = 1) {
  stopifnot( is(prep_df, "data.frame") )

  set.seed(seed)
  n_de <- as.integer(ceiling(nrow(prep_df) * prop_de))
  n_null <- as.integer(nrow(prep_df) - n_de)

  condition <- factor(c(rep("A", n_a), rep("B", n_b)))
  X <- model.matrix(~ condition)

  row_idx_de <- sample(which(prep_df$sim_filt), n_de, replace = TRUE)

  # choose which will be DE
  prep_df <- prep_df %>%
    mutate(is_de = FALSE)

  prep_df$is_de[row_idx_de] <- TRUE

  # assign log fold change using a normal centered at zero
  log_fc <- rnorm(n_de, mean = 0, sd = log_fc_sd)
  prep_df <- prep_df %>%
    mutate(log_fc = 0)
  prep_df$log_fc[row_idx_de] <- log_fc

  sim <- lapply(1:n_sim,
    function(i) {
      s <- make_sim(prep_df, X)
      colnames(s)
      colnames(s) <- paste0(condition, c(1:n_a, 1:n_b))
      rownames(s) <- prep_df$target_id
      s
    })

  list(counts = sim, prep_df = prep_df, condition = condition,
    size_factors = rep(1, nrow(X)))
}

debugonce(simulate_counts)
debugonce(make_sim)

sim1 <- simulate_counts(prep_fin, n_sim = 10, seed = 42)

mean_a <- apply(sim1$counts[[1]][, sim1$condition %in% "A"], 1, mean)
mean_b <- apply(sim1$counts[[1]][, sim1$condition %in% "B"], 1, mean)

lfc <- log((mean_b + 0.5)/ (mean_a + 0.5))
lfc[is.nan(lfc)] <- 0

cor(lfc, sim1$prep_df$log_fc)

cor(lfc[which(sim1$prep_df$is_de)], filter(sim1$prep_df, is_de)$log_fc)

plot(lfc[which(sim1$prep_df$is_de)], filter(sim1$prep_df, is_de)$log_fc)
