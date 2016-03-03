source("benchmark_methods.R")
source("gene_common.R")

library("cowplot")
library("mamabear")

###
# simulation specific constants
###
sim_name <- 'isoform_3_3_5_1_1'
sim_name <- 'isoform_3_3_20_1_1'
sims_dir <- file.path('../sims', sim_name)
base_dir <- file.path('../sims', sim_name)
n_a <- 3
n_b <- 3
n <- n_a + n_b


###
# generalization to several simulations
###

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

# @param sim_name a simulation name such as 'isoform_3_3_20_1_1'
# @param which_sample which sample (replication) to load (an integer from 1 to N)
# @param method_filtering if \code{TRUE}, use the methods own filtering.
# Otherwise, use the filtering provided from sleuth.
# NOTE: cuffdiff uses a filter method internally and it cannot be changed
load_results <- function(sim_name, which_sample, method_filtering = FALSE,
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
  sir$so <- NULL

  gene_methods <- list(
    DESeq2 = runDESeq2,
    edgeR = runEdgeR,
    limmaVoom = runVoom
    # edgerRobust = runEdgeRRobust,
    # EBSeq = runEBSeq
    )
  isoform_results <- lapply(gene_methods,
    function(f) {
      f(isoform_cds, FALSE, method_filtering)
    })
  all_results <- c(Filter(is.data.frame, sir) , isoform_results)

  all_results
}

debugonce(load_results)
all_results <- lapply(1:5,
  function(i) {
    cat('Sample: ', i, '\n')
    load_results(sim_name, i)
  })

# library('BiocParallel')
# mcp <- MulticoreParam(workers = 6)
transcript_gene_mapping <- get_human_gene_names()

all_results_10 <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_results(sim_name, i, min_reads = 10)
  }, mc.cores = 10)

saveRDS(all_results_10, file = '../sims/isoform_3_3_20_1_1/tmp.rds')
all_results_10 <- readRDS('../sims/isoform_3_3_20_1_1/tmp.rds')

all_results_10_filter <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_results(sim_name, i, method_filtering = TRUE, min_reads = 10)
  }, mc.cores = 10)
saveRDS(all_results_10_filter, file = '../sims/isoform_3_3_20_1_1/all_results_10_filter.rds')
all_results_10_filter <- readRDS('../sims/isoform_3_3_20_1_1/all_results_10_filter.rds')

all_bench_10 <- lapply(seq_along(all_results_10),
  function(i) {
    res <- all_results_10[[i]]
    sim_info <- get_de_info(sim_name, i, transcript_gene_mapping)
    bench <- new_de_benchmark(res, names(res), sim_info$de_info)
    bench
  })

fdr_nde_plot(all_bench_10) +
  xlim(0, 3000) +
  ylim(0, 0.10) +
  theme_cowplot(25) +
  theme(legend.position = c(0.15, 0.80))

all_results_5 <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_results(sim_name, i, min_reads = 5)
  }, mc.cores = 10)

saveRDS(all_results_5, file = '../sims/isoform_3_3_20_1_1/5.rds')
all_results_5 <- readRDS('../sims/isoform_3_3_20_1_1/5.rds')

all_bench_10_filter <- lapply(seq_along(all_results_10_filter),
  function(i) {
    res <- all_results_10_filter[[i]]
    sim_info <- get_de_info(sim_name, i, transcript_gene_mapping)
    bench <- new_de_benchmark(res, names(res), sim_info$de_info)
    bench
  })

all_bench_10_filter <- lapply(seq_along(all_results_10_filter),
  function(i) {
    res <- all_results_10_filter[[i]]
    sim_info <- get_de_info(sim_name, i, transcript_gene_mapping)

    # bench <- new_de_benchmark(res, names(res), sim_info$de_info)
    bench <- compare_reference(res, sim_info$de_info,
      c('sleuth.wt', 'sleuth.lrt'))
    bench
  })

all_bench_10_filter_pairwise <- Map(function(m) {
    lapply(all_bench_10_filter, function(x) x[[m]])
  },
  c('DESeq2', 'edgeR', 'limmaVoom'))

p <- lapply(all_bench_10_filter_pairwise,
  function(x) {
    fdr_nde_plot(x) +
      xlim(0, 3000) +
      ylim(0, 0.10) +
      theme_cowplot(25) +
      theme(legend.position = c(0.15, 0.80))
  })

plot_grid(plotlist = p)

fdr_nde_plot(all_bench_10_filter) +
  xlim(0, 3000) +
  ylim(0, 0.10) +
  theme_cowplot(25) +
  theme(legend.position = c(0.15, 0.80))

###
# pairwise comparison filtering at 5
###

all_results_5_filter <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_results(sim_name, i, method_filtering = TRUE, min_reads = 5, min_prop = 0.3)
  }, mc.cores = 20)
saveRDS(all_results_5_filter, file = '../sims/isoform_3_3_20_1_1/all_results_5_filter.rds')

all_results_5_filter <- readRDS('../sims/isoform_3_3_20_1_1/all_results_5_filter.rds')
all_bench_5_filter <- lapply(seq_along(all_results_5_filter),
  function(i) {
    res <- all_results_5_filter[[i]]
    sim_info <- get_de_info(sim_name, i, transcript_gene_mapping)

    # bench <- new_de_benchmark(res, names(res), sim_info$de_info)
    bench <- compare_reference(res, sim_info$de_info,
      c('sleuth.wt', 'sleuth.lrt'))
    bench
  })

all_bench_5_filter_pairwise <- Map(function(m) {
    lapply(all_bench_5_filter, function(x) x[[m]])
  },
  c('DESeq2', 'edgeR', 'limmaVoom'))

p <- lapply(all_bench_5_filter_pairwise,
  function(x) {
    fdr_nde_plot(x) +
      xlim(0, 3000) +
      ylim(0, 0.10) +
      theme_cowplot(25) +
      theme(legend.position = c(0.15, 0.80))
  })

plot_grid(plotlist = p)
# all_results_10 <- lapply(1,
#   function(i) {
#     cat('Sample: ', i, '\n')
#     load_results(sim_name, i, min_reads = 10)
#   })

all_bench_5 <- lapply(seq_along(all_results_5),
  function(i) {
    res <- all_results_5[[i]]
    sim_info <- get_de_info(sim_name, i, transcript_gene_mapping)
    bench <- new_de_benchmark(res, names(res), sim_info$de_info)
    bench
  })

fdr_nde_plot(all_bench_5) +
  xlim(0, 3000) +
  ylim(0, 0.10) +
  theme_cowplot(25) +
  theme(legend.position = c(0.15, 0.80))

###
# end here for now
###

all_bench <- lapply(seq_along(all_results),
  function(i) {
    res <- all_results[[i]]
    sim_info <- get_de_info(sim_name, i, transcript_gene_mapping)
    bench <- new_de_benchmark(res, names(res), sim_info$de_info)
    bench
  })

debugonce(average_bench_fdr)

tmp <- filter(all_fdr, method == 'sleuth.lrt')

# figure out why getting NAs in the sd columns
debugonce(mamabear:::average_bench_fdr)
all_fdr <- mamabear:::average_bench_fdr(all_bench)

fdr_efdr_plot(all_bench)

fdr_nde_plot(all_bench) +
  xlim(0, 3000) +
  ylim(0, 0.10) +
  theme_cowplot(25) +
  theme(legend.position = c(0.15, 0.80))

ggplot(all_fdr, aes(nde, mean_tFDR)) +
  geom_line(aes(color = method, linetype = method)) +
  xlim(0, 3000) +
  ylim(0, 0.10)

fdr_nde_plot(all_bench[[1]], FALSE) +
  xlim(0, 3000) +
  ylim(0, 0.10) +
  theme_cowplot(25) +
  theme(legend.position = c(0.15, 0.80)) +
  xlab('number of isoforms called DE') +
  ylab('FDR')

###
# more constants
###

sim_info <- get_de_info(sim_name, 1, transcript_gene_mapping)
de_info <- sim_info$de_info
de_genes <- sim_info$de_genes

kal_dirs <- file.path(base_dir, "exp_1", 1:n, "kallisto")
sample_to_condition <- get_sample_to_condition(n_a, n_b, kal_dirs)

cr <- get_cuffdiff(file.path(base_dir, 'exp_1', 'results', 'cuffdiff'))

###
# try to estimate the best filtering
###
library('parallel')
minimum_reads <- c(0, 1, 3, 5, 10, 20)
filter_sleuth <- mclapply(minimum_reads,
  function(mr) {
    res <- run_sleuth(sample_to_condition, min_reads = mr)
    res$so <- NULL
    res
  },
  mc.cores = length(minimum_reads),
  mc.preschedule = FALSE)

saveRDS(filter_sleuth, file = '../results/filter_sleuth.rds')

filter_sleuth <- readRDS('../results/filter_sleuth.rds')

# trying to generalize the fdr


all_bench <- lapply(seq_along(filter_sleuth),
  function(i) {
    tmp <- filter_sleuth[[i]]$sleuth.lrt
    filter_sleuth[[i]]$sleuth.lrt <- tmp[complete.cases(tmp), ]

    tmp <- filter_sleuth[[i]]$sleuth.wt
    filter_sleuth[[i]]$sleuth.wt <- tmp[complete.cases(tmp), ]

    new_de_benchmark(filter_sleuth[[i]],
      paste0(names(filter_sleuth[[i]]), ' ', minimum_reads[i]),
      de_info)
  })

tmp <- mamabear:::calculate_fdr(all_bench[[1]])


fdr_plots <- lapply(all_bench,
  function(bench) {
    fdr_nde_plot(bench, FALSE) +
      xlim(0, 3000) +
      ylim(0, 0.10) +
      theme_cowplot(25) +
      theme(legend.position = c(0.15, 0.80)) +
      xlab('number of isoforms called DE') +
      ylab('FDR')
  })

plot_grid(plotlist = fdr_plots, nrow = 2)

##################################################
# testing

so <- sgr$so
rownames(so$sample_to_covariates) <- so$sample_to_variance$sample
sample_to_condition <- get_sample_to_condition(3, 3, 'hello')
tmp <- get_filtered_isoform_cds(so, so$sample_to_covariates)
tmp <- get_filtered_isoform_cds(so, sample_to_condition)




##################################################


###
# run the isoform analysis
###
sir <- run_sleuth(sample_to_condition, lift_genes = FALSE)

isoform_cds <- get_filtered_isoform_cds(sir$so, sample_to_condition)

gene_methods <- list(
  DESeq2 = runDESeq2,
  edgeR = runEdgeR,
  limmaVoom = runVoom
  # edgerRobust = runEdgeRRobust,
  # EBSeq = runEBSeq
  )
isoform_results <- lapply(gene_methods,
  function(f) {
    f(isoform_cds, FALSE)
  })

all_results <- c(Filter(is.data.frame, sir) , isoform_results,
  list(Cuffdiff2 = cr$isoform))

all_results <- c(Filter(is.data.frame, sir) , isoform_results)

library('colorspace')


de_bench <- new_de_benchmark(all_results, names(all_results), de_info,
  rainbow_hcl(length(all_results), c = 100))

fdr_nde_plot(de_bench, FALSE) +
  xlim(0, 3000) +
  ylim(0, 0.10) +
  theme(legend.position = c(0.1, 0.80)) +
  xlab('number of isoforms called DE') +
  ylab('FDR') +
  theme_cowplot(25)

fdr_tpr_plot(de_bench) +
  xlim(0, 0.25) +
  theme_cowplot(25)

fdr_efdr_plot(de_bench) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0, 0.25) +
  ylim(0, 0.25) +
  theme_cowplot(25)

# gene lifting
sgr <- run_sleuth(sample_to_condition, min_reads = 5, lift_genes = TRUE)

gene_counts <- get_gene_counts(sim_name)
gene_filter <- apply(gene_counts, 1, basic_filter)
gene_counts_filtered <- gene_counts[gene_filter, ]
gene_cds <- make_count_data_set(gene_counts_filtered, sample_to_condition)

gene_methods_count <- list(
  DESeq2 = runDESeq2,
  edgeR = runEdgeR,
  limmaVoom = runVoom
  # edgerRobust = runEdgeRRobust,
  # EBSeq = runEBSeq
  )
gene_results <- lapply(gene_methods,
  function(f) {
    f(gene_cds, FALSE)
  })

all_gene_results <- c(Filter(is.data.frame, sgr), gene_results)

de_gene_bench <- new_de_benchmark(all_gene_results, names(all_gene_results),
  de_genes)

fdr_nde_plot(de_gene_bench, FALSE) +
  xlim(0, 4200) +
  ylim(0, 0.10) +
  theme(legend.position = c(0.1, 0.80)) +
  xlab('number of genes called DE') +
  ylab('FDR') +
  theme_cowplot(25)

fdr_tpr_plot(de_gene_bench) +
  xlim(0, 0.25) +
  theme_cowplot(25)

fdr_efdr_plot(de_gene_bench) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0, 0.25) +
  ylim(0, 0.25) +
  theme_cowplot(25)
