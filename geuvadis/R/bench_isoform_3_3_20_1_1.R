source("benchmark_methods.R")
source("gene_common.R")

library("parallel")
library("cowplot")
library("mamabear")

sim_name <- 'isoform_3_3_20_1_1'

n_cpu <- 20

reference_methods <- c('sleuth.wt', 'sleuth.lrt')

transcript_gene_mapping <- get_human_gene_names()

# compute all of the results using each method independently
all_results_5_filter <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_isoform_results(sim_name, i, method_filtering = FALSE, min_reads = 5,
      min_prop = 0.3)
  }, mc.cores = n_cpu)

# compare each to a reference, in this case sleuth.
# afterwards, `all_bench_5_filter` contains a list of 20 lists,
# each of which contains comparisons of a method to the reference
all_bench_5_filter <- mclapply(seq_along(all_results_5_filter),
  function(i) {
    res <- all_results_5_filter[[i]]
    sim_info <- get_de_info(sim_name, i, transcript_gene_mapping)

    bench <- compare_reference(res, sim_info$de_info, reference_methods)

    bench
  }, mc.cores = n_cpu)

other_methods <- names(all_bench_5_filter[[1]])

# arrange the results so that all of the replicates are under the method name
# e.g. all_bench_5_filter_pairwise[['Cuffdiff2']] contains N samples of the same
# simulation
all_bench_5_filter_pairwise <- Map(
  function(m) {
    lapply(all_bench_5_filter, function(x) x[[m]])
  }, other_methods)

tmp_file_name <- paste0('../sims/', sim_name, '5_filter_pairwise.rds')

saveRDS(all_bench_5_filter_pairwise, file = tmp_file_name)

all_bench_5_filter_pairwise <- readRDS(tmp_file_name)

tmp <- lapply(all_bench_5_filter_pairwise, function(x) {
    lapply(x, function(y) {
      a <- lapply(y, function(z) {
        if (is(z, 'data.frame')) {
          z[complete.cases(z), ]
        } else {
          z
        }
      })
      class(a) <- 'de_benchmark'
      a
    })
  })

# we can now make aggregate plots like below
p <- lapply(all_bench_5_filter_pairwise,
  function(x) {
    fdr_nde_plot(x) +
      xlim(0, 3000) +
      ylim(0, 0.10) +
      theme_cowplot(25) +
      theme(legend.position = c(0.15, 0.80))
  })

plot_grid(plotlist = p)

# we can now make aggregate plots like below
p <- lapply(tmp,
  function(x) {
    fdr_nde_plot(x) +
      xlim(0, 3000) +
      ylim(0, 0.10) +
      theme_cowplot(25) +
      theme(legend.position = c(0.15, 0.80))
  })

plot_grid(plotlist = p)

###
# trying to optimize limma
###

# TODO: use edgeR filter with limma
i <- 1
l_test <- load_isoform_results(sim_name, i, method_filtering = TRUE, min_reads = 5,
    min_prop = 0.3)
sim_info <- get_de_info(sim_name, i, transcript_gene_mapping)

l_bench <- new_de_benchmark(l_test, names(l_test), sim_info$de_info)
fdr_nde_plot(l_bench) +
  xlim(0, 3000) +
  ylim(0, 0.10) +
  theme_cowplot(25) +
  theme(legend.position = c(0.15, 0.80))

# TODO: use no filter but remove all zeros from limma

i <- 1
sim_info <- get_de_info(sim_name, i, transcript_gene_mapping)

l_test <- load_isoform_results_intersect(sim_name, i, method_filtering = TRUE, min_reads = 5,
    min_prop = 0.3)
l_test <- load_isoform_results_intersect(sim_name, i, method_filtering = TRUE)

l_bench <- new_de_benchmark(l_test$limma, names(l_test$limma), sim_info$de_info)
fdr_nde_plot(l_bench) +
  xlim(0, 3000) +
  ylim(0, 0.10) +
  theme_cowplot(25) +
  theme(legend.position = c(0.15, 0.80))

###
# doing the pairwise comparison
###

l_pairwise <- load_isoform_results_intersect(sim_name, i, 'limmaVoom', limma_filter_and_run)

l_pairwise_bench <- new_de_benchmark(l_pairwise, names(l_pairwise), sim_info$de_info)

fdr_nde_plot(l_pairwise_bench) +
  xlim(0, 3000) +
  ylim(0, 0.10) +
  theme_cowplot(25) +
  theme(legend.position = c(0.15, 0.80))

debugonce(DESeq2_filter_and_run)
d_pairwise <- load_isoform_results_intersect(sim_name, i, 'DESeq2', DESeq2_filter_and_run)

d_pairwise_bench <- new_de_benchmark(d_pairwise, names(d_pairwise), sim_info$de_info)

a <- fdr_nde_plot(d_pairwise_bench) +
  xlim(0, 3000) +
  ylim(0, 0.10) +
  theme_cowplot(25) +
  theme(legend.position = c(0.15, 0.80))

###
# see which version of DESeq does better
###

d2_pairwise <- load_isoform_results_intersect(sim_name, i, 'DESeq2', DESeq2_filter_and_run_intersect)

d2_pairwise_bench <- new_de_benchmark(d2_pairwise, names(d2_pairwise), sim_info$de_info)

b <- fdr_nde_plot(d2_pairwise_bench) +
  xlim(0, 3000) +
  ylim(0, 0.10) +
  theme_cowplot(25) +
  theme(legend.position = c(0.15, 0.80))

plot_grid(a, b, ncol = 1)

###
# new refactoring
###

all_results <- list()

all_results$limmaVoom <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_isoform_results_intersect(sim_name, i, 'limmaVoom',
      limma_filter_and_run)
  }, mc.cores = n_cpu)

all_results$DESeq2 <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_isoform_results_intersect(sim_name, i, 'DESeq2',
      DESeq2_filter_and_run_intersect)
  }, mc.cores = n_cpu)

all_results$Cuffdiff2 <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_isoform_results_intersect(sim_name, i, 'Cuffdiff2',
      cuffdiff_filter_and_run)
  }, mc.cores = n_cpu)


all_results$edgeR <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_isoform_results_intersect(sim_name, i, 'edgeR',
      edgeR_filter_and_run)
  }, mc.cores = n_cpu)

all_results$EBSeq <- mclapply(1:20,
  function(i) {
    cat('Sample: ', i, '\n')
    load_isoform_results_intersect(sim_name, i, 'EBSeq',
      EBSeq_isoform_filter_and_run)
  }, mc.cores = n_cpu)


all_benchmarks <- lapply(all_results,
  function(result) {
    mclapply(seq_along(result),
      function(i) {
        x <- result[[i]]
        sim_info <- get_de_info(sim_name, i, transcript_gene_mapping)
        new_de_benchmark(x, names(x), sim_info$de_info)
      }, mc.cores = n_cpu)
  })

all_nde_plot <- mclapply(all_benchmarks,
  function(bench) {
    fdr_nde_plot(bench) +
      xlim(0, 3000) +
      ylim(0, 0.10) +
      theme_cowplot(25) +
      theme(legend.position = c(0.15, 0.80))
  })

saveRDS(all_nde_plot, file = '../sims/temporary_all_nde_plot.rds')

all_nde_plot <- readRDS('../sims/temporary_all_nde_plot.rds')

# put the isoform information for EBSeq in global variables
NG_LIST <- GetNg(transcript_gene_mapping$target_id,
  transcript_gene_mapping$ens_gene)

plot_grid(plotlist = all_nde_plot)


# debugonce(cuffdiff_filter_and_run)
tmp <- load_isoform_results_intersect(sim_name, 1, 'edgeR',
      edgeR_filter_and_run)

tmp <- get_sample_to_condition(3, 3, rep("a", 6))
temporary_2 <- so$bs_summary$obs_counts
colnames(temporary_2) <- tmp$sample
cds <- make_count_data_set(temporary_2, tmp)

# debugonce(EBSeq_isoform)
tmp_EBSeq <- load_isoform_results_intersect(sim_name, 1, 'EBSeq',
  EBSeq_isoform_filter_and_run)

sim_info <- get_de_info(sim_name, i, transcript_gene_mapping)

de_temporary <- new_de_benchmark(tmp_EBSeq, names(tmp_EBSeq), sim_info$de_info)

fdr_nde_plot(de_temporary) +
  xlim(0, 3000) +
  ylim(0, 0.10) +
  theme_cowplot(25) +
  theme(legend.position = c(0.15, 0.80))
