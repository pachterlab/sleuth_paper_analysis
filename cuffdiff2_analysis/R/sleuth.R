source('load_info.R')
source('../../geuvadis/R/benchmark_methods.R')
source('../../simulation_core/R/simulate_de.R')

library('sleuth')

info <- dplyr::mutate(info,
  condition = ifelse(condition == 'scramble', 'A', 'B'))

###
# get the union of sleuth results using lifting
###

sgr <- run_sleuth(info, max_bootstrap = 100, gene_mode = 'lift')
sr <- dplyr::full_join(
  dplyr::filter(sgr$sleuth.lrt, qval < 0.05),
  dplyr::filter(sgr$sleuth.wt, qval < 0.05),
  by = 'target_id')

sr <- dplyr::select(ungroup(sr),
  -starts_with("ens_gene"),
  -starts_with('num_transcripts'),
  -starts_with('all_target_ids')
  )

###
# merge the results
###

all_results <- sr
all_results <- dplyr::select(all_results, ens_gene = target_id)

###
# get the fold change
###

significant_genes <- dplyr::inner_join(
  transcript_gene_mapping,
  all_results,
  by = 'ens_gene')

# get the normalized counts for everything (even the things that we filtered out)
counts <- kallisto_table(sgr$so)

###
# this is so we can match the ranks later. we don't want something with incredibly
# high counts to get a very large fold change from something that had very low counts
###
count_ranks <- dplyr::group_by(counts, target_id)
count_ranks <- dplyr::summarize(count_ranks, est_counts = mean(est_counts))

count_ranks <- dplyr::inner_join(count_ranks, transcript_gene_mapping,
  by = 'target_id')
count_ranks <- dplyr::group_by(count_ranks, ens_gene)
count_ranks <- dplyr::mutate(count_ranks, rank = rank(est_counts, ties = 'first'))
count_ranks <- dplyr::ungroup(count_ranks)
count_ranks <- dplyr::select(count_ranks, -ens_gene)

fc <- dplyr::summarize(group_by(counts, target_id, condition),
  average_counts = mean(est_counts))
fc <- dplyr::summarize(fc,
  fold_change =
    average_counts[condition == 'B'] / average_counts[condition == 'A'])
fc <- dplyr::inner_join(fc, significant_genes, by = 'target_id')

# note that in R
# 0 / 0 == NaN
# 3 / 0 == Inf
fc <- dplyr::mutate(fc, fold_change = ifelse(is.nan(fold_change),
  1, fold_change))

# If it is infinite, give it a relatively large fold_change
large_fold_change <- quantile(fc$fold_change[is.finite(fc$fold_change)],
  probs = 0.75)

# dplyr::mutate(data.frame(valid = c(TRUE, FALSE, TRUE, TRUE), value = c(NA, 3, NA, NA)),
#   new_value = ifelse(valid, truncated_normal(1, large_fold_change), value))
n_infinite <- sum(is.infinite(fc$fold_change))
infinite_replacement <- abs(truncated_normal(n_infinite, large_fold_change))

fc <- dplyr::mutate(fc, is_infinite = is.infinite(fold_change))

# TODO: verify that things are actually getting replaced correctly
fc <- dplyr::mutate(fc, fold_change = ifelse(is.infinite(fold_change),
  infinite_replacement, fold_change))

n_large <- sum(fc$fold_change > 5)
large_replacement <- abs(truncated_normal(n_large, large_fold_change))

fc <- dplyr::mutate(fc, fold_change = ifelse(fold_change > 5, large_replacement,
  fold_change))

fc <- dplyr::inner_join(fc,
  count_ranks, by = 'target_id')

fc <- dplyr::mutate(group_by(fc, ens_gene),
  n_transcripts = length(ens_gene))

small_fold_change <- quantile(fc$fold_change[fc$fold_change > 0], probs = 0.05)
fc <- dplyr::mutate(fc, fold_change = pmax(fold_change, 0.001))

fc <- dplyr::mutate(fc, log_fc = log(fold_change))

saveRDS(fc, '../results/fc.rds')

# the distribution of transcripts per gene that we are looking at
n_transcripts_valid <- dplyr::distinct(fc, ens_gene) %>%
  group_by(n_transcripts) %>%
  summarize(length(n_transcripts)) %>%
  print(., n = nrow(.))
n_transcripts_valid <- n_transcripts_valid[['n_transcripts']]
