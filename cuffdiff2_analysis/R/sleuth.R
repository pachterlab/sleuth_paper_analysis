source('load_info.R')
source('../../geuvadis/R/benchmark_methods.R')

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

# debugonce(sleuth_to_matrix)
# counts <- sleuth_to_matrix(sgr$so, 'obs_norm', 'est_counts')$data
counts <- kallisto_table(sgr$so)

fc <- dplyr::summarize(group_by(counts, target_id, condition),
  average_counts = mean(est_counts))
fc <- dplyr::summarize(fc,
  fold_change =
    average_counts[condition == 'B'] / average_counts[condition == 'A'])
fc <- dplyr::inner_join(fc, significant_genes, by = 'target_id')

# remove genes with only 1 transcript
fc <- dplyr::mutate(group_by(fc, ens_gene),
  n_transcripts = length(ens_gene))
fc <- dplyr::filter(fc, n_transcripts != 1)

# deal with infinite values
