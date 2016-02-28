source('load_info.R')
source('../../geuvadis/R/benchmark_methods.R')

library('sleuth')

info <- dplyr::mutate(info,
  condition = ifelse(condition == 'scramble', 'A', 'B'))
sir <- run_sleuth(info, max_bootstrap = 100, gene_mode = 'lift')
