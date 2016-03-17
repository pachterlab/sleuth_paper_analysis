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

saveRDS(sgr, file = '../results/sgr.rds')
