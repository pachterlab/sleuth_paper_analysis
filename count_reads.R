simulation_files <- system("find . -name 'sims.rds'", intern = TRUE)

count_reads_from_simulation <- function(sim_file) {
  sim <- readRDS(sim_file)

  all_counts <- sapply(sim,
    function(s) {
      sum(s$counts)
    })

  sum(all_counts)
}

reads_per_simulation <- sapply(simulation_files, count_reads_from_simulation)
sum(reads_per_simulation)
# 13752129731
# 13,752,129,731

# make a table

sum_columns <- function(sim_file) {
  sim <- readRDS(sim_file)

  all_counts <- sapply(sim,
    function(s) {
      apply(s$counts, 2, sum)
    })

  t(all_counts)
}

sum_columns(simulation_files[1])

total_reads <- lapply(simulation_files, sum_columns)

simulation_labels <- lapply(simulation_files,
  function(x) {
    res <- sub('./geuvadis/sims/', '', x)
    sub('/sims.rds', '', res)
  })

names(total_reads) <- simulation_labels

lapply(seq_along(total_reads), function(i) {
  write.csv(total_reads[[i]],
    file = file.path('results', paste0(simulation_labels[[i]], '.csv')),
    quote = FALSE)
  invisible(NULL)
})
