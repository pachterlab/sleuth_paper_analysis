sim_name <- 'gfr_3_3_20_42_2'
sim_number <- 1
benchmarks <- paste0("../benchmark/", sim_name, "/exp_", sim_number, "/", 1:6, "/BitSeq_cli/getExpression.json")

library('jsonlite')
timings <- fromJSON(benchmarks[1])
