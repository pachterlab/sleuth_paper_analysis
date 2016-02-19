library("BitSeq")

base_path <- "../results/3_3_1_1/sample_1"

bitseq_fnames <- file.path(base_path, 1:6, "bitseq", "bitseq.rpkm")

cond_a <- 1:3
cond_b <- 4:6
all_conditions <- list(bitseq_fnames[cond_a], bitseq_fnames[cond_b])

bitseq_res <- getDE(all_conditions, outPrefix = file.path(base_path, "bitseq_de"))
saveRDS(bitseq_res, file = "../results/3_3_1_1/sample_1/bitseq_de.rds")
