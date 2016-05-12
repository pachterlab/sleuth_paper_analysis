args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop('Usage: Rscript BitSeq_getExpression.R BAM TRANSCRIPTOME_FA OUTPUT_PREFIX')
}

alignment_file <- '../sims/gfr_3_3_20_42_2/exp_1/3/bowtie.bam'

alignment_file <- '../sims/gfr_3_3_20_42_2/exp_1/5/bowtie.bam'
transcriptome <- '../../annotation/Homo_sapiens.GRCh38.cdna.all.rel80.fa'
out <- '../sims/gfr_3_3_20_42_2/exp_1/5/BitSeq'

alignment_file <- args[1]
transcriptome <- args[2]
out <- args[3]

seed <- 42

suppressPackageStartupMessages(library('BitSeq'))

system.time(res <- getExpression(alignment_file, transcriptome, out,
  seed = seed))

###
# attempting to parallelize
# IMPORTANT: the code below does not run. it also segfaults (see the Snakefile)
###

# library('parallel')
#
# options(mc.cores = 6)
# alignment_base <- '../sims/gfr_3_3_20_42_2/exp_1'
# transcriptome <- '../../annotation/Homo_sapiens.GRCh38.cdna.all.rel80.fa'
# out <- '../sims/gfr_3_3_20_42_2/exp_1/5/BitSeq'
#
# all_quant <- mclapply(1:6,
#   function(i) {
#     alignment_file <- file.path(alignment_base, i, 'bowtie.bam')
#     print(alignment_file)
#     out_base <- file.path(alignment_base, i, 'BitSeq')
#     getExpression(alignment_file, transcriptome, out, seed = seed)
#   })
