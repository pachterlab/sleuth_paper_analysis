args <- commandArgs(trailingOnly = TRUE)

if ( length(args) != 4 ) {
    stop("Usage: Rscript run_bitseq.R varname bam-path")
}

varname <- args[1]
out_prefix <- args[2]
bam <- args[3]
fasta <- args[4]

library("BitSeq")

x <- getExpression(bam,
    fasta,
    outPrefix = out_prefix,
    log = TRUE,
    seed = 47
    )

#assign(varname, x)
#rm(x)
#gc()
#save(list = varname, file = out_fname)
