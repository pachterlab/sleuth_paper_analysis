source("https://bioconductor.org/biocLite.R")
biocLite("tximport")

library('DESeq2')
library('tximport')


###
# import transcripts
###

counts <- tximport(file.path(info$path, 'abundance.tsv'),
  type = 'kallisto', tx2gene = transcript_gene_mapping)
