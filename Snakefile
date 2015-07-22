ANNO = "annotation"
ANNO_PRE = "Homo_sapiens.GRCh38.rel80.cdna.all"

ANNO_FA = "{0}/{1}.fa".format(ANNO, ANNO_PRE)
BWT_IDX = "index/{0}".format( ANNO_PRE )

#RSEM_ANNO = "{ANNO}/Homo_sapiens.GRCh38.rel80.cdna.all_rsem/ref"

rule all:
    input:
        BWT_IDX + '.1.ebwt'

rule bowtie_idx:
    input:
        ANNO_FA
    output:
        expand(BWT_IDX + '.{i}.ebwt', i = range(1, 5)),
        expand(BWT_IDX + '.rev.{i}.ebwt', i = range(1, 3))
    shell:
        'bowtie-build '
        '--offrate 1 '
        '--seed 37 '
        #'--ftabchars 15 ' # use about a 4GB for lookup table
        '--ntoa '
        '{input} ' +
        BWT_IDX
