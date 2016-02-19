# BASE = '/home/hjp/sleuth_paper_analysis'
# BASE = '/Users/hjp/analysis/sleuth_paper'
# TODO: make this a absolute path before submission
from os.path import expanduser
HOME = expanduser('~')
BASE = HOME + '/sleuth_paper_analysis'

# software
KALLISTO = BASE + '/software/kallisto_linux/kallisto'
RSEM_SIMULATE = BASE + '/software/rsem_simulate/rsem-simulate-reads'

# functions
def source_r(base, fname):
    return 'Rscript --vanilla --default-packages=methods,stats,utils -e \'setwd("{0}")\' -e \'source("{1}")\''.format(base, fname)

# annotations
TRANSCRIPTOME_NAME = 'Homo_sapiens.GRCh38.cdna.all.rel80'
TRANSCRIPTOME_FA = BASE + '/annotation/' + TRANSCRIPTOME_NAME + '.fa.gz'
TRANSCRIPTOME_GTF = BASE + '/annotation/' + TRANSCRIPTOME_NAME + '.gtf.gz'

GENOME_NAME = 'Homo_sapiens.GRCh38.dna.primary_assembly.rel80'
GENOME_FA = BASE + '/annotation/' + GENOME_NAME + '.fa.gz'
