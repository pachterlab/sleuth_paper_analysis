###
# configuration that should be handled manually and will differ between systems
###

# BASE = '/home/hjp/sleuth_paper_analysis'
# BASE = '/Users/hjp/analysis/sleuth_paper'
# TODO: make this a absolute path before submission
from os.path import expanduser
from os import getenv

HOME = expanduser('~')
BASE = HOME + '/sleuth_paper_analysis'

N_THREADS = 35

###
# software
###
BIN = BASE + '/software/bin'

KALLISTO = BIN + '/kallisto'

RSEM_PATH = BASE + '/software/rsem_simulate'
RSEM_SIMULATE = RSEM_PATH + '/rsem-simulate-reads'
RSEM_SIMULATE = BIN + '/rsem-simulate-reads'

# UPDATED_PATH = 'PATH=' + ':'.join([
#     RSEM_PATH,
#     getenv('PATH')
#     ])
UPDATED_PATH = 'PATH=' + BIN + ':$PATH'
HISAT = BIN + '/hisat2'

# import os

# os.environ['PATH'] = RSEM_PATH + ':' + os.environ['PATH']
# print(RSEM_PATH + ':' + os.environ['PATH'])
# os.environ['PATH'] = RSEM_PATH + ':' + os.environ['PATH']
# print(type(os.environ['PATH']))

###
# annotations
###
TRANSCRIPTOME_NAME = 'Homo_sapiens.GRCh38.cdna.all.rel80'
TRANSCRIPTOME_FA = BASE + '/annotation/' + TRANSCRIPTOME_NAME + '.fa'
TRANSCRIPTOME_GTF = BASE + '/annotation/' + TRANSCRIPTOME_NAME + '.gtf'

GENOME_NAME = 'Homo_sapiens.GRCh38.dna.primary_assembly.rel80'
GENOME_FA = BASE + '/annotation/' + GENOME_NAME + '.fa'

# STAR_DIRECTORY = BASE + '/index/star_' + TRANSCRIPTOME_NAME
# STAR_INDEX = STAR_DIRECTORY + '/Genome'
HISAT_INDEX = BASE + '/index/' + GENOME_NAME

RSEM_ANNOTATION_DIR = '/'.join([
    BASE,
    'annotation',
    TRANSCRIPTOME_NAME + '_rsem'])
RSEM_ANNOTATION = RSEM_ANNOTATION_DIR + '/ref'
RSEM_MODEL = BASE + '/geuvadis/results/rsem/HG00365_7/out.stat/out.model'

###
# functions
###
def source_r(base, fname):
    return 'Rscript --vanilla --default-packages=methods,stats,utils -e \'setwd("{0}")\' -e \'source("{1}")\''.format(base, fname)
