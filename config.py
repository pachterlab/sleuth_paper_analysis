BASE = '/home/hjp/sleuth_paper_analysis'
# BASE = '/Users/hjp/analysis/sleuth_paper'

# software
KALLISTO = BASE + '/software/kallisto_linux/kallisto'
RSEM_SIMULATE = BASE + '/software/rsem_simulate/rsem-simulate-reads'

# functions
def source_r(base, fname):
    return 'Rscript --vanilla --default-packages=methods,stats,utils -e \'setwd("{0}")\' -e \'source("{1}")\''.format(base, fname)
