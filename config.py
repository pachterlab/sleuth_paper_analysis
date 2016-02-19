BASE = '/home/hjp/sleuth_paper_analysis'

# software
KALLISTO = BASE + '/software/kallisto_linux/kallisto'

# functions
def source_r(base, fname):
    return 'Rscript --vanilla --default-packages=methods,stats,utils -e \'setwd("{0}")\' -e \'source("{1}")\''.format(base, fname)
