###
# load the permutations and metadata we will benchmark on
###

transcript_gene_mapping <- get_mouse_gene_names()
transcript_gene_mapping <- dplyr::select(transcript_gene_mapping, target_id,
  ext_gene, ens_gene)

# put the isoform information for EBSeq in global variables
NG_LIST <- GetNg(transcript_gene_mapping$target_id,
  transcript_gene_mapping$ens_gene)

metadata <- read.csv('../metadata/experiment.csv', stringsAsFactors = FALSE)
training_sets <- readRDS('../metadata/permutations.rds')

metadata <- dplyr::select(metadata, sample = run_accession, strain,
  library_name)

metadata <- dplyr::mutate(metadata,
  path = file.path('..', 'results', 'single', sample, 'kallisto', 'abundance.h5'))
metadata <- dplyr::mutate(metadata,
  condition = factor(strain, labels = c('A', 'B')))
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$sample

# add factors to every single configuration as well as paths
training_sets <- lapply(training_sets,
  function(x){
    x <- dplyr::rename(x, sample = run_accession)
    x <- dplyr::mutate(x,
      # Moved to generate_resampling.R
      # condition = factor(strain, levels = sort(unique(strain)),
      #   labels = c('A', 'B')),
      path = file.path('..', 'results', 'single', sample, 'kallisto',
        'abundance.h5')
      )
    x <- data.frame(x, stringsAsFactors = FALSE)
    rownames(x) <- x$sample

    x
  })

validation_sets <- lapply(training_sets,
  function(x) {
    validation_sample <- setdiff(metadata$sample, x$sample)
    validation_sample <- data.frame(sample = validation_sample,
      stringsAsFactors = FALSE)
    dplyr::inner_join(metadata, validation_sample, by = 'sample')
  })

validation_sets <- lapply(validation_sets,
  function(validation) {
    validation <- data.frame(validation)
    rownames(validation) <- validation$sample
    validation
  })
