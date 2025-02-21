# load_data.R

library(DESeq2)
library(dplyr)

# Función para procesar Read Counts
process_readcounts <- function(file) {
  read.csv(file, header = TRUE, row.names = 1)
}

# Función para procesar Sample Sheet
process_samplesheet <- function(file) {
  read.csv(file, header = TRUE, row.names = 1)
}

# Función para procesar DESeq Results
process_deseqres <- function(file) {
  read.csv(file, header = TRUE, row.names = 1)
}
