# visualization.R
library(ggplot2)

# Función para visualizar el tamaño de las bibliotecas
plot_library_sizes <- function(deseq_dataset) {
  library_sizes <- colSums(counts(deseq_dataset))
  library_sizes_df <- data.frame(Sample = names(library_sizes), Library_Size = library_sizes)
  
  ggplot(library_sizes_df, aes(x = Sample, y = Library_Size)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = "Tamaño de las bibliotecas", x = "Muestras", y = "Número de lecturas") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Función para visualizar la expresión de un gen
plot_gene_expression <- function(deseq_dataset, gene) {
  plotCounts(deseq_dataset, gene = gene, intgroup = "condition")
}
