# deseq_analysis.R

# Función para ejecutar el análisis DESeq2
run_deseq_analysis <- function(readcounts, samplesheet) {
  deseq_dataset <- DESeqDataSetFromMatrix(countData = readcounts, colData = samplesheet, design = ~ condition)
  deseq_dataset <- DESeq(deseq_dataset)
  
  res <- results(deseq_dataset)
  resSig <- res[which(res$padj < 0.05), ]  # Filtrar genes con p-ajustada < 0.05
  
  return(list(deseq_dataset = deseq_dataset, res = res, resSig = resSig))
}
