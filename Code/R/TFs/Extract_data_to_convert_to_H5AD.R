library(Seurat)
library(tidyverse)
library(Matrix)


#write matrix data (gene expression counts)
counts_matrix <- GetAssayData(
  my_data,
  assay = 'RNA',
  slot = 'counts',
)


writeMM(
  counts_matrix,
  file = paste0(
    file = '/Users/derek/Dropbox (University of Oregon)/lab_stuffs/RNAseq/seurat_to_h5ad//matrix.mtx'
  )
)


#write dimensional reduction matrix (PCA)
write.csv(
  my_data@reductions$pca@cell.embeddings,
  file = '/Users/derek/Dropbox (University of Oregon)/lab_stuffs/RNAseq/seurat_to_h5ad//pca.csv',
  quote = FALSE,
  row.names = FALSE
)


#write gene names
write.table(
  data.frame('gene' = rownames(counts_matrix)),
  file = '/Users/derek/Dropbox (University of Oregon)/lab_stuffs/RNAseq/seurat_to_h5ad//gene_names.csv',
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


view(my_data@meta.data)


#save metadata table
my_data$barcode <- colnames(my_data)


my_data$UMAP_1 <- my_data@reductions$umap@cell.embeddings[,1]


my_data$UMAP_2 <- my_data@reductions$umap@cell.embeddings[,2]


write.csv(
  my_data@meta.data,
  file = '/Users/derek/Dropbox (University of Oregon)/lab_stuffs/RNAseq/seurat_to_h5ad//metadata.csv',
  quote = FALSE,
  row.names = FALSE
)
