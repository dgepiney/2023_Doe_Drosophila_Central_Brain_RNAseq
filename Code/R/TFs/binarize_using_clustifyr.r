library(clustifyr)


# Read in single-cell data
mat <- ReadMtx(
    "/Users/derek/Dropbox (University of Oregon)/lab_stuffs/RNAseq/seurat_to_h5ad/matrix.mtx",
    "/Users/derek/Dropbox (University of Oregon)/lab_stuffs/RNAseq/cells.csv",
    "/Users/derek/Dropbox (University of Oregon)/lab_stuffs/RNAseq/seurat_to_h5ad/gene_names.csv",
    cell.column = 1,
  feature.column = 1,
  cell.sep = "",
  feature.sep = "",
  skip.cell = 1,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)


# Binaize genes
binary <- binarize_expr(mat, n = 1000, cut = 0)


# Write binary output to csv file
write.csv(
  binary,
  "/Users/derek/Dropbox (University of Oregon)/lab_stuffs/RNAseq/hmgbtfs_binary_genes.csv", 
  row.names = TRUE
)

# Read csv file that was just created
binary <- read.csv(
  "/Users/derek/Dropbox (University of Oregon)/lab_stuffs/RNAseq/hmgbtfs_binary_genes.csv",
  header = TRUE
)


# Creates list of genes to keep 
keep <- c(hmgbtfs$SYMBOL)


# Goes through genes and removes any not in keep list
binary <- binary[binary[,1] %in% keep,]


#
write.csv(
  binary,
  "/Users/derek/Dropbox (University of Oregon)/lab_stuffs/RNAseq/hmgbtfs_binary_genes.csv", 
  row.names = TRUE
)
