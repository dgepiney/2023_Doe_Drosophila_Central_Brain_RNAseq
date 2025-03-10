library(tidyverse)
library(Seurat)
library(patchwork)
library(Matrix)
library(cowplot)

## T2 atlas

# Read Pip-seq (first round of sequencing)
PIPseq = Read10X(data.dir='/Volumes/DoeLab65TB/DOE_lab_sequencing_data/Sen-Lin/Fluent-PIPseq/Alignment/pip (gfp+unsorted)/raw_matrix_gz', gene.column=2)
PIPseq = CreateSeuratObject(PIPseq, project='PIPseq', min.cells=3, 
                            min.features=200)
PIPseq[['percent.mt']] = PercentageFeatureSet(PIPseq, patter='^mt:')
PIPseq = subset(PIPseq, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
PIPseq = subset(PIPseq, unc84sfGFP>0|flpD5>0, slot='counts')
PIPseq$sample = 'unknown'
PIPseq$method = 'PIPseq'

# Read Parse first round of sequencing
for (i in 1:3) {
  path = paste0('/Volumes/DoeLab65TB/DOE_lab_sequencing_data/Sen-Lin/Parse/T2_Lineages/output/2023_Jun/All_combined_output/sample', i,'/DGE_unfiltered')
  
  cell_meta = read.delim(paste0(path, '/cell_metadata.csv'), 
                         stringsAsFactors = F, sep = ',')
  genes = read.delim(paste0(path, '/all_genes.csv'), 
                     stringsAsFactors = F, sep = ',')
  mat = readMM(paste0(path, '/DGE.mtx'))
  
  cell_meta$bc_wells = make.unique(cell_meta$bc_wells, sep = '_dup')
  rownames(cell_meta) = cell_meta$bc_wells
  
  genes$gene_name = make.unique(genes$gene_name, sep = ',')
  colnames(mat) = genes$gene_name
  rownames(mat) = rownames(cell_meta)
  
  mat = t(mat)
  mat = mat[(rownames(mat) != ''),]
  mat = as(mat, 'dgCMatrix')
  
  seurat_obj = CreateSeuratObject(mat, min.features = 200, min.cells = 3, 
                                  meta.data = cell_meta)
  seurat_obj[['percent.mt']] = PercentageFeatureSet(seurat_obj, pattern = '^mt:')
  seurat_obj = subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  seurat_obj$method = paste0('parse.1.',i)
  seurat_obj$sample = 'unknown'
  seurat_obj = subset(seurat_obj, redstinger>0|unc84sfGFP>0|flpD5>0, slot='counts')
  assign(paste0('parse.1.',i), seurat_obj)
}

# Read Parse second round of sequencing (P4&S4) Sept 2024
path = '/Volumes/DoeLab65TB/DOE_lab_sequencing_data/Sen-Lin/Parse/T2_Lineages_2/output/GC3F-SL-9731/ParseAnalyses/P4+S4/20241006_X4_1t8.combined/all-sample/DGE_unfiltered'

cell_meta = read.delim(paste0(path, '/cell_metadata.csv'), stringsAsFactors = F, sep = ',')
genes = read.delim(paste0(path, '/all_genes.csv'), stringsAsFactors = F, sep = ',')
mat = readMM(paste0(path, '/count_matrix.mtx'))

cell_meta$bc_wells = make.unique(cell_meta$bc_wells, sep = '_dup')
rownames(cell_meta) = cell_meta$bc_wells

genes$gene_name = make.unique(genes$gene_name, sep = ',')
colnames(mat) = genes$gene_name
rownames(mat) = rownames(cell_meta)

mat = t(mat)
mat = mat[(rownames(mat) != ''),]
mat = as(mat, 'dgCMatrix')

seurat_obj = CreateSeuratObject(mat, min.features = 200, min.cells = 3, 
                                meta.data = cell_meta)
seurat_obj[['percent.mt']] = PercentageFeatureSet(seurat_obj, pattern = '^mt:')
seurat_obj = subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
seurat_obj$method = 'parse.2.0'
parse.2.0 = subset(seurat_obj, redstinger>0|flpD5>0, slot='counts')

# Read Parse S4 sublibrary 10
path = '/Volumes/DoeLab65TB/DOE_lab_sequencing_data/Sen-Lin/Parse/T2_Lineages_2/output/GC3F-SL-9731/ParseAnalyses/S4/20240924_S4_sublib_10/all-sample/DGE_unfiltered'

cell_meta = read.delim(paste0(path, '/cell_metadata.csv'), stringsAsFactors = F, sep = ',')
genes = read.delim(paste0(path, '/all_genes.csv'), stringsAsFactors = F, sep = ',')
mat = readMM(paste0(path, '/count_matrix.mtx'))

cell_meta$bc_wells = make.unique(cell_meta$bc_wells, sep = '_dup')
rownames(cell_meta) = cell_meta$bc_wells

genes$gene_name = make.unique(genes$gene_name, sep = ',')
colnames(mat) = genes$gene_name
rownames(mat) = rownames(cell_meta)

mat = t(mat)
mat = mat[(rownames(mat) != ''),]
mat = as(mat, 'dgCMatrix')

seurat_obj = CreateSeuratObject(mat, min.features = 200, min.cells = 3, 
                                meta.data = cell_meta)
seurat_obj[['percent.mt']] = PercentageFeatureSet(seurat_obj, pattern = '^mt:')
seurat_obj = subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
seurat_obj$method = 'parse.2.10'
parse.2.10 = subset(seurat_obj, redstinger>0|flpD5>0, slot='counts')

#Integration of different rounds
T2.atlas = merge(PIPseq, y = c(parse.1.1, parse.1.2, parse.1.3, parse.2.0, parse.2.10))
T2.atlas = NormalizeData(T2.atlas)
T2.atlas = FindVariableFeatures(T2.atlas)
T2.atlas = ScaleData(T2.atlas, features = rownames(T2.atlas))
T2.atlas = RunPCA(T2.atlas)
T2.atlas = IntegrateLayers(object = T2.atlas, method = RPCAIntegration, 
                           orig.reduction = 'pca', new.reduction = 'integrated.rpca', verbose = T)
T2.atlas[['RNA']] = JoinLayers(T2.atlas[['RNA']])
T2.atlas = FindNeighbors(T2.atlas, reduction = 'integrated.rpca', dims = 1:50)
T2.atlas = FindClusters(T2.atlas, resolution = 12)
T2.atlas = RunUMAP(T2.atlas, dims=1:50, reduction='integrated.rpca', reduction.name='umap.rpca')

#Visualization
DimPlot(T2.atlas, reduction='umap.rpca', label = T, repel = T)+NoLegend()
FeaturePlot(T2.atlas, features = 'repo')

#Subseting and reclustering only neurons
T2.atlas.neurons = subset(T2.atlas, idents = c(78, 9, 104, 121, 176, 16, 137, 10, 
                                               196, 35, 144, 99, 194, 28, 94, 
                                               184, 134, 116, 63, 44, 145, 22, 
                                               26, 1, 110, 44, 57, 171, 125, 8, 
                                               17, 23, 57, 37, 11, 159, 61, 178, 
                                               179), invert = T)

T2.atlas.neurons = RunUMAP(T2.atlas.neurons, dims=1:50, reduction='integrated.rpca', reduction.name='umap.rpca')
