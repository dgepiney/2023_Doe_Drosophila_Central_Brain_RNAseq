library(tidyverse)
library(Seurat)
library(patchwork)
library(harmony)
library(Matrix)
library(cowplot)
library(scCustomize)

#Generate T1-T2 atlas
#Import sequenced unsorted sample (=sample1)
path_unsorted.1='~/T2snRNAseq/unsorted_DGE_1/'
unsorted.mat.1 = readMM(paste0(path_unsorted.1,'DGE.mtx'))
unsorted.cell_meta.1 = read.delim(paste0(path_unsorted.1,'cell_metadata.csv'), stringsAsFactors = F, sep=',')
unsorted.genes.1 = read.delim(paste0(path_unsorted.1,'all_genes.csv'), stringsAsFactors = F, sep=',')
unsorted.cell_meta.1$bc_wells = make.unique(unsorted.cell_meta.1$bc_wells, sep='_dup')
rownames(unsorted.cell_meta.1) = unsorted.cell_meta.1$bc_wells
unsorted.genes.1$gene_name = make.unique(unsorted.genes.1$gene_name, sep=',')
colnames(unsorted.mat.1) = unsorted.genes.1$gene_name
rownames(unsorted.mat.1) = rownames(unsorted.cell_meta.1)
unsorted.mat_t.1 = t(unsorted.mat.1)
unsorted.mat_t.1 = unsorted.mat_t.1[(rownames(unsorted.mat_t.1) != ''),]
unsorted.mat_t.1 = as(unsorted.mat_t.1, 'dgCMatrix')
unsorted.1 = CreateSeuratObject(unsorted.mat_t.1, min.features=200, min.cells=3, meta.data=unsorted.cell_meta.1)
unsorted.1[['percent.mt']]=PercentageFeatureSet(unsorted.1, patter='^mt:')
unsorted.1 <- subset(unsorted.1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

path_unsorted.2='~/T2snRNAseq/unsorted_DGE_2/'
unsorted.mat.2 = readMM(paste0(path_unsorted.2,'DGE.mtx'))
unsorted.cell_meta.2 = read.delim(paste0(path_unsorted.2,'cell_metadata.csv'), stringsAsFactors = F, sep=',')
unsorted.genes.2 = read.delim(paste0(path_unsorted.2,'all_genes.csv'), stringsAsFactors = F, sep=',')
unsorted.cell_meta.2$bc_wells = make.unique(unsorted.cell_meta.2$bc_wells, sep='_dup')
rownames(unsorted.cell_meta.2) = unsorted.cell_meta.2$bc_wells
unsorted.genes.2$gene_name = make.unique(unsorted.genes.2$gene_name, sep=',')
colnames(unsorted.mat.2) = unsorted.genes.2$gene_name
rownames(unsorted.mat.2) = rownames(unsorted.cell_meta.2)
unsorted.mat_t.2 = t(unsorted.mat.2)
unsorted.mat_t.2 = unsorted.mat_t.2[(rownames(unsorted.mat_t.2) != ''),]
unsorted.mat_t.2 = as(unsorted.mat_t.2, 'dgCMatrix')
unsorted.2 = CreateSeuratObject(unsorted.mat_t.2, min.features=200, min.cells=3, meta.data=unsorted.cell_meta.2)
unsorted.2[['percent.mt']]=PercentageFeatureSet(unsorted.2, patter='^mt:')
unsorted.2 <- subset(unsorted.2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Integrate unsorted samples
unsorted.list <- list(unsorted.1, unsorted.2)
names(unsorted.list)<-c('unsorted.1','unsorted.2')
unsorted.list = lapply (X = unsorted.list , FUN = function(x) {
  x = NormalizeData(x)
  x = FindVariableFeatures(x, selection.method='vst', nfeatures=2000)
})
integrate.features=SelectIntegrationFeatures(unsorted.list)
unsorted.anchors=FindIntegrationAnchors(unsorted.list, anchor.features=integrate.features)
unsorted.integrated=IntegrateData(anchorset=unsorted.anchors)
DefaultAssay(unsorted.integrated)<-'integrated'
all.genes.unsorted=rownames(unsorted.integrated)
unsorted.integrated=ScaleData(unsorted.integrated, feature=all.genes.unsorted)
unsorted.integrated=RunPCA(unsorted.integrated, npcs=50)
unsorted.integrated=RunUMAP(unsorted.integrated, dims=1:50)
unsorted.integrated=FindNeighbors(unsorted.integrated, dims=1:50)
unsorted.integrated=FindClusters(unsorted.integrated, resolution=12)
DimPlot(unsorted.integrated, label=T, repel=T)+NoLegend()
Idents(unsorted.integrated)<-'T1'
Idents(unsorted.integrated, WhichCells(unsorted.integrated, expression = (flpD5>0 | redstinger>0 | unc84sfGFP>0)))<-'T2'
unsorted.integrated$Lineage <- Idents(unsorted.integrated)
Idents(unsorted.integrated)<-unsorted.integrated$seurat_clusters
DimPlot(unsorted.integrated, group.by = 'Lineage', cols=c('darkgreen','cyan'), raster=F)
unsorted.integrated.ids=read.csv('~/T2snRNAseq/Unsorted_identity_by_clusters.csv', header=T)
unsorted.integrated.ids.list=unsorted.integrated.ids$Identity
names(unsorted.integrated.ids.list)<-levels(unsorted.integrated)
unsorted.integrated <- RenameIdents(unsorted.integrated, unsorted.integrated.ids.list)
unsorted.integrated$Cell.Identity<-Idents(unsorted.integrated)
Idents(unsorted.integrated)<-'seurat_clusters'
DimPlot(unsorted.integrated, label=T, repel=T, group.by = 'Cell.Identity')+NoLegend() #Figure 1B
DimPlot(unsorted.integrated, group.by = 'Lineage', cols=c('darkgreen','cyan')) #Figure 1C
unsorted.integrated.markers=FindAllMarkers(unsorted.integrated)
T2.clusters=c(43, 71, 78, 82, 96, 99, 126, 142, 145, 146, 148)
T2.enriched.clusters=subset(unsorted.integrated, idents=T2.clusters)
DimPlot(T2.enriched.clusters, label=T)+NoLegend() #Figure 1C'
T2.enriched.clusters.markers=unsorted.integrated.markers[unsorted.integrated.markers$cluster %in% T2.clusters,]
T2.enrcihed.clusters.markers.top5=Extract_Top_Markers(T2.enriched.clusters.markers, num_genes = 5)
DotPlot(T2.enriched.clusters, features=unique(T2.enrcihed.clusters.markers.top5), cols=c('white','black'), col.min=0)+RotatedAxis() #Figure 1D

#generated T1 neuron only atlas form unsorted samples 
unsorted.T1=subset(unsorted.integrated, rna_flpD5==0&rna_redstinger==0&rna_unc84sfGFP==0)
DefaultAssay(unsorted.T1)<-'RNA'
unsorted.T1=NormalizeData(unsorted.T1)
all.genes.T1=rownames(unsorted.T1)
unsorted.T1=ScaleData(unsorted.T1, features=all.genes.T1)
unsorted.T1=RunPCA(unsorted.T1, features=VariableFeatures(unsorted.T1), assay='RNA', reduction.name = 'pca', npcs = 50)
unsorted.T1=RunUMAP(unsorted.T1, dims=1:50)
unsorted.T1=FindNeighbors(unsorted.T1, dims=1:50)
unsorted.T1=FindClusters(unsorted.T1, resolution=12)
Idents(unsorted.T1)<-unsorted.T1$Cell.Identity
unsorted.T1.neurons=subset(unsorted.T1, idents=c('Glia','Hemocytes','OL'), invert=T)
unsorted.T1.neurons=RunUMAP(unsorted.T1.neurons, dims=1:50)
DimPlot(unsorted.T1.neurons, group.by = 'Cell.Identity', label=T, repel=T)+NoLegend() #Figure 2A

#Heatmap of T1 neuron markers (Figure 2B)
library(ComplexHeatmap)
library(viridis)
Idents(unsorted.T1.neurons)<-unsorted.T1.neurons$seurat_clusters
unsorted.T1.neurons.markers=FindAllMarkers(unsorted.T1.neurons)
unsorted.T1.neurons.markers.top10=Extract_Top_Markers(unsorted.T1.neurons, num_genes = 10)
unsorted.T1.neurons.markers.top10.exp=DotPlot(unsorted.T1.neurons, features=unique(unsorted.T1.neurons.markers.top10), group.by = 'Cell.Identity')$data
unsorted.T1.neurons.markers.top10.exp.mat = unsorted.T1.neurons.markers.top10.exp %>% select(-pct.exp, -avg.exp) %>% pivot_wider(names_from = id, values_from = avg.exp.scaled, values_fill = 0) %>% as.data.frame()
row.names(unsorted.T1.neurons.markers.top10.exp.mat)=unsorted.T1.neurons.markers.top10.exp.mat$features.plot
unsorted.T1.neurons.markers.top10.exp.mat = unsorted.T1.neurons.markers.top10.exp.mat %>% select(-features.plot)
Heatmap(unsorted.T1.neurons.markers.top10.exp.mat, cluster_rows = F, cluster_columns = F, row_names_gp = gpar(fontsize=2), row_names_side = 'left', column_names_gp = gpar(fontsize=5), col=viridis(40), heatmap_legend_param=list(title=c('Scaled average expression'), at=c(-3, 0, 3), title_position = 'leftcenter-rot'))

#Generate T2 neuron atlas (Figure 2C)
#Import sequenced results from FACS
#RFP (=sample2)
path_rfp.1='~/T2snRNAseq/rfp_DGE_1/'
rfp.mat.1 = readMM(paste0(path_rfp.1,'DGE.mtx'))
rfp.cell_meta.1 = read.delim(paste0(path_rfp.1,'cell_metadata.csv'), stringsAsFactors = F, sep=',')
rfp.genes.1 = read.delim(paste0(path_rfp.1,'all_genes.csv'), stringsAsFactors = F, sep=',')
rfp.cell_meta.1$bc_wells = make.unique(rfp.cell_meta.1$bc_wells, sep='_dup')
rownames(rfp.cell_meta.1)=rfp.cell_meta.1$bc_wells
rfp.genes.1$gene_name = make.unique(rfp.genes.1$gene_name, sep=',')
colnames(rfp.mat.1) = rfp.genes.1$gene_name
rownames(rfp.mat.1) = rownames(rfp.cell_meta.1)
rfp.mat_t.1 = t(rfp.mat.1)
rfp.mat_t.1 = rfp.mat_t.1[(rownames(rfp.mat_t.1) != ''),]
rfp.mat_t.1 = as(rfp.mat_t.1, 'dgCMatrix')
rfp.1 = CreateSeuratObject(rfp.mat_t.1, min.features=200, min.cells=3, meta.data=rfp.cell_meta.1)
rfp.1[['percent.mt']]=PercentageFeatureSet(rfp.1, patter='^mt:')
rfp.1 <- subset(rfp.1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#GFP (=sample3)
path_gfp.1='~/T2snRNAseq/gfp_DGE_1/'
gfp.mat.1 = readMM(paste0(path_gfp.1,'DGE.mtx'))
gfp.cell_meta.1 = read.delim(paste0(path_gfp.1,'cell_metadata.csv'), stringsAsFactors = F, sep=',')
gfp.genes.1 = read.delim(paste0(path_gfp.1,'all_genes.csv'), stringsAsFactors = F, sep=',')
gfp.cell_meta.1$bc_wells = make.unique(gfp.cell_meta.1$bc_wells, sep='_dup')
rownames(gfp.cell_meta.1) = gfp.cell_meta.1$bc_wells
gfp.genes.1$gene_name = make.unique(gfp.genes.1$gene_name, sep=',')
colnames(gfp.mat.1) = gfp.genes.1$gene_name
rownames(gfp.mat.1) = rownames(gfp.cell_meta.1)
gfp.mat_t.1 = t(gfp.mat.1)
gfp.mat_t.1 = gfp.mat_t.1[(rownames(gfp.mat_t.1) != ''),]
gfp.mat_t.1 = as(gfp.mat_t.1, 'dgCMatrix')
gfp.1 = CreateSeuratObject(gfp.mat_t.1, min.features=200, min.cells=3, meta.data=gfp.cell_meta.1)
gfp.1[['percent.mt']] = PercentageFeatureSet(gfp.1, patter='^mt:')
gfp.1 <- subset(gfp.1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Load pipseq data
pipseq.data = Read10X(data.dir='~/T2snRNAseq/pip_raw_matrix_gz/', gene.column=2)
pipseq = CreateSeuratObject(pipseq.data, project='pipseq', min.cells=3, min.features=200)
pipseq[['percent.mt']] = PercentageFeatureSet(pipseq, patter='^mt:')
pipseq <- subset(pipseq, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

T2.unsorted.1=subset(unsorted.1, rna_unc84sfGFP>0|rna_redstinger>0|rna_flpD5>0, slot='counts')
T2.gfp.1=subset(gfp.1, unc84sfGFP>0|flpD5>0, slot='counts')
T2.rfp.1=subset(rfp.1, redstinger>0|flpD5>0, slot='counts')
T2.pip=subset(pipseq, unc84sfGFP>0|flpD5>0, slot='counts')
T2.unsorted.1$dataset<-'parse-unsorted'
T2.unsorted.1$methods<-'Parse'
T2.gfp.1$dataset<-'parse-unc84sfGFP+'
T2.gfp.1$methods<-'Parse'
T2.rfp.1$dataset<-'parse-redstinger+'
T2.rfp.1$methods<-'Parse'
T2.pip$dataset<-'pipseq'
T2.pip$methods<-'Fluent'

T2.list <- list(T2.unsorted.1, T2.gfp.1, T2.rfp.1, T2.pip)
names(T2.list)<-c('T2.unsorted.1', 'T2.gfp.1', 'T2.rfp.1', 'pip')
T2.list = lapply(X = T2.list, FUN = function (x) {
  DefaultAssay(x) <- 'RNA'
  x = NormalizeData(x)
  x = FindVariableFeatures(x, selection.method='vst', nfeatures=2000)
})
integrate.features.T2 = SelectIntegrationFeatures(T2.list)
T2.anchors = FindIntegrationAnchors(T2.list, anchor.features=integrate.features.T2)
T2.integrated = IntegrateData(anchorset = T2.anchors)

#RFP (=sample2)
path_rfp.2='~/T2snRNAseq/rfp_DGE_2/'
rfp.mat.2 = readMM(paste0(path_rfp.2,'DGE.mtx'))
rfp.cell_meta.2 = read.delim(paste0(path_rfp.2,'cell_metadata.csv'), stringsAsFactors = F, sep=',')
rfp.genes.2 = read.delim(paste0(path_rfp.2,'all_genes.csv'), stringsAsFactors = F, sep=',')
rfp.cell_meta.2$bc_wells = make.unique(rfp.cell_meta.2$bc_wells, sep='_dup')
rownames(rfp.cell_meta.2)=rfp.cell_meta.2$bc_wells
rfp.genes.2$gene_name = make.unique(rfp.genes.2$gene_name, sep=',')
colnames(rfp.mat.2) = rfp.genes.2$gene_name
rownames(rfp.mat.2) = rownames(rfp.cell_meta.2)
rfp.mat_t.2 = t(rfp.mat.2)
rfp.mat_t.2 = rfp.mat_t.2[(rownames(rfp.mat_t.2) != ''),]
rfp.mat_t.2 = as(rfp.mat_t.2, 'dgCMatrix')
rfp.2 = CreateSeuratObject(rfp.mat_t.2, min.features=200, min.cells=3, meta.data=rfp.cell_meta.2)
rfp.2[['percent.mt']]=PercentageFeatureSet(rfp.2, patter='^mt:')
rfp.2 <- subset(rfp.2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#GFP (=sample3)
path_gfp.2='~/T2snRNAseq/gfp_DGE_2/'
gfp.mat.2 = readMM(paste0(path_gfp.2,'DGE.mtx'))
gfp.cell_meta.2 = read.delim(paste0(path_gfp.2,'cell_metadata.csv'), stringsAsFactors = F, sep=',')
gfp.genes.2 = read.delim(paste0(path_gfp.2,'all_genes.csv'), stringsAsFactors = F, sep=',')
gfp.cell_meta.2$bc_wells = make.unique(gfp.cell_meta.2$bc_wells, sep='_dup')
rownames(gfp.cell_meta.2) = gfp.cell_meta.2$bc_wells
gfp.genes.2$gene_name = make.unique(gfp.genes.2$gene_name, sep=',')
colnames(gfp.mat.2) = gfp.genes.2$gene_name
rownames(gfp.mat.2) = rownames(gfp.cell_meta.2)
gfp.mat_t.2 = t(gfp.mat.2)
gfp.mat_t.2 = gfp.mat_t.2[(rownames(gfp.mat_t.2) != ''),]
gfp.mat_t.2 = as(gfp.mat_t.2, 'dgCMatrix')
gfp.2 = CreateSeuratObject(gfp.mat_t.2, min.features=200, min.cells=3, meta.data=gfp.cell_meta.2)
gfp.2[['percent.mt']] = PercentageFeatureSet(gfp.2, patter='^mt:')
gfp.2 <- subset(gfp.2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

T2.unsorted.2=subset(unsorted.2, unc84sfGFP>0|redstinger>0|flpD5>0, slot='counts')
T2.gfp.2=subset(gfp.2, unc84sfGFP>0|flpD5>0, slot='counts')
T2.rfp.2=subset(rfp.2, redstinger>0|flpD5>0, slot='counts')
T2.unsorted.2$dataset<-'parse-unsorted'
T2.unsorted.2$methods<-'Parse_2'
T2.gfp.2$dataset<-'parse-unc84sfGFP+'
T2.gfp.2$methods<-'Parse_2'
T2.rfp.2$dataset<-'parse-redstinger+'
T2.rfp.2$methods<-'Parse_2'

T2.list <- list(T2.unsorted.2, T2.gfp.2, T2.rfp.2)
names(T2.list.2)<-c('T2.unsorted.2', 'T2.gfp.2', 'T2.rfp.2')
T2.list.2 = lapply(X = T2.list.2, FUN = function (x) {
  DefaultAssay(x) <- 'RNA'
  x = NormalizeData(x)
  x = FindVariableFeatures(x, selection.method='vst', nfeatures=2000)
})
integrate.features.T2.2 = SelectIntegrationFeatures(T2.list.2)
T2.anchors.2 = FindIntegrationAnchors(T2.list.2, anchor.features=integrate.features.T2)
T2.integrated.2 = IntegrateData(anchorset = T2.anchors.2)
DefaultAssay(T2.integrated.2) <- 'integrated'
all.genes_T2.2 = rownames(T2.integrated.2)
T2.integrated.2 = ScaleData(T2.integrated.2, features=all.genes_T2.2)
T2.integrated.2 = RunPCA(T2.integrated.2, npcs=50)
T2.integrated.2 = RunUMAP(T2.integrated.2, dims=1:50)
T2.integrated.2 = FindNeighbors(T2.integrated.2, dims=1:50)
T2.integrated.2 = FindClusters(T2.integrated.2, resolution=12)

T2_all.list <- list(T2.integrated, T2.integrated.2)
names(T2_all.list)<-c('T2.integrated', 'T2.integrated.2')
T2_all.list = lapply(X = T2_all.list, FUN = function (x) {
  DefaultAssay(x) <- 'RNA'
  x = NormalizeData(x)
  x = FindVariableFeatures(x, selection.method='vst', nfeatures=2000)
})
integrate.features.T2_all = SelectIntegrationFeatures(T2_all.list)
T2_all.anchors = FindIntegrationAnchors(T2_all.list, anchor.features=integrate.features.T2_all)
T2_all.integrated = IntegrateData(anchorset = T2_all.anchors)
DefaultAssay(T2_all.integrated)<-'integrated'
all.genes_T2_all=rownames(T2_all.integrated)
T2_all.integrated=ScaleData(T2_all.integrated, features=all.genes_T2_all)
T2_all.integrated=RunPCA(T2_all.integrated, npcs=50)
T2_all.integrated=RunUMAP(T2_all.integrated, dims=1:50)
T2_all.integrated=FindNeighbors(T2_all.integrated, dims=1:50)
T2_all.integrated=FindClusters(T2_all.integrated, resolution=12)

T2_all.integrated.neurons=subset(T2_all.integrated, idents=c(38,142,20,31,7,46,70,147,13,23,34,17,80,18,15,14,174,146,87,92), invert=T)
T2.neuron.ids=read.csv('~/T2snRNAseq/T2_neuron_identity_by_cluster.csv', header=T)
T2.neuron.ids.list=T2.neuron.ids$Identity
names(T2.neuron.ids.list)<-levels(T2_all.integrated.neurons)
T2_all.integrated.neurons <- RenameIdents(T2_all.integrated.neurons, T2.neuron.ids.list)
T2_all.integrated.neurons$Cell.Identity<-Idents(T2_all.integrated.neurons)
Idents(T2_all.integrated.neurons)<-'seurat_clusters'
DimPlot(T2_all.integrated.neurons, group.by = 'Cell.Identity', label=T, repel=T)+NoLegend() #Figure 2C

#Heatmap of T2 neuron markers (Figure 2D)
T2_all.integrated.neurons.markers=FindAllMarkers(T2_all.integrated.neurons)
T2_all.integrated.neurons.markers.top10=Extract_Top_Markers(T2_all.integrated.neurons, num_genes = 10)
T2_all.integrated.neurons.markers.top10.exp=DotPlot(T2_all.integrated.neurons, features=unique(T2_all.integrated.neurons.markers.top10), group.by = 'Cell.Identity')$data
T2_all.integrated.neurons.markers.top10.exp.mat = T2_all.integrated.neurons.markers.top10.exp %>% select(-pct.exp, -avg.exp) %>% pivot_wider(names_from = id, values_from = avg.exp.scaled, values_fill = 0) %>% as.data.frame()
row.names(T2_all.integrated.neurons.markers.top10.exp.mat)=T2_all.integrated.neurons.markers.top10.exp.mat$features.plot
T2_all.integrated.neurons.markers.top10.exp.mat = T2_all.integrated.neurons.markers.top10.exp.mat %>% select(-features.plot)
Heatmap(T2_all.integrated.neurons.markers.top10.exp.mat, cluster_rows = F, cluster_columns = F, row_names_gp = gpar(fontsize=2), row_names_side = 'left', column_names_gp = gpar(fontsize=5), col=viridis(40), heatmap_legend_param=list(title=c('Scaled average expression'), at=c(-3, 0, 3), title_position = 'leftcenter-rot'))
