library(tidyverse)
library(Seurat)
library(patchwork)
library(harmony)
library(Matrix)
library(cowplot)
library(scCustomize)

####Figure 1. Cell atlas of central brain with single nuclei RNAseq
#Figure 1B
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


#####Figure 2. Cell atlas of neruons produced by T2 NBs
#Figure 2A
T2.atlas.neurons=read_rds('/T2.atlas.neurons.gz')
DimPlot(T2.atlas.neurons) 
#Figure 2B
T2.atlas.neuronsmarkers=FindAllMarkers(unsorted.T1.neurons)
T2.atlas.neuronsmarkers.top10=Extract_Top_Markers(unsorted.T1.neurons, num_genes = 10)
T2.atlas.neuronsmarkers.top10.exp=DotPlot(unsorted.T1.neurons, features=unique(T2.atlas.neuronsmarkers.top10), group.by = 'Cell.Identity')$data
T2.atlas.neuronsmarkers.top10.exp.mat = T2.atlas.neuronsmarkers.top10.exp %>% select(-pct.exp, -avg.exp) %>% pivot_wider(names_from = id, values_from = avg.exp.scaled, values_fill = 0) %>% as.data.frame()
row.names(T2.atlas.neuronsmarkers.top10.exp.mat)=T2.atlas.neuronsmarkers.top10.exp.mat$features.plot
T2.atlas.neuronsmarkers.top10.exp.mat = T2.atlas.neuronsmarkers.top10.exp.mat %>% select(-features.plot)
Heatmap(T2.atlas.neuronsmarkers.top10.exp.mat, cluster_rows = F, cluster_columns = F, row_names_gp = gpar(fontsize=2), row_names_side = 'left', column_names_gp = gpar(fontsize=5), col=viridis(40), heatmap_legend_param=list(title=c('Scaled average expression'), at=c(-3, 0, 3), title_position = 'leftcenter-rot'))

####Supplemental Figure 2. Cell atlas of neruons produced by T1 NBs
#Figure Supp 2A - generated T1 neuron only atlas form unsorted samples 
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
DimPlot(unsorted.T1.neurons, group.by = 'Cell.Identity', label=T, repel=T)+NoLegend() #(Supplemental Figure 2A)

#Supp Figure 2B - Heatmap of T1 neuron markers (Supplemental Figure 2B)
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

##### Figure 8. Mapping T2-derived neurons to UMAP clusters
#Figure 8A
library(Matrix)
library(scater)
library(Seurat)
library(patchwork)
library(SeuratDisk)
library(SeuratData)
library(SCopeLoomR)
library(parallel)
library(pheatmap)
library(preprocessCore)
library(nnls)
library(pheatmap)
library(ComplexHeatmap)
library(viridis)

#Converting Seurat Object to Loom
DefaultAssay(T2.atlas.neurons)<-'RNA'
TypeII_loom <- as.loom(T2.atlas.neurons, filename = "~/Desktop/TypeII.loom", verbose = T)

# Access single cell data. This will create three lists: one for all cell IDs, 
# one for all genes and 
cell_ids <- get_cell_ids(TypeII_loom)
scGenes <- get_genes(TypeII_loom)

# Getting a data frame of cell barcodes x clusters
cell_clusters <- T2.atlas.neurons@meta.data[c(19)]
cell_clusters

#load entire sc gene x cell matrix
print("Getting scMat")
scMat.loom <- get_dgem(TypeII_loom) 

#Normalizing their value to CPM
print("Normalizing to CPM")
scMat.original <- sweep(scMat.loom, 2, colSums(scMat.loom), `/`) * 1E6
scMat = scMat.original

# Prep TAPIN-seq data
library(readxl)
bulkEmat.PEG=as.data.frame(read_xlsx('/Turner-Evans_BulkRNASeq_updated.xlsx', sheet=1))
row.names(bulkEmat.PEG)=bulkEmat.PEG$gene_symbol
bulkEmat.PEG$gene_symbol=NULL
bulkEmat.PEG$EPG.avg=(bulkEmat.PEG$EPG_1+bulkEmat.PEG$EPG_2+bulkEmat.PEG$EPG_3+bulkEmat.PEG$EPG_4)/4
bulkEmat.PEG$PEG.avg=bulkEmat.PEG$PEG_1
cor.PEG=cor(bulkEmat.PEG[,1:6])
Heatmap(cor.PEG)

fb=read.csv('/FB_bulk.csv')
row.names(fb)=fb$gene_symbol
fb$gene_symbol=NULL
fb.mat=as(t(as.matrix(fb)), 'sparseMatrix')
fb.seurat=CreateSeuratObject(counts=fb.mat)
fb.seurat=NormalizeData(fb.seurat)
fb.exp=as.data.frame(GetAssayData(fb.seurat, assay='RNA', layer='data'))
fb.exp$cell.types=rownames(fb.exp)
fb.exp=fb.exp %>% filter(!grepl('ER5', cell.types))
fb.exp=fb.exp %>% filter(!grepl('ExR', cell.types))
fb.exp$cell.types=NULL
bulkEmat.fb=as.data.frame(t(fb.exp))
bulkEmat.fb$hDeltaK.avg=rowMeans(bulkEmat.fb[1:5])
bulkEmat.fb$FB6A.avg=rowMeans(bulkEmat.fb[6:7])
bulkEmat.fb$FB7A.avg=rowMeans(bulkEmat.fb[10:12])
bulkEmat.fb$`FB2I-a-b.avg`=rowMeans(bulkEmat.fb[c(15, 17:20)])
bulkEmat.fb$`FB6A-2.avg`=rowMeans(bulkEmat.fb[25:26])
cor.fb=cor(bulkEmat.fb[,1:27])
Heatmap(cor.fb, cluster_rows = F, cluster_columns = F)

sharedGenes.bulkEmat=intersect(rownames(bulkEmat.fb), rownames(bulkEmat.PEG))
bulkEmat.fb=bulkEmat.fb[sharedGenes.bulkEmat,]
bulkEmat.PEG=bulkEmat.PEG[sharedGenes.bulkEmat,]
bulkEmat=cbind(bulkEmat.fb, bulkEmat.PEG)
bulkEmat.fb=bulkEmat.fb[sharedGenes.bulkEmat,c('hDeltaK-r1','hDeltaK-r2','hDeltaK-r3','hDeltaK-r4','hDeltaK-r5','FB6A-r1','FB6A-r2','FB7A-r1', 'FB7A-r2','FB7A-r3','FB2I-a-b-r2','FB2I-a-b-r4','FB2I-a-b-r5','FB2I-a-b-r6','FB2I-a-b-r7','FB6A-2-r3','FB6A-2-r4')]
bulkEmat.PEG=bulkEmat.PEG[sharedGenes.bulkEmat,c('EPG_1','EPG_2','EPG_3','EPG_4','PEG_1')]

# Saves the clusters names in a list and then counts the # of clusters
scClusters <- unique(as.character(cell_clusters[,1]))
numSCclusters <- length(scClusters)

#Creates an empty data frame of genes x clusters
scClusterMat <- matrix(nrow=nrow(scMat), ncol=numSCclusters, data = 0)
rownames(scClusterMat) <- rownames(scMat)
colnames(scClusterMat) <- scClusters

#Averages gene expression by cell cluster. Last step we can do before having bulk RNAseq darta, next step is to normalize both. 
scClusterSummary <- 'mean'
if (scClusterSummary == "mean") {
  for (scCluster in scClusters) {
    scClusterMat[,scCluster] <- rowMeans(scMat[,rownames(cell_clusters)[cell_clusters[,1] == scCluster]])
  }
}

#Saved the shared genes between all genes in the scdata with the ones in bulk RNAseq
sharedGenes <- intersect(scGenes, rownames(bulkEmat))

print(paste0("Genes in scMat=",length(scGenes),
             ", bulkMat=",nrow(bulkEmat),
             ". shared=",length(sharedGenes)))

# Get sc and bulk matrices in register
bulkEmat     <- bulkEmat[sharedGenes,]
scMat        <- scMat[sharedGenes,]
scClusterMat <- scClusterMat[sharedGenes,]

# Prep matrices & quantile normalized version
renlMats <- list()
renlMats$bulkEmat <- as.matrix(bulkEmat,
                               dimnames=list(rownames(bulkEmat),
                                             colnames(bulkEmat)))

renlMats$scClusterMat  <- as.matrix(scClusterMat,
                                    dimnames=list(rownames(scClusterMat),
                                                  colnames(scClusterMat)))

for (matType in names(renlMats)) {
  newMatType <- paste0(matType,".qnl")
  t.colnames <- colnames(renlMats[[matType]])
  t.rownames <- rownames(renlMats[[matType]])
  renlMats[[newMatType]] <- normalize.quantiles(renlMats[[matType]])
  colnames(renlMats[[newMatType]]) <- t.colnames
  rownames(renlMats[[newMatType]]) <- t.rownames
}

for (matType in names(renlMats)) {
  renlMats[[matType]] <- log1p(renlMats[[matType]])
}


# Find all required marker gene sets
nnls.options = list(c(
    direction         = 'sumOfSc', #ScAsSumOfBulk
    nMarkers          = 200,
    markerType        = "both", #or bulkEmat, scClusterMat
    intercept         = TRUE,
    nl.hmap           = "column", #row, column, none, both
    qnl.bulk          = FALSE,
    qnl.sc            = FALSE
))

tx.nMarkers <- unique(unlist(lapply(nnls.options, function(x) {return(x['nMarkers'])})))
markerLists <- list()
for (nMarker in tx.nMarkers) {
  print(paste0("Looking for top-",nMarker," markers"))
  markerLists[[paste0("top",nMarker)]] <- list()
  for (matType in names(renlMats)) {
    print(paste0("-> matrix ",matType))
    curMarkers <- c()
    for (i in 1:ncol(renlMats[[matType]])) {
      print(paste0("on column ",i))
      
      fc_v_med <- renlMats[[matType]][,i] -
        apply(renlMats[[matType]][,-1 * i],1,median)
      
      fc_v_max <- renlMats[[matType]][,i] -
        apply(renlMats[[matType]][,-1 * i],1,max)
      
      colMarkers <- union(
        head(rownames(renlMats[[matType]][
          order(fc_v_med, decreasing=TRUE),]),
          n=as.numeric(nMarker)),
        head(rownames(renlMats[[matType]][
          order(fc_v_max, decreasing=TRUE),]),
          n=as.numeric(nMarker)))
      
      curMarkers <- unique(c(curMarkers, colMarkers))
    }
    markerLists[[paste0("top",nMarker)]][[matType]] <- curMarkers
  }
}

#runNLLS function
runNNLS <- function(x, y, intercept, method) {
  coef.rows <- colnames(x)
  coef.cols <- colnames(y)
  if (intercept) { x <- cbind(1, x) }
  if (method == "nnls") {
    curFit <- list(fits = list(),
                   coef = matrix(nrow=length(coef.rows),
                                 ncol=length(coef.cols),
                                 data=0,
                                 dimnames=list(coef.rows,
                                               coef.cols)))
    for (i in 1:ncol(y)) {
      curFit$fits[[i]] <- nnls(x, y[,i])
      if (intercept) {
        curFit$coef[,i] <- coef(curFit$fits[[i]])[-1]
      } else {
        curFit$coef[,i] <- coef(curFit$fits[[i]])
      }
    }
  } else if (method == "fcnnls") {
    curFit <- list()
    curFit$fits <- .fcnnls(x,y)
    return(curFit)
    if (intercept) {
      curFit$coef <- curFit$fits$coef[-1,]
    } else {
      curFit$coef <- curFit$fits$coef
    }
  }
  return(curFit)
}

# Start NNLS regression analyses
nnlsFits <- list()
for (optionInd in 1:length(nnls.options)) {
  print(paste0("Running regression #",optionInd))
  curOptions <- c(nnls.options[[optionInd]])
  
  library(NMF)
  nmf.options(cores=4)
  
  # set what matrix to use: quantileNormalize or not
  curMatTypes <- list()
  curMatTypes$bulkEmat <- "bulkEmat"
  curMatTypes$scClusterMat <- "scClusterMat"
  
  if (curOptions["qnl.bulk"]) {
    curMatTypes$bulkEmat <- "bulkEmat.qnl"}
  
  if (curOptions["qnl.sc"]) {
    curMatTypes$scClusterMat <- "scClusterMat.qnl" }
  
  curMats <- list()
  curMats$bulkEmat <- renlMats[[curMatTypes$bulkEmat]]
  curMats$scClusterMat<- renlMats[[curMatTypes$scClusterMat]]
  
  curNmarker <- paste0("top",curOptions["nMarkers"])
}  

# define markers
markerGenes <- list(bulkEmat = c(), scClusterMat = c())
for (matType in names(curMatTypes)) {
  print(paste0("on matrix ",matType))
  markerGenes[[matType]] <- markerLists[[curNmarker]][[curMatTypes[[matType]]]]
}
markerGenes$both <- unique(unlist(markerGenes))



# run NNLS on each marker set
print(paste0("Running nnls using markers ", as.list(curOptions)$markerType))
t.bMat  <- curMats$bulkEmat#[markerGenes[[curOptions$markerType]],]
t.scMat <- curMats$scClusterMat#[markerGenes[[curOptions$markerType]],]
t.scMat=as.matrix(t.scMat) %>% replace(is.na(.),0)
t.bMat=as.matrix(t.bMat) %>% replace(is.na(.),0)



# Start NNLS regression analyses
nnls.options = list(c(
  direction         = 'sumOfSc', #ScAsSumOfBulk
  nMarkers          = 200,
  markerType        = "both", #or bulkEmat, scClusterMat
  intercept         = TRUE,
  nl.hmap           = "column", #row, column, none, both
  qnl.bulk          = FALSE,
  qnl.sc            = FALSE
))

nnlsFits <- list()
library(NMF)
nmf.options(cores=4)

# set what matrix to use: quantileNormalize or not
curMatTypes <- list()
curMatTypes$bulkEmat <- "bulkEmat"
curMatTypes$scClusterMat <- "scClusterMat"

curMats <- list()
curMats$bulkEmat <- renlMats[[curMatTypes$bulkEmat]]
curMats$scClusterMat<- renlMats[[curMatTypes$scClusterMat]]
View(curOptions)
curNmarker <- paste0("top200")

if (curOptions['direction'] == "sumOfSc") {
nnlsMethod='nnls'
curFit <- runNNLS(x            = as.matrix(t.scMat),
                  y            = as.matrix(t.bMat),
                  intercept    = curOptions['intercept'],
                  method       = nnlsMethod)

outMat <- as.data.frame(curFit$coef)
outMat=outMat[rowMeans(outMat)!=0,c('EPG.avg','PEG.avg','hDeltaK-r1','FB6A.avg','FB6A-2.avg','FB7A-r3')]
fb.PEG.row_order=c('EPG.avg','PEG.avg','hDeltaK-r1','FB6A.avg','FB6A-2.avg','FB7A-r3')
outMat=outMat[c('102','75','95','100','76','164','162','40'),]
row_labels=structure(c('E-PG','P-EG','hDeltaK','FB6A(1)','FB6A(2)', 'FB7A'), names=c('EPG.avg','PEG.avg','hDeltaK-r1','FB6A.avg','FB6A-2.avg','FB7A-r3'))
pdf('~/all.bulk.pdf', width=5, height=5)
Heatmap(t(outMat), cluster_rows = F, row_order = fb.PEG.row_order, col=viridis(3), row_labels=row_labels[colnames(outMat)])
dev.off()
}

#Figure 8B
splits=DotPlot(T2.atlas.neurons, features=c('Wnt10','Lmpt','Dscam3','Ms','Cyp49a1','mtgo','ds','Octalpha2R','Dr','dnc','erm','klg','DAT'))$data
splits.exp=splits %>% select(-pct.exp, -avg.exp) %>% pivot_wider(names_from = id, values_from = avg.exp.scaled, values_fill = 0) %>% as.data.frame()
rownames(splits.exp)=splits.exp$features.plot
splits.exp$features.plot = NULL
splits.exp=t(splits.exp)
splits.exp=splits.exp[c('149','162','126','168','82','139','64','197','42','189','157','170'),]
Heatmap(splits.exp, cluster_rows = F, cluster_columns = F, row_names_side = 'left', column_names_side = 'top')

#Figure 8C
DotPlot(T2.atlas.neurons, features=c('toy','run','rho','shakB','Ggamma30A','Tdc2','ct','Lim1','Rx'), idents=c(105, 38, 123, 46, 66, 130), col.min=0, cols=c('white','black'))

