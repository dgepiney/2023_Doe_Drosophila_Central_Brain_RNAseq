library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(preprocessCore)
library(tidyverse)
library(ComplexHeatmap)
library(viridis)
library(circlize)

t2 = readRDS('~/Desktop/T2RNAseq/Adult/T2.atlas.neurons.gz')
DefaultAssay(t2) = 'RNA'
t2.data = as.data.frame(GetAssayData(t2, assay = 'RNA', layer='scale.data'))
tf = read.delim('/Users/gonzalomoraleschaya/University of Oregon Dropbox/Gonzalo Morales Chaya/23 T2 snRNAseq/_new T2rnaseq paper/other stuff/Data Plots/NP overview/Panel F/TFs.txt')
np = read.delim('/Users/gonzalomoraleschaya/University of Oregon Dropbox/Gonzalo Morales Chaya/23 T2 snRNAseq/_new T2rnaseq paper/other stuff/Data Plots/NP overview/Panel F/neuropeptides.txt')
np = as.character(np[,"SYMBOL"])
np = np[!np %in% c("Pburs", "SP")]
np = FindAllMarkers(object = t2, assay = "RNA", test.use = "wilcox",
                      features = np, only.pos = TRUE, 	
                      return.thresh = 0.01, logfc.threshold = 0.1, 
                      min.pct = 0.25, verbose = TRUE)
tf.data=t2.data[rownames(t2.data)%in%tf$SYMBOL,]
np.data=t2.data[rownames(t2.data)%in%unique(np$gene),]
cor.np.tf=cor(t(tf.data), t(np.data))
cor.np.tf=na.omit(cor.np.tf)
cor.np.tf=as.data.frame(t(cor.np.tf))
top.TFs = c()
for (np in rownames(cor.np.tf)) {
  row.values = cor.np.tf[np, ]
  row.values = as.numeric(row.values)
  names(row.values) = colnames(cor.np.tf)
  sorted.TFs = sort(row.values, decreasing = TRUE)
  this.top.5.TFs = names(head(sorted.TFs, 5))
  top.TFs = unique(c(top.TFs, this.top.5.TFs))
}
cor.np.tf = cor.np.tf[c("NPF","AstA", "CCHa2", "sNPF", "Ms", "AstC", "Dh44","spab", "Dh31", "Tk", "Mip", "Nplp1", "Proc"), ]
#Panel F
Heatmap(t(cor.np.tf[,top.TFs]), cluster_columns = FALSE, clustering_method_rows = "ward.D", width =nrow(cor.np.tf[,top.TFs])*unit(20, 'points'), column_dend_side = "bottom", column_names_side = "top")
#Supplemental Figure
np.list = read.delim('/Users/gonzalomoraleschaya/University of Oregon Dropbox/Gonzalo Morales Chaya/23 T2 snRNAseq/_new T2rnaseq paper/other stuff/Data Plots/NP overview/Panel F/neuropeptides.txt')
np.list = as.character(np.list[,"SYMBOL"])
np.list = np.list[!np.list %in% c("Pburs", "SP")]
np.list.data=t2.data[rownames(t2.data)%in%unique(np.list),]
cor.np.list.tf=cor(t(tf.data), t(np.list.data))
cor.np.list.tf=na.omit(cor.np.list.tf)
cor.np.list.tf=as.data.frame(t(cor.np.list.tf))
Heatmap(cor.np.list.tf, cluster_rows = T, cluster_columns = T, width=ncol(cor.np.list.tf)*unit(3,'points'), height=nrow(cor.np.list.tf)*unit(3, 'points'), 
        row_names_gp=gpar(fontsize=3), column_names_gp=gpar(fontsize=3))
#Supplemental Table
cor.np.list.tf <- cor.np.list.tf %>%
  rownames_to_column(var = "RowNames")
write_csv(cor.np.list.tf,"/Users/gonzalomoraleschaya/University of Oregon Dropbox/Gonzalo Morales Chaya/23 T2 snRNAseq/_new T2rnaseq paper/_current figs, tables, movies/sup figs & tables/Supplemental Table 9 Neuropeptides Neurotransmitter expression by cluster.csv")
