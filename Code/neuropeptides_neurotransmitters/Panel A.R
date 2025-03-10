#General View of neuropeptides
library(scCustomize)
library(BiocManager)
library(ComplexHeatmap)
library(Seurat)
library(ggplot2)
T2.atlas.neurons = readRDS("~/Desktop/T2.atlas.neurons.gz")
neuropeptides = read.delim("/Users/gonzalomoraleschaya/University of Oregon Dropbox/Gonzalo Morales Chaya/23 T2 snRNAseq/_new T2rnaseq paper/other stuff/Data Plots/NP overview/Panel A/neuropeptides.txt")
neuropeptides = as.character(neuropeptides$SYMBOL)
neuropeptides = neuropeptides[!neuropeptides %in% c("Pburs", "SP")]
cluster.defining.NP = FindAllMarkers(object = T2.atlas.neurons, assay = "RNA", 
                                     test.use = "wilcox",features = neuropeptides, 
                                     only.pos = T, return.thresh = 0.01, 
                                     logfc.threshold = 0.1, min.pct = 0.25, 
                                     verbose = T)
T2.neuropeptides = unique(cluster.defining.NP$gene)
T2.atlas.neurons.peptidergic = subset(T2.atlas.neurons, idents = unique(cluster.defining.NP$cluster))
Clustered_DotPlot(seurat_object = T2.atlas.neurons.peptidergic, features = T2.neuropeptides, colors_use_exp = c("grey90","grey40","black"))
neuropeptides.ordered = c("Proc", "Nplp1", "Mip", "Tk", "Dh31", "spab", "Dh44", "AstC", "Ms", "sNPF", "CCHa2", "AstA", "NPF")
neurotransmitters = c("ple", "SerT", "Tbh", "Tdc2", "VGlut", "VAChT", "Gad1")
all_features = c(neurotransmitters, neuropeptides.ordered)
clusters.ordered = c("100","128","86","197","195","124","160","157","140","52","48","119","133","41","82","64","120","174","192","155","169",
                      "76","163","65","49","50","150","27","173","93","170","162","113","12","56","20","13","115","111","98","127","25","84",
                      "188","156","152","40","168","54","147","106","108","166","149","186","165","58","77","191","47","143","187","74","167",
                      "88","71","75","62","46","66","122","130","190","97","154","158","126","183")
Idents(T2.atlas.neurons.peptidergic) = factor(Idents(T2.atlas.neurons.peptidergic), levels = clusters.ordered)
#Panel A
DotPlot(T2.atlas.neurons.peptidergic, features = all_features, cols = c("white","black")) + coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  theme(axis.text.y = element_text(size = 10)) + NoLegend()
#Dot Plot of all genes by all clusters
DotPlot(T2.atlas.neurons, features = c(neuropeptides, neurotransmitters), cols = c("white","black"))+ coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  theme(axis.text.y = element_text(size = 10)) + NoLegend()
#Table
Supp.Table = DotPlot(T2.atlas.neurons, features = c(neuropeptides, neurotransmitters))$data
write.csv(Supp.Table, "/Users/gonzalomoraleschaya/University of Oregon Dropbox/Gonzalo Morales Chaya/23 T2 snRNAseq/_new T2rnaseq paper/other stuff/Data Plots/NP overview/Supplemental Table.csv")
