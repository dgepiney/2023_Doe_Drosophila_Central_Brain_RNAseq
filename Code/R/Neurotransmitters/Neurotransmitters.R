#Load the data
load("/Users/gonzalomoraleschaya/Dropbox (University of Oregon)/23 T2 snRNAseq/z methods/Fluent-Parse integrated/Seurat/pip.parse1.parse2.integrated.short-ALL-T2.RData")

#Feature Plots Fast-acting neurotransmitters
FeaturePlot(object = T2_all.integrated.neurons, features = "VGlut", cols = c('grey95', 'darkgreen'), min.cutoff = 0)
FeaturePlot(object = T2_all.integrated.neurons, features = "VAChT", cols = c('grey95', 'darkred'), min.cutoff = 0)
FeaturePlot(object = T2_all.integrated.neurons, features = "Gad1", cols = c('grey95', 'darkblue'), min.cutoff = 0)

ggsave('../Desktop/VGlut.jpg', width = 13, height = 10)
ggsave('../Desktop/VAChT.jpg', width = 13, height = 10)
ggsave('../Desktop/Gad1.jpg', width = 13, height = 10)


#Venn Diagram
vGlut_cells <- WhichCells(T2_all.integrated.neurons, expression = VGlut > 2)
vAChT_cells <- WhichCells(T2_all.integrated.neurons, expression = VAChT > 2)
Gad1_cells <- WhichCells(T2_all.integrated.neurons, expression = Gad1 > 2)
all_cells <- T2_all.integrated.neurons$barcode


venn_data <- list(
  "Glutamate" = vGlut_cells,
  "GABA" = Gad1_cells,
  "Acetylcholine" = vAChT_cells)

plot <- venn.diagram(
  x = venn_data,
  category.names = c("Glutamate", "GABA", "Acetylcholine"), filename = NULL, scaled =TRUE)


