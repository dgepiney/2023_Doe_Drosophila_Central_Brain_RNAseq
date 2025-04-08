library(ComplexHeatmap)
library(eulerr)

#Feature Plots Fast-acting neurotransmitters
#Panels A-A''
FeaturePlot(object = T2.atlas.neurons, features = "VGlut", 
            cols = c('grey95', 'darkgreen'), min.cutoff = 0, max.cutoff = 5)
ggsave('~/Desktop/VGlut.svg', width = 13, height = 10)
FeaturePlot(object = T2.atlas.neurons, features = "VAChT", 
            cols = c('grey95', 'darkred'), min.cutoff = 0, max.cutoff = 5)
ggsave('~/Desktop/VAChT.svg', width = 13, height = 10)
FeaturePlot(object = T2.atlas.neurons, features = "Gad1", 
            cols = c('grey95', 'darkblue'), min.cutoff = 0, max.cutoff = 5)
ggsave('~/Desktop/Gad1.svg', width = 13, height = 10)
FeaturePlot(object = T2.atlas.neurons, features = "ple", 
            cols = c('grey95', 'darkorange4'), min.cutoff = 0, max.cutoff = 5)
ggsave('~/Desktop/ple.svg', width = 13, height = 10)
FeaturePlot(object = T2.atlas.neurons, features = "SerT", 
            cols = c('grey95', 'black'), min.cutoff = 0, max.cutoff = 5)
ggsave('~/Desktop/SerT.svg', width = 13, height = 10)
FeaturePlot(object = T2.atlas.neurons, features = "Tbh", 
            cols = c('grey95', 'deeppink4'), min.cutoff = 0, max.cutoff = 5)
ggsave('~/Desktop/Tbh.svg', width = 13, height = 10)
FeaturePlot(object = T2.atlas.neurons, features = "Tdc2", 
            cols = c('grey95', 'turquoise4'), min.cutoff = 0, max.cutoff = 5)
ggsave('~/Desktop/Tdc2.svg', width = 13, height = 10)

#Neurotransmitter List
vGlut_cells = WhichCells(T2.atlas.neurons, expression = VGlut > 2)
vAChT_cells = WhichCells(T2.atlas.neurons, expression = VAChT > 2)
Gad1_cells = WhichCells(T2.atlas.neurons, expression = Gad1 > 2)
Ple_cells = WhichCells(T2.atlas.neurons, expression = ple > 2)
SerT_cells = WhichCells(T2.atlas.neurons, expression = SerT > 2)
Tbh_cells = WhichCells(T2.atlas.neurons, expression = Tbh > 2)
Tdc2_cells = WhichCells(T2.atlas.neurons, expression = Tdc2 > 2)
all_cells = rownames(T2.atlas.neurons@meta.data)

#Create a list of all NT expressing cells
NT.list = list(
  "Glutamatergic (VGlut)" = vGlut_cells,
  "GABAergic (Gad1)" = Gad1_cells,
  "Cholinergic (VAChT)" = vAChT_cells,
  "Dopaminergic (Ple)" = Ple_cells,
  "Serotonergic (SerT)" = SerT_cells,
  "Tyraminergic (Tbh)" = Tbh_cells,
  "Octopaminergic (Tdc2)" = Tdc2_cells
  )

#Create UpSet plot. Panel H
combination.mtx = make_comb_mat(NT.list)
combination.mtx = combination.mtx[comb_size(combination.mtx) >= 28]
ss = set_size(combination.mtx)
cs = comb_size(combination.mtx)
set_colors = c("darkred", "darkgreen", "darkblue", "deeppink4", "darkorange4", "black", "turquoise4",rep("grey",20))
comb_index = order(comb_degree(combination.mtx), -cs)
comb_col = set_colors[seq_along(comb_index)]
plot = UpSet(combination.mtx, right_annotation = upset_right_annotation(combination.mtx, add_numbers = TRUE, ylim = c(0, 30)), set_order = order(-ss), 
             comb_order = order(comb_degree(combination.mtx), -cs), comb_col = "black",
             top_annotation = upset_top_annotation(combination.mtx, add_numbers = TRUE, ylim = c(0, 20)))

env = environment(plot@top_annotation@anno_list[["intersection_size"]]@fun@fun)
env[["value"]] = round((env[["value"]] / total_cells) * 100, 1)
env2 <- plot@right_annotation@anno_list[["set_size"]]@fun@var_env
env2[["value"]] = round((env2[["value"]] / total_cells) * 100, 1)


draw(plot)