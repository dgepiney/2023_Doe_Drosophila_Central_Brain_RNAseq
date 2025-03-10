library(Seurat)
library(ggplot2)

# read RNAseq data
my_data <- readRDS('/Users/derek/University of Oregon Dropbox/Derek Epiney/lab_stuffs/r_data/t2_rnaseq.rds')

# list of non-unique clusters
clusters_to_keep <- c(24, 30, 34, 55, 62, 75, 95, 65, 76, 66, 130)

# read list of top 5 markers for non-uniq clusters
df <- read.csv("/Users/derek/Desktop/RNAseq/hdtf_nonunique_clusters_top_five_markers.csv")

# ensure no repeats of markers in list
markers <- unique(df$gene)

print(markers)

# subset data
my_data <- subset(my_data, idents = clusters_to_keep)

# rearrange clusters
Idents(my_data) <- factor(Idents(my_data), levels=clusters_to_keep)

# create dotplot
print(
    DotPlot(
        my_data, idents=c(24, 30, 34, 55, 62, 75, 95, 65, 76, 66, 130),
        features = markers, cols = c("white", "black"),
        assay = "RNA", col.min = 0) + coord_flip() + 
        theme(axis.text.x = element_text(angle = 90, size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18)
    )
)

# save dotplot
ggsave(filename = file.path("/Users/derek/Dropbox (University of Oregon)/lab_stuffs/plots/hdtf_top_five_markers_dotplot.pdf"), plot = last_plot(), 
        height = 3400, width = 2200, units = "px", dpi = 300, 
        limitsize = FALSE, )