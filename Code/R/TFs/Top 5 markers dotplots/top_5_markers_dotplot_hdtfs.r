top_5 <- function(file_list) {
    path <- "/Users/derek/Dropbox (University of Oregon)/lab_stuffs/r_data/"
    
# read output of "FindMarkers" function and take top 5 markers for each cluster
# based on adjusted p value
    cumulative_result <- data.frame(X = numeric(0))
    
    for (clusters in file_list) {
        cluster <- read.csv(paste0(path, clusters, ".csv"))
        cluster_ordered <- cluster[order(cluster$p_val_adj), ]
        cluster_top_5 <- head(cluster_ordered$X, 5)
        
        # Append the top 5 values of the "X" column to the cumulative result
        cumulative_result <- rbind(cumulative_result, data.frame(X = cluster_top_5))
    }
    
    # Reset row names to avoid duplicates
    rownames(cumulative_result) <- NULL
    
    return(cumulative_result)
}


files_list <- c(
  "cluster_003", "cluster_004", "cluster_004", "cluster_030", "cluster_109",
  "cluster_115", "cluster_112", "cluster_027", "cluster_127", "cluster_140", 
  "cluster_128", "cluster_045", "cluster_130", "cluster_156", "cluster_149", 
  "cluster_090", "cluster_151", "cluster_152", "cluster_016", "cluster_024",
  "cluster_168", "cluster_076", "cluster_035", "cluster_097", "cluster_040",
  "cluster_060", "cluster_043", "cluster_057", "cluster_054", "cluster_091",
  "cluster_002", "cluster_029"
)
result <- unique(top_5(files_list))
print(result)


# cluster order based on clusters with same binarized gene expression
cluster_order <- c(
    109, 115, 112, 27, 127, 140, 128, 45, 130, 156, 149, 90, 151, 152, 16, 24,
    168, 76, 2, 29, 3, 30, 4, 5, 35, 97, 40, 60, 43, 57, 54, 91
)


# subset data to only clusters being plotted
my_data <- subset(x = my_data, idents = cluster_order)


# change cluster order for dotplot
Idents(my_data) <- factor(x = my_data@active.ident, levels = cluster_order)


DotPlot(my_data, 
    features = result$X, 
    cols = c("white", "black"), 
    assay = "RNA"
) + coord_flip() + theme(
    axis.text.x=element_text(angle=45, hjust = 1, size=18), 
    axis.text.y=element_text(size=22)
)

# save plot with custom ggplot function
sv_plt("hdtf_top_5_marker_dotplot.pdf", 25, 30)