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
  "cluster_000", "cluster_101", "cluster_001", "cluster_010", "cluster_012",
  "cluster_036", "cluster_051", "cluster_054", "cluster_067", "cluster_068",
  "cluster_086", "cluster_113", "cluster_002", "cluster_009", "cluster_021",
  "cluster_029", "cluster_041", "cluster_064", "cluster_077", "cluster_088",
  "cluster_123", "cluster_003", "cluster_059", "cluster_004", "cluster_047",
  "cluster_128", "cluster_005", "cluster_008", "cluster_011", "cluster_045",
  "cluster_057", "cluster_082", "cluster_089", "cluster_006", "cluster_030",
  "cluster_025", "cluster_037", "cluster_040", "cluster_075", "cluster_083",
  "cluster_032", "cluster_127", "cluster_131", "cluster_035", "cluster_043",
  "cluster_103", "cluster_055", "cluster_069", "cluster_063", "cluster_118",
  "cluster_065", "cluster_137", "cluster_081", "cluster_158", "cluster_100",
  "cluster_117"
)
result <- unique(top_5(files_list))
print(result)

# cluster order based on clusters with same binarized gene expression
cluster_order <- c(
    0, 101, 1, 10, 12, 16, 24, 28, 36, 51, 54, 67, 68, 86, 113, 2, 9, 21, 29,
    41, 64, 77, 88, 123, 3, 59, 4, 47, 128, 5, 8, 11, 45, 57, 82, 89, 6, 30, 25,
    37, 40, 75, 83, 32, 127, 131, 35, 43, 60, 48, 71, 85, 94, 111, 49, 96, 52,
    66, 103, 55, 69, 63, 118, 65, 137, 81, 158, 100, 117
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
sv_plt("bdtfs_top_5_marker_dotplot.pdf", 25, 55)