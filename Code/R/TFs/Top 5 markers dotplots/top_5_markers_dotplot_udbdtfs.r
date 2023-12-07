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
  "cluster_000", "cluster_002", "cluster_009", "cluster_101", "cluster_001",
  "cluster_016", "cluster_024", "cluster_062", "cluster_067", "cluster_113",
  "cluster_172", "cluster_003", "cluster_004", "cluster_005", "cluster_089",
  "cluster_006", "cluster_131", "cluster_143", "cluster_155", "cluster_008",
  "cluster_030", "cluster_075", "cluster_096", "cluster_097", "cluster_158",
  "cluster_011", "cluster_055", "cluster_022", "cluster_059", "cluster_110",
  "cluster_121", "cluster_127", "cluster_133", "cluster_151", "cluster_152",
  "cluster_160", "cluster_025", "cluster_093", "cluster_112", "cluster_115",
  "cluster_026", "cluster_027", "cluster_028", "cluster_047", "cluster_048",
  "cluster_056", "cluster_058", "cluster_074", "cluster_076", "cluster_078",
  "cluster_083", "cluster_091", "cluster_098", "cluster_107", "cluster_129",
  "cluster_134", "cluster_149", "cluster_157", "cluster_166", "cluster_029",
  "cluster_064", "cluster_077", "cluster_088", "cluster_032", "cluster_035",
  "cluster_043", "cluster_054", "cluster_072", "cluster_039", "cluster_040",
  "cluster_071", "cluster_116", "cluster_042", "cluster_081", "cluster_090",
  "cluster_099", "cluster_137", "cluster_084", "cluster_104", "cluster_111",
  "cluster_122", "cluster_124", "cluster_130", "cluster_095", "cluster_123",
  "cluster_106", "cluster_139", "cluster_165", "cluster_114", "cluster_159",
  "cluster_118", "cluster_132", "cluster_140"
)
result <- unique(top_5(files_list))
print(result)

# cluster order based on clusters with same binarized gene expression
cluster_order <- c(
    0, 2, 9, 101, 1, 16, 24, 62, 67, 113, 172, 3, 4, 5, 89, 6, 131, 143, 155, 8,
    30, 75, 96, 97, 158, 11, 55, 22, 59, 110, 121, 127, 133, 151, 152, 160, 25,
    93, 112, 115, 26, 27, 28, 47, 48, 56, 58, 74, 76, 78, 83, 91, 98, 107, 129,
    134, 149, 157, 166, 29, 64, 77, 88, 32, 35, 43, 54, 72, 39, 40, 71, 116, 42,
    171, 49, 68, 50, 153, 53, 167, 57, 144, 61, 154, 162, 63, 135, 148, 168, 65,
    81, 90, 99, 137, 84, 104, 111, 122, 124, 130, 95, 123, 106, 139, 165, 114,
    159, 118, 132, 140
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
sv_plt("udbdtfs_top_5_marker_dotplot.pdf", 25, 55)