top_5 <- function(file_list) {
    path <- "/Users/derek/Dropbox (University of Oregon)/lab_stuffs/r_data/"
    
    # Initialize an empty dataframe to store the cumulative result
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
  "cluster_001", "cluster_024", "cluster_002", "cluster_009", "cluster_029",
  "cluster_036", "cluster_064", "cluster_088", "cluster_113", "cluster_003",
  "cluster_112", "cluster_114", "cluster_133", "cluster_140", "cluster_008",
  "cluster_057", "cluster_010", "cluster_012", "cluster_051", "cluster_077",
  "cluster_016", "cluster_067", "cluster_021", "cluster_095", "cluster_123",
  "cluster_025", "cluster_093", "cluster_026", "cluster_033", "cluster_045",
  "cluster_048", "cluster_053", "cluster_090", "cluster_094", "cluster_117",
  "cluster_128", "cluster_143", "cluster_156", "cluster_043", "cluster_054",
  "cluster_109", "cluster_047", "cluster_170", "cluster_049", "cluster_050",
  "cluster_055", "cluster_056", "cluster_078", "cluster_084", "cluster_130",
  "cluster_131", "cluster_052", "cluster_083", "cluster_063", "cluster_082",
  "cluster_167", "cluster_168", "cluster_065", "cluster_108", "cluster_071",
  "cluster_116", "cluster_072", "cluster_086", "cluster_073", "cluster_103",
  "cluster_079", "cluster_102", "cluster_115", "cluster_132", "cluster_139",
  "cluster_081", "cluster_135", "cluster_091", "cluster_106", "cluster_104",
  "cluster_124", "cluster_107", "cluster_141", "cluster_118", "cluster_122",
  "cluster_120", "cluster_151", "cluster_160", "cluster_154", "cluster_173",
  "cluster_161", "cluster_171"
)
result <- unique(top_5(files_list))
print(result)

cluster_order <- c(
	1, 24, 2, 9, 29, 36, 64, 88, 113, 3, 6, 30, 137, 155, 4, 5, 11, 19, 37, 61,
	112, 114, 133, 140, 8, 57, 10, 12, 51, 77, 16, 67, 21, 95, 123, 25, 93, 26,
	33, 45, 48, 53, 90, 94, 117, 128, 143, 156, 43, 54, 109, 47, 170, 49, 50,
	55, 56, 78, 84, 130, 131, 52, 83, 63, 82, 167, 168, 65, 108, 71, 116, 72,
	86, 73, 103, 79, 102, 115, 132, 139, 81, 135, 91, 106, 104, 124, 107, 141,
	118, 122, 120, 151, 160, 154, 173, 161, 171
)

my_data <- subset(x = my_data, idents = cluster_order)

Idents(my_data) <- factor(x = my_data@active.ident, levels = cluster_order)

DotPlot(my_data, 
    features = result$X, 
    cols = c("white", "black"), 
    assay = "RNA"
) + coord_flip() + theme(
    axis.text.x=element_text(angle=45, hjust = 1, size=18), 
    axis.text.y=element_text(size=22)
)

sv_plt("hthtfs_top_5_marker_dotplot.pdf",
  25,
  55
)