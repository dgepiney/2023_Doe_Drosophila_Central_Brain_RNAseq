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
  "cluster_003",
  "cluster_004",
  "cluster_004",
  "cluster_030",
  "cluster_109",
  "cluster_115",
  "cluster_112",
  "cluster_027",
  "cluster_127",
  "cluster_140",
  "cluster_128",
  "cluster_045",
  "cluster_130",
  "cluster_156",
  "cluster_149",
  "cluster_090",
  "cluster_151",
  "cluster_152",
  "cluster_016",
  "cluster_024",
  "cluster_168",
  "cluster_076",
  "cluster_035",
  "cluster_097",
  "cluster_040",
  "cluster_060",
  "cluster_043",
  "cluster_057",
  "cluster_054",
  "cluster_091",
  "cluster_002",
  "cluster_029"
)
result <- unique(top_5(files_list))
print(result)

DotPlot(my_data, features = result, assay = "RNA") + coord_flip()