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
  "cluster_000", "cluster_001", "cluster_002", "cluster_009", "cluster_010",
  "cluster_011", "cluster_012", "cluster_016", "cluster_021", "cluster_022",
  "cluster_026", "cluster_028", "cluster_029", "cluster_032", "cluster_036",
  "cluster_042", "cluster_048", "cluster_049", "cluster_050", "cluster_051",
  "cluster_053", "cluster_054", "cluster_055", "cluster_056", "cluster_059",
  "cluster_063", "cluster_064", "cluster_065", "cluster_067", "cluster_069",
  "cluster_073", "cluster_076", "cluster_077", "cluster_078", "cluster_082",
  "cluster_084", "cluster_086", "cluster_094", "cluster_095", "cluster_101",
  "cluster_103", "cluster_104", "cluster_105", "cluster_107", "cluster_109",
  "cluster_110", "cluster_113", "cluster_114", "cluster_117", "cluster_118",
  "cluster_123", "cluster_124", "cluster_127", "cluster_128", "cluster_130",
  "cluster_131", "cluster_133", "cluster_135", "cluster_137", "cluster_138",
  "cluster_139", "cluster_143", "cluster_148", "cluster_149", "cluster_150",
  "cluster_152", "cluster_153", "cluster_157", "cluster_160", "cluster_161",
  "cluster_163", "cluster_165", "cluster_166", "cluster_168", "cluster_172",
  "cluster_003", "cluster_004", "cluster_005", "cluster_006", "cluster_008",
  "cluster_042", "cluster_048", "cluster_049", "cluster_050", "cluster_051",
  "cluster_019", "cluster_024", "cluster_027", "cluster_030", "cluster_035",
  "cluster_037", "cluster_045", "cluster_047", "cluster_052", "cluster_057",
  "cluster_061", "cluster_062", "cluster_071", "cluster_079", "cluster_085",
  "cluster_089", "cluster_090", "cluster_093", "cluster_100", "cluster_102",
  "cluster_108", "cluster_112", "cluster_119", "cluster_121", "cluster_129",
  "cluster_132", "cluster_136", "cluster_140", "cluster_151", "cluster_155",
  "cluster_156", "cluster_162", "cluster_164", "cluster_171", "cluster_025",
  "cluster_125", "cluster_154", "cluster_040", "cluster_060", "cluster_075",
  "cluster_081", "cluster_145", "cluster_043", "cluster_088", "cluster_097",
  "cluster_083", "cluster_134", "cluster_167", "cluster_091", "cluster_096",
  "cluster_098", "cluster_111", "cluster_099", "cluster_120", "cluster_122",
  "cluster_126", "cluster_141", "cluster_159", "cluster_116", "cluster_158",
  "cluster_170"
)
result <- unique(top_5(files_list))
print(result)

cluster_order <- c(
    0, 1, 2, 9, 10, 11, 12, 16, 21, 22, 26, 28, 29, 32, 36, 42, 48, 49, 50, 51,
    53, 54, 55, 56, 59, 63, 64, 65, 67, 69, 73, 76, 77, 78, 82, 84, 86, 94, 95,
    101, 103, 104, 105, 107, 109, 110, 113, 114, 117, 118, 123, 124, 127, 128,
    130, 131, 133, 135, 137, 138, 139, 143, 148, 149, 150, 152, 153, 157, 160,
    161, 163, 165, 166, 168, 172, 3, 4, 5, 6, 8, 19, 24, 27, 30, 35, 37, 45,
    47, 52, 57, 61, 62, 71, 79, 85, 89, 90, 93, 100, 102, 108, 112, 119, 121,
    129, 132, 136, 140, 151, 155, 156, 162, 164, 171, 25, 66, 115, 33, 106, 39,
    125, 154, 40, 60, 75, 81, 145, 43, 88, 97, 44, 58, 68, 72, 74, 83, 134, 167,
    91, 96, 98, 111, 99, 120, 122, 126, 141, 159, 116, 158, 170
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

sv_plt("hmgbtfs_top_5_marker_dotplot.pdf", 25, 55)