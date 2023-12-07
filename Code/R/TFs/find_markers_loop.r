library(Seurat)

# For loop to find markers in each cluster
marker_loop <- function(genes){
    for(i in 148:173){
        temp_markers <<- FindMarkers(
            sub_my_data,
            ident.1 = i,
            features = intersect(rownames(my_data), genes$SYMBOL),
            min.pct = 0.25,
            only.pos = TRUE 
        )
        write.csv(temp_markers, sprintf("/Users/derek/Dropbox (University of Oregon)/lab_stuffs/r_data/cluster_%03d.csv", i))
    }
}

