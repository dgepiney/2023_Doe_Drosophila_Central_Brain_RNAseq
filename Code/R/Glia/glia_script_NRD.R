# scRNA-seq analysis of adult central brain
# Exploring differences in T1 and T2 glia
# Author: Noah R Dillon 
# Date Oct-27-2023

# Load libraries 
{
  library(dplyr) # version 1.0.7
  library(Seurat) # version 4.0.4
  library(patchwork) # version 1.1.1
  library(glmGamPoi) # version 1.2.0
  library(ggplot2) # version 3.3.5
  library(sctransform) # version 0.3.2
  library(cowplot) # version 1.1.1
  library(plyr) # version 1.8.6
}

# Load whole brain atlas
whole_brain <- load("Whole_brain_atlas/unsorted.integrated.RData")

# subset and re-cluster glia cell types
{
  # Markers to start sub-setting clusters 
  whole_atlas.all_repo <- FindAllMarkers(object = unsorted.integrated, assay = "RNA", features = "repo",
                                         logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.1,
                                         min.diff.pct = 0.1,verbose = TRUE, only.pos = TRUE, return.thresh = 0.01)
  write.csv(whole_atlas.all_repo, file = "csv_files/whole_atlas_all_repo.csv")
  df <- read.csv("csv_files/whole_atlas_all_repo.csv")
  repo_clusters.list <- as.list(df["cluster"])
  
  # Subset and re-cluster
  repo_unlabeled <- subset(unsorted.integrated, idents = repo_clusters.list$cluster)
  DefaultAssay(repo_unlabeled) <- "integrated"
  repo_unlabeled <- RunUMAP(repo_unlabeled, reduction = "pca", dims = 1:50)
  repo_unlabeled <- FindNeighbors(repo_unlabeled, reduction = "pca", dims = 1:50)
  repo_unlabeled <- FindClusters(repo_unlabeled, resolution = .3) 
  DefaultAssay(repo_unlabeled) <- "RNA"
  
  # find cluster defining genes
  glia_markers.all <- FindAllMarkers(object = repo_unlabeled, assay = "RNA",
                                     logfc.threshold = 0.1,test.use = "wilcox", 
                                     min.pct = 0.1, min.diff.pct = 0.1,
                                     verbose = TRUE, only.pos = TRUE, 
                                     return.thresh = 0.001)
  write.csv(glia_markers.all, file = "csv_files/glia_subset_all_markers.csv")
}

# Assigning cluster identities 
{
  repo_labeled <- repo_unlabeled
  new.cluster.ids <- c("Cortex - 0", "Perineurial - 1", "Astrocyte - 2", "Ensheathing - 3",
                       "Ensheathing - 4", "Surface - 5", "Astrocyte - 6", "Subperineurial - 7",
                       "Subperineurial - 8")
  names(new.cluster.ids) <- levels(repo_labeled)
  repo_labeled <- RenameIdents(repo_labeled, new.cluster.ids)
  cluster_ident.list <- c("Surface - 5", "Subperineurial - 7","Subperineurial - 8",
                          "Perineurial - 1", "Ensheathing - 4", "Ensheathing - 3",
                          "Cortex - 0","Astrocyte - 2","Astrocyte - 6")
  repo_labeled@active.ident <- factor(x = repo_labeled@active.ident, levels = cluster_ident.list)
  
  p1 <- DimPlot(repo_labeled, reduction = "umap", label = FALSE) + NoAxes() +NoLegend()
  
  ggsave(filename = "plots/glia_subtype_UMAP.svg", plot = p1, device = "svg",
         scale = 1, width = 20, height = 20, units = c("cm"), dpi = 300, 
         limitsize = TRUE, bg = "white")
  
  df <- read.csv("glia_subtype_markers_curated.csv")
  glia_subtype_curated_markers.list <- as.list(df["Gene"])
  p2 <- DotPlot(repo_labeled, features = glia_subtype_curated_markers.list$Gene, 
                cols = c("white", "black"),col.min = 0, cluster.idents = FALSE, dot.min = 0, 
                scale.by = "size", dot.scale = 6) + RotatedAxis()
  
  ggsave(filename = "plots/glia_subtype_dotplot.svg", plot = p2, device = "svg",
         scale = 1, width = 25, height = 12, units = c("cm"), dpi = 300, 
         limitsize = TRUE, bg = "white")
}

# Comparing T1 vs T2 glia subtypes
{
  # plot T1 and T2 glia UMAPs
  {
    p1<- DimPlot(repo_labeled, group.by = 'Lineage', cols=c("cyan", "dark red"), 
            order = c("T2"), raster=FALSE) + NoAxes() +NoLegend()
    
    ggsave(filename = "plots/T1vsT2_glia_umap.svg", plot = p1, device = "svg",
           scale = 1, width = 20, height = 20, units = c("cm"), dpi = 300, 
           limitsize = TRUE, bg = "white")
  }
  
  # find differential expression of genes between T1 and T2 
  {
    # Do T1 and T2 have conserved glial marker expression?
    Idents(repo_labeled) <- factor(Idents(repo_labeled), levels = cluster_ident.list)
    p1 <- DotPlot(repo_labeled, features = rev(glia_subtype_curated_markers.list), cols = c("cyan", "dark red"), dot.scale = 6, 
            split.by = "Lineage") + RotatedAxis()
    
    ggsave(filename = "plots/T1vsT2_glia_dotplot.svg", plot = p1, device = "svg",
           scale = 1, width = 25, height = 20, units = c("cm"), dpi = 300, 
           limitsize = TRUE, bg = "white")
    
    # differential expression of genes 
    repo_unlabeled_test <- repo_unlabeled
    repo_unlabeled_test$celltype.Lineage <- paste(Idents(repo_unlabeled_test), 
                                                  repo_unlabeled_test$Lineage, sep = "_")
    repo_unlabeled_test$celltype <- Idents(repo_unlabeled_test)
    Idents(repo_unlabeled_test) <- "celltype.Lineage"
    for (i in 1:9){
      cluster = i-1 
      T1.vs.T2 <- FindMarkers(repo_unlabeled_test, assay = "RNA", ident.1 = (sprintf("%d_T1", cluster)),
                              ident.2 = (sprintf("%d_T2", cluster)), logfc.threshold = 0.1,test.use = "wilcox", 
                              min.pct = 0.1, min.diff.pct = 0.1, only.pos = TRUE, return.thresh = 0.01, verbose = TRUE)
      write.csv(T1.vs.T2, file = (sprintf("csv_files/T1_vs_T2_cluster_%d.csv", cluster)))
      i = +1
    }
  }
}

# Glia in T2
{
  # subset out repo+ clusters
  {
    # Markers to start sub-setting clusters 
    T2_atlas.all_repo <- FindAllMarkers(object = T2_all.integrated, assay = "RNA", features = "repo",
                                        logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.1,
                                        min.diff.pct = 0.1,verbose = TRUE, only.pos = TRUE, return.thresh = 0.01)
    write.csv(T2_atlas.all_repo, file = "csv_files/T2_atlas_all_repo.csv")
    df <- read.csv("csv_files/T2_atlas_all_repo.csv")
    repo_clusters.list <- as.list(df["cluster"])
    
    # Subset and re-cluster
    repo_unlabeled <- subset(T2_all.integrated, idents = repo_clusters.list$cluster)
    
    DefaultAssay(repo_unlabeled) <- "integrated"
    repo_unlabeled <- RunUMAP(repo_unlabeled, reduction = "pca", dims = 1:50)
    repo_unlabeled <- FindNeighbors(repo_unlabeled, reduction = "pca", dims = 1:50)
    repo_unlabeled <- FindClusters(repo_unlabeled, resolution = .3) 
    DefaultAssay(repo_unlabeled) <- "RNA"
    
    DimPlot(repo_unlabeled, reduction = "umap", label = TRUE) + NoLegend() + NoAxes()
    
    df <- read.csv("csv_files/glia_subtype_markers_curated.csv")
    glia_subtype_curated_markers.list <- as.list(df["Gene"])
    
    DotPlot(repo_unlabeled, features = glia_subtype_curated_markers.list$Gene, 
            cols = c("white", "black"),col.min = 0, cluster.idents = FALSE, dot.min = 0, 
            scale.by = "size", dot.scale = 6) + RotatedAxis()
    
    glia_markers.all <- FindAllMarkers(object = repo_unlabeled, assay = "RNA",
                                       logfc.threshold = 0.1,test.use = "wilcox", 
                                       min.pct = 0.1, min.diff.pct = 0.1,
                                       verbose = TRUE, only.pos = TRUE, 
                                       return.thresh = 0.001)
    write.csv(glia_markers.all, file = "csv_files/glia_subset_all_markers.csv")
  }
  
  # Assigning cluster identities 
  {
    repo_labeled <- repo_unlabeled
    new.cluster.ids <- c("Neurons/glia mix - 0", "Ensheathing - 1", "Cortex - 2", "Perineurial - 3",
                         "Subperineurial - 4", "Astrocyte - 5", "Neurons/glia mix - 6")
    names(new.cluster.ids) <- levels(repo_labeled)
    repo_labeled <- RenameIdents(repo_labeled, new.cluster.ids)
    cluster_ident.list <- c("Neurons/glia mix - 6", "Neurons/glia mix - 0",
                            "Subperineurial - 4", "Perineurial - 3", "Ensheathing - 1",
                            "Cortex - 2", "Astrocyte - 5")
    repo_labeled@active.ident <- factor(x = repo_labeled@active.ident, levels = cluster_ident.list)
    
    p1 <- DimPlot(repo_labeled, reduction = "umap", label = FALSE) + NoAxes() + NoLegend()
    
    ggsave(filename = "plots/glia_subtype_UMAP.svg", plot = p1, device = "svg",
           scale = 1, width = 25, height = 20, units = c("cm"), dpi = 300, 
           limitsize = TRUE, bg = "white")
    
    p2 <- DimPlot(repo_labeled, group.by = 'dataset', cols=c("cyan", "cyan", "cyan","dark red"), 
                  order = c("parse-unsorted", "parse-redstinger+","parse-unc84sfGFP+", "pipseq"), 
                  raster=FALSE) + NoAxes() +NoLegend()
    
    ggsave(filename = "plots/glia_subtype_UMAP_data_set.svg", plot = p2, device = "svg",
           scale = 1, width = 20, height = 20, units = c("cm"), dpi = 300, 
           limitsize = TRUE, bg = "white")
    
    df <- read.csv("csv_files/glia_subtype_markers_curated.csv")
    glia_subtype_curated_markers.list <- as.list(df["Gene"])
    
    p3 <- DotPlot(repo_labeled, features = glia_subtype_curated_markers.list$Gene, 
                  cols = c("white", "black"),col.min = 0, cluster.idents = FALSE, dot.min = 0, 
                  scale.by = "size", dot.scale = 6) + RotatedAxis()
    
    ggsave(filename = "plots/glia_subtype_dotplot.svg", plot = p3, device = "svg",
           scale = 1, width = 25, height = 12, units = c("cm"), dpi = 300, 
           limitsize = TRUE, bg = "white")
    
    saveRDS(repo_labeled, file = "T2_only_atlas/repo_labeled.rds")
  }
}
