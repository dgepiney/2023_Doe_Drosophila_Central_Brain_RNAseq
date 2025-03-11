#------------------------------------------------------------------------------- 
# Preliminaries:
# Close graphics
graphics.off()
# Clear R memory
rm(list=ls()) 
# Set directory 
setwd("/Users/noahdillon/University of Oregon Dropbox/Noah Dillon/23 T2 snRNAseq/_new T2rnaseq paper/Data Plots/T2_glia_and_neuron_sex_DE_NRD/") 
#------------------------------------------------------------------------------- 

# Load libraries 
{
  library(tidyverse)
  library(Seurat)
  library(patchwork)
  library(Matrix)
  library(cowplot)
  library(viridis)
  library(ggplot2)
  library(ggrepel)
  library(BiocManager)
  library(multtest)
  library(metap)
  library(DESeq2)
  library(Scillus)
  library(scCustomize)
}

load("T2.atlas.RData")

# Subset and re-cluster glial cells (Repo+)
{
  # Identify Repo+ clusters 
  T2.atlas.all_repo <- FindAllMarkers(object = T2.atlas, assay = "RNA", features = "repo",
                                      logfc.threshold = 0.1, test.use = "wilcox", min.pct = 0.2,
                                      verbose = TRUE, only.pos = TRUE, return.thresh = 0.01)
  write.csv(T2.atlas.all_repo, file = "DE_lists/T2_atlas_all_repo.csv")
  df <- read.csv("DE_lists/T2_atlas_all_repo.csv")
  repo_clusters.list <- as.list(df["cluster"])
  
  # Subset and re-cluster 
  T2.atlas.glia = subset(T2.atlas, idents = repo_clusters.list$cluster)
  T2.atlas.glia <- RunUMAP(T2.atlas.glia, reduction = 'integrated.rpca', dims = 1:50,
                           reduction.name='umap.rpca')
  T2.atlas.glia <- FindNeighbors(T2.atlas.glia, reduction = 'integrated.rpca', dims = 1:50)
  T2.atlas.glia <- FindClusters(T2.atlas.glia, resolution = .2) 
  
  # saveRDS(T2.atlas.glia, file = "RDS/T2_glia_unlabeled.rds")
  # read RDS file
  T2.atlas.glia <- readRDS("RDS/T2_glia_unlabeled.rds")
  
  # Exploration 
  p1 <- DimPlot(T2.atlas.glia, reduction='umap.rpca', label = T, repel = T) # Clusters 
  
  # Sub type identification of DE genes
  df <- read.csv("Marker_lists/glia_subtype_markers_filtered.csv")
  glia_markers.list <- as.list(df["Gene"])
  
  # Find DE genes 
  glia_markers <- FindAllMarkers(object = T2.atlas.glia, test.use = "wilcox", 
                                 features = glia_markers.list$Gene, only.pos = TRUE, 	
                                 return.thresh = 0.01, logfc.threshold = 0.1, 
                                 min.pct = 0.25, verbose = TRUE)
  write.csv(glia_markers, file = "DE_lists/glia_markers_DE.csv")
  glia_markers_filtered <- glia_markers[glia_markers$p_val_adj < 0.05, ]
  write.csv(glia_markers_filtered, file = "DE_lists/glia_markers_DE_filtered.csv")
  write.csv(glia_markers_filtered, file = "Supp_tables/Supp_table_glia_markers.csv")
  glia_markers_filtered <- glia_markers_filtered[!duplicated(glia_markers_filtered$gene), ]
  
  # Identify cluster sub type ident
  {
    T2.atlas.glia_labeled <- T2.atlas.glia
    new.cluster.ids <- c("Ensheathing - 0", "Ensheathing/Cortex - 1", "Astrocyte - 2", "Astrocyte-like - 3",
                         "Perineurial - 4", "Chiasm - 5", "Subperineurial - 6", "Ensheathing - 7",
                         "Astrocyte - 8", "Ensheathing - 9", "Ensheathing - 10", 
                         "Ensheathing/Astrocyte - 11", "Chiasm - 12", "Astrocyte - 13")
    names(new.cluster.ids) <- levels(T2.atlas.glia_labeled)
    T2.atlas.glia_labeled <- RenameIdents(T2.atlas.glia_labeled, new.cluster.ids)
    cluster_ident.list <- c("Ensheathing - 0", "Ensheathing - 10", "Ensheathing - 9",
                            "Ensheathing/Cortex - 1", "Ensheathing - 7", "Ensheathing/Astrocyte - 11",
                            "Astrocyte - 2", "Astrocyte - 8", "Astrocyte - 13", "Astrocyte-like - 3",
                            "Chiasm - 12", "Chiasm - 5", "Perineurial - 4", "Subperineurial - 6")
    T2.atlas.glia_labeled@active.ident <- factor(x = T2.atlas.glia_labeled@active.ident, levels = cluster_ident.list)
    T2.atlas.glia_labeled$CellType <- Idents(T2.atlas.glia_labeled)
    # T2.atlas.glia_labeled <- saveRDS(T2.atlas.glia_labeled, file = "RDS/T2_glia_labeled.rds")
    
    glia_DE <- read.csv("DE_lists/glia_markers_DE_filtered.csv")
    glia_markers_filtered <- read.csv("Marker_lists/glia_subtype_markers_filtered.csv")
    glia_markers_filtered <- glia_DE[glia_DE$gene %in% glia_markers_filtered$Gene, ]
    glia_markers_filtered <- glia_markers_filtered[!duplicated(glia_markers_filtered$gene), ]
    glia_markers_filtered <- as.list(glia_markers_filtered["gene"])
    
    ## Figure panels
    p1 <- DimPlot(T2.atlas.glia_labeled, reduction='umap.rpca', label = F, repel = T) # Clusters 
    ggsave(filename = "Plots/Glia_clusters.pdf", plot = p1, device = "pdf",
           scale = 1, width = 18, height = 10, units = c("cm"), dpi = 300, limitsize = TRUE)
    
    p1 <- DotPlot(T2.atlas.glia_labeled, features = glia_markers_filtered$gene, cols = c("white", "black"),
                  col.min = 0, cluster.idents = F, dot.min = 0, scale.by = "size", dot.scale = 6) + RotatedAxis()
    ggsave(filename = "Plots/Glia_markers_dotplot_labeled.pdf", plot = p1, device = "pdf",
           scale = 1, width = 25, height = 12, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
  
  # saveRDS(T2.atlas.glia_labeled, file = "RDS/T2_glia_labeled.rds")
  # read RDS file
  T2.atlas.glia_labeled <- readRDS("RDS/T2_glia_labeled.rds")
}

# Explore Sex differences in glia
{
  # number of nuclei for each sample
  table(T2.atlas.glia_labeled@meta.data$sample)
  
  # Female vs Male
  {
    cluster.population = as.data.frame.matrix(table(T2.atlas.glia$seurat_clusters, T2.atlas.glia$sample))
    n_male <- length(WhichCells(T2.atlas.glia, expression = sample == "T2-male"))
    n_female <- length(WhichCells(T2.atlas.glia, expression = sample == "T2-female"))
    global.ratio <- n_male / n_female
    
    cluster.population = cluster.population %>% mutate( male.adjusted = `T2-male` / global.ratio)
    
    cluster.ratios <- cluster.population %>% mutate(ratio = male.adjusted / (`T2-female` + 1)) 
    cluster.ratios <- cluster.ratios %>%
      mutate(log2.ratio = log2(ratio),
             p.value = pnorm(abs(log2.ratio), lower.tail = F) * 2)
    
    ## Figure panel
    p1 <- ggplot(cluster.ratios, aes(x = log2.ratio, y = -log10(p.value))) + 
      geom_point(aes(color = p.value < 0.05), size = 3) + ggtitle("Male vs female biased clusters") +
      geom_text_repel(label = ifelse(cluster.ratios$p.value < 0.05, rownames(cluster.ratios),""), box.padding = 0.75, size = 6) +
      theme_classic() + labs(x = "Log2(Male/Female Ratio)", y = "-log10(p-value)", color = "Biased Cluster") 
    
    ggsave(filename = "Plots/Male_female_enrichment.tiff", plot = p1, device = "tiff",
           scale = 1, width = 8, height = 8, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
  
  # Visualize Sex differences in clusters 
  ## Figure panel
  p1 <- DimPlot(T2.atlas.glia_labeled, reduction='umap.rpca', label = F, split.by = "sample")
  ggsave(filename = "Plots/Glia_clusters_sample.pdf", plot = p1, device = "pdf",
         scale = 1, width = 32, height = 10, units = c("cm"), dpi = 300, limitsize = TRUE)
  
  # pseudo-bulk and DE between all male and female nuclei
  {
    bulk <- AggregateExpression(T2.atlas.glia, group.by = c("sample", "orig.ident"), return.seurat = T)
    Idents(bulk) <- "sample"
    de_markers <- FindMarkers(bulk, ident.1 = "T2-male", ident.2 = "T2-female", 
                              slot = "counts", test.use = "DESeq2", verbose = T)
    de_markers$gene <- rownames(de_markers)
    write.csv(de_markers, file = "Supp_tables/Supp_table_glia_sex_DE.csv")
    p1 <- ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.4, alpha = 0.5) + 
      theme_classic() + ylab("-log10(unadjusted p-value)") + xlab("avg log2FC Male vs Female") +
      geom_point(aes(colour=cut(p_val_adj, c(-Inf, 0.01, Inf)))) + theme(legend.position="none") +
      scale_color_manual(name = "p_val_adj", values = c("(-Inf,0.01]" = "magenta","(0.01,Inf]" = "black")) +
      geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,"")), colour = "black", max.overlaps = 6,
                      size = 2) + 
      labs(title= "Male vs Female expression")
    ggsave(filename = "Plots/Sex_DE/sex_differences_all_.pdf", plot = p1, device = "pdf",
           scale = 1, width = 10, height = 10, units = c("cm"), dpi = 300, limitsize = TRUE)
  }
  
  # Heat map of top genes in sex-enriched clusters 
  {
    # Subset male and female glia
    T2.atlas.glia_M_F = subset(T2.atlas.glia, subset = sample == c("T2-male","T2-female"))
    
    # Find significant clusters 
    cluster.population = as.data.frame.matrix(table(T2.atlas.glia$seurat_clusters, T2.atlas.glia$sample))
    n_male <- length(WhichCells(T2.atlas.glia, expression = sample == "T2-male"))
    n_female <- length(WhichCells(T2.atlas.glia, expression = sample == "T2-female"))
    global.ratio <- n_male / n_female
    
    cluster.population = cluster.population %>% mutate( male.adjusted = `T2-male` / global.ratio)
    
    cluster.ratios <- cluster.population %>% mutate(ratio = male.adjusted / (`T2-female` + 1)) 
    cluster.ratios <- cluster.ratios %>%
      mutate(log2.ratio = log2(ratio),
             p.value = pnorm(abs(log2.ratio), lower.tail = F) * 2)
    cluster.ratios <- cluster.ratios[cluster.ratios$p.value > 0, ]
    
    signif_clusters <- cluster.ratios[cluster.ratios$p.value < 0.05, ]
    signif_clusters <- rownames(signif_clusters)
    
    # Subset significant clusters 
    T2.atlas.glia_M_F <- subset(T2.atlas.glia_M_F, idents = signif_clusters)
    
    # Find significantly DE genes within clusters 
    bulk <- AggregateExpression(T2.atlas.glia_M_F, group.by = c("sample", "seurat_clusters", "orig.ident"), 
                                return.seurat = T)
    de_markers_all = data.frame()
    for(i in signif_clusters) {
      x <- subset(bulk, seurat_clusters == i)
      Idents(x) <- "sample"
      de_markers <- FindMarkers(x, ident.1 = "T2-male", ident.2 = "T2-female", 
                                slot = "counts", test.use = "DESeq2", verbose = T)
      de_markers <- de_markers[de_markers$p_val_adj < 0.05, ]
      de_markers$gene <- rownames(de_markers)
      de_markers_all <- rbind(de_markers_all,de_markers)
    }
    
    # Remove duplicated values 
    de_markers_unique <- de_markers_all[!duplicated(de_markers_all$gene), ]
    de_markers_unique <- na.omit(de_markers_unique)
    
    # Heat map of DE between M and F in clusters 
    color_cluster <- DiscretePalette_scCustomize(num_colors = 3, palette = "varibow")
    p1 <- plot_heatmap(dataset = T2.atlas.glia_M_F, markers = de_markers_unique$gene,
                       sort_var = c("seurat_clusters","sample"), anno_var = c("seurat_clusters","sample"),
                       anno_colors = list(color_cluster,c("green","magenta")), 
                       hm_colors = c("#440154","#21918c","#fde725"))
    # MANUALLY SAVE HEATMAP
  }
}

# T2 neuron sex differences
{
  # number of nuclei for each sample
  table(T2.atlas.neurons@meta.data$sample)
  
  # Female vs Male
  {
    cluster.population = as.data.frame.matrix(table(T2.atlas.neurons$seurat_clusters, T2.atlas.neurons$sample))
    n_male <- length(WhichCells(T2.atlas.neurons, expression = sample == "T2-male"))
    n_female <- length(WhichCells(T2.atlas.neurons, expression = sample == "T2-female"))
    global.ratio <- n_male / n_female
    
    cluster.population = cluster.population %>% mutate( male.adjusted = `T2-male` / global.ratio)
    
    cluster.ratios <- cluster.population %>% mutate(ratio = male.adjusted / (`T2-female` + 1)) 
    cluster.ratios <- cluster.ratios %>%
      mutate(log2.ratio = log2(ratio),
             p.value = pnorm(abs(log2.ratio), lower.tail = F) * 2)
    cluster.ratios <- cluster.ratios[cluster.ratios$p.value > 0, ]
    
    ## Figure panel
    p1 <- ggplot(cluster.ratios, aes(x = log2.ratio, y = -log10(p.value))) + 
      geom_point(aes(color = p.value < 0.05), size = 3) + ggtitle("Male vs female biased clusters") +
      geom_text_repel(label = ifelse(cluster.ratios$p.value < 0.05, rownames(cluster.ratios),""), 
                      box.padding = 0.75, size = 5, max.overlaps = 9) +
      theme_classic() + labs(x = "Log2(Male/Female Ratio)", y = "-log10(p-value)", color = "Biased Cluster") 
    
    ggsave(filename = "Plots/Male_female_enrichment_neurons.tiff", plot = p1, device = "tiff",
           scale = 1, width = 8, height = 8, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
  
  # Visualize Sex differences in clusters 
  ## Figure panel
  p1 <- DimPlot(T2.atlas.neurons, reduction='umap.rpca', label = F, split.by = "sample") + NoLegend()
  ggsave(filename = "Plots/Glia_clusters_sample_neurons.tiff", plot = p1, device = "tiff",
         scale = 1, width = 28, height = 10, units = c("cm"), dpi = 200, limitsize = TRUE)
  
  # pseudo-bulk and DE between all male and female nuclei
  {
    bulk <- AggregateExpression(T2.atlas.neurons, group.by = c("sample", "orig.ident"), return.seurat = T)
    Idents(bulk) <- "sample"
    de_markers <- FindMarkers(bulk, ident.1 = "T2-male", ident.2 = "T2-female", 
                              slot = "counts", test.use = "DESeq2", verbose = T)
    de_markers$gene <- rownames(de_markers)
    write.csv(de_markers, file = "Supp_tables/Supp_table_neuron_sex_DE.csv")
    p1 <- ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.4, alpha = 0.5) + 
      theme_classic() + ylab("-log10(unadjusted p-value)") + xlab("avg log2FC Male vs Female") +
      geom_point(aes(colour=cut(p_val_adj, c(-Inf, 0.01, Inf)))) + theme(legend.position="none") +
      scale_color_manual(name = "p_val_adj", values = c("(-Inf,0.01]" = "magenta","(0.01,Inf]" = "black")) +
      geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,"")), colour = "black", max.overlaps = 4, 
                      box.padding	= 0.1 , point.padding = 0.2, size = 2) + 
      labs(title= "Male vs Female expression")
    ggsave(filename = "Plots/Sex_DE/sex_differences_all_neurons.pdf", plot = p1, device = "pdf",
           scale = 1, width = 10, height = 10, units = c("cm"), dpi = 300, limitsize = TRUE)
  }
  
  # Heat map of top genes in sex-enriched clusters 
  {
    # Subset male and female neurons
    T2.atlas.neuron_M_F = subset(T2.atlas.neurons, subset = sample == c("T2-male","T2-female"))
    
    # Find significant clusters 
    cluster.population = as.data.frame.matrix(table(T2.atlas.neurons$seurat_clusters, T2.atlas.neurons$sample))
    n_male <- length(WhichCells(T2.atlas.neurons, expression = sample == "T2-male"))
    n_female <- length(WhichCells(T2.atlas.neurons, expression = sample == "T2-female"))
    global.ratio <- n_male / n_female
    
    cluster.population = cluster.population %>% mutate( male.adjusted = `T2-male` / global.ratio)
    
    cluster.ratios <- cluster.population %>% mutate(ratio = male.adjusted / (`T2-female` + 1)) 
    cluster.ratios <- cluster.ratios %>%
      mutate(log2.ratio = log2(ratio),
             p.value = pnorm(abs(log2.ratio), lower.tail = F) * 2)
    cluster.ratios <- cluster.ratios[cluster.ratios$p.value > 0, ]
    
    signif_clusters <- cluster.ratios[cluster.ratios$p.value < 0.05, ]
    signif_clusters <- rownames(signif_clusters)
    
    # Subset significant clusters 
    T2.atlas.neuron_M_F <- subset(T2.atlas.neuron_M_F, idents = signif_clusters)
    
    # Find significantly DE genes within clusters 
    bulk <- AggregateExpression(T2.atlas.neuron_M_F, group.by = c("sample", "seurat_clusters", "orig.ident"), 
                                return.seurat = T)
    de_markers_all = data.frame()
    for(i in signif_clusters) {
      x <- subset(bulk, seurat_clusters == i)
      Idents(x) <- "sample"
      de_markers <- FindMarkers(x, ident.1 = "T2-male", ident.2 = "T2-female", 
                                slot = "counts", test.use = "DESeq2", verbose = T)
      de_markers <- de_markers[de_markers$p_val_adj < 0.05, ]
      de_markers$gene <- rownames(de_markers)
      de_markers_all <- rbind(de_markers_all,de_markers)
    }
    
    # Remove duplicated values 
    de_markers_unique <- de_markers_all[!duplicated(de_markers_all$gene), ]
    de_markers_unique <- na.omit(de_markers_unique)
    
    # Heat map of DE between M and F in clusters 
    color_cluster <- DiscretePalette_scCustomize(num_colors = 14, palette = "varibow")
    p1 <- plot_heatmap(dataset = T2.atlas.neuron_M_F, markers = de_markers_unique$gene,
                       sort_var = c("seurat_clusters","sample"), anno_var = c("seurat_clusters","sample"),
                       anno_colors = list(color_cluster,c("green","magenta")), 
                       hm_colors = c("#440154","#21918c","#fde725"))
    # MANUALLY SAVE HEATMAP
  }
}