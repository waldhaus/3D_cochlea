## Use LIGER to align E12.5 and E14.5 datasets to validate the Seurat clustering 
library(rliger)
library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)

## Sec1: Align E12 and E14 datasets only 
E12_seurat = readRDS("~/Dropbox/Tonotopy_project/Data/Tonotopy_analysis_rerun/E12_WT_seurat_rerun.RDS")
E14_seurat = readRDS("~/Dropbox/Tonotopy_project/Data/Tonotopy_analysis_rerun/E14_WT_seurat_rerun.RDS")

E12_seurat@meta.data$seurat_anno3 = E12_seurat@meta.data$seurat_clusters
E12_seurat@meta.data$seurat_anno3 = factor(E12_seurat@meta.data$seurat_clusters,levels = 0:8,labels = c(
  "CochleaApex","CochleaRoof","HB","SGN","CochleaBase","HB","Schwann","SGN","Mes"))
E14_seurat@meta.data$seurat_anno3 = E14_seurat@meta.data$seurat_clusters
E14_seurat@meta.data$seurat_anno3 = factor(E14_seurat@meta.data$seurat_clusters,levels = 0:9,labels = c(
  "Schwann","CochleaRoof","CochleaBase","CochleaApex","HB","Schwann","HB","SGN","Mes","Keratin"))

# Select data
E12_data = E12_seurat@assays$RNA@data
E14_data = E14_seurat@assays$RNA@data

e1214_liger = createLiger(list(E12 = E12_data, E14 = E14_data))

# Start LIGER pipeline
e1214_liger@norm.data$E12 = e1214_liger@raw.data$E12
e1214_liger@norm.data$E14 = e1214_liger@raw.data$E14
e1214_liger = selectGenes(e1214_liger)
e1214_liger = scaleNotCenter(e1214_liger)

e1214_liger = optimizeALS(e1214_liger, k = 20,rand.seed = 1)
e1214_liger = quantile_norm(e1214_liger)
e1214_liger = louvainCluster(e1214_liger, resolution = 0.25)
e1214_liger = runUMAP(e1214_liger, distance = 'cosine', n_neighbors = 20, min_dist = 0.3)
liger_plot = plotByDatasetAndCluster(e1214_liger, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = T)
liger_plot[[1]]
liger_plot[[2]]

# Save the LIGER object
saveRDS(e1214_liger,"~/Dropbox/Tonotopy_project/Data/Tonotopy_analysis/E12_E14_LIGER_alignment.RDS")
