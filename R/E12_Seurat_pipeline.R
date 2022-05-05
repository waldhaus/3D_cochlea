## E12.5 Seurat scRNA-seq analysis pipeline
library(dplyr)
library(Seurat)
library(tidyr)
library(ggplot2)

# Load dataset (The dataset includes mutant cells and wildtype cells)
data = Read10X(data.dir = "~/Dropbox/Tonotopy_project/Data/aggr_3048_2916_2674/outs/filtered_feature_bc_matrix/")
set.seed(10)

gene_name = rownames(data)
cell_name = colnames(data)
egfp_index = which('Egfp'==gene_name)
eyfp_index = which('Eyfp'==gene_name)
creert2_index = which('Creert2' == gene_name)

# Select E12.5 cells
stage_names = c("E10_Pax2Cre","E12_YFP_Base","E12_YFP_Apex","E14_Apex","E14_Base","E13_Apex","E13_Base")
target_cell_idx = which(grepl('2', cell_name) | grepl('3', cell_name))

Egfp_target = data[egfp_index,target_cell_idx]
Eyfp_target = data[eyfp_index,target_cell_idx]
Creert2_target = data[creert2_index, target_cell_idx]

# Select E12.5 wildtype cells (remove cells with either YFP positive or Creert2 positive)
EGFP_positive_cell = target_cell_idx[which(Egfp_target > 1 & Eyfp_target == 0)]
EGFP_positive_idx = cell_name[EGFP_positive_cell]

target_12 = data[,EGFP_positive_idx]

# Seurat pipeline starts
df_12 = CreateSeuratObject(counts = target_12, project = "corti")
df_12[["percent.mt"]] = PercentageFeatureSet(df_12,pattern = "mt-")
VlnPlot(df_12, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 = FeatureScatter(df_12, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(df_12, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Dimension Reduction
df_12 = subset(df_12, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 8)
df_12 = NormalizeData(df_12, normalization.method = "LogNormalize", scale.factor = 10000)
df_12 = FindVariableFeatures(df_12, selection.method = "vst", nfeatures = 2000)

df_12 = ScaleData(df_12, features = rownames(df_12))
df_12 = RunPCA(df_12, features = VariableFeatures(object = df_12),verbose = FALSE)
DimPlot(df_12, reduction = "pca")

# df_12 = JackStraw(df_12)
# df_12 = ScoreJackStraw(df_12, dims = 1:20)
ElbowPlot(df_12,ndims = 50)

df_12 = FindNeighbors(df_12, dims = 1:30,k.param = 20)
df_12 = FindClusters(df_12, resolution = 0.45)

# Visualization
df_12 = RunUMAP(df_12, dims = 1:12)
DimPlot(df_12, reduction = "umap") 

# Save the E12.5 Seurat object
saveRDS(df_12,"~/Dropbox (University of Michigan)/Tonotopy_project/Data/Tonotopy_analysis/E12_WT_seurat.RDS")

# DE analysis
de_genes = FindAllMarkers(df_12,min.pct = 0.25,logfc.threshold = 0.25,only.pos = T, test.use = "wilcox")
de_genes = subset(de_genes,de_genes$p_val_adj < 0.05)
# write.table(de_genes,"~/Dropbox/Tonotopy_project/Data/Tonotopy_analysis/E12_WT_de_genes_all_clusters.txt",sep = "\t",col.names = T,row.names = F, quote = F)
