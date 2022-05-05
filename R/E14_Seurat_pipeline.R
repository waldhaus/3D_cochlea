## E14.5 Seurat scRNA-seq analysis pipeline
library(dplyr)
library(Seurat)
library(tidyr)
library(ggplot2)

# Load dataset (The dataset includes mutant cells and wildtype cells)
data = Read10X(data.dir = "~/Dropbox (University of Michigan)/Tonotopy_project/Data/aggr_3048_2916_2674/outs/filtered_feature_bc_matrix/")
set.seed(0)

gene_name = rownames(data)
cell_name = colnames(data)
egfp_index = which('Egfp'==gene_name)
eyfp_index = which('Eyfp'==gene_name)
creert2_index = which('Creert2' == gene_name)

# Select E14.5 cells
stage_names = c("E10_Pax2Cre","E12_YFP_Base","E12_YFP_Apex","E14_Apex","E14_Base","E13_Apex","E13_Base")
target_cell_idx = which(grepl('4', cell_name) | grepl('5', cell_name))

Egfp_target = data[egfp_index,target_cell_idx]
Eyfp_target = data[eyfp_index,target_cell_idx]
Creert2_target = data[creert2_index, target_cell_idx]

EGFP_positive_cell = target_cell_idx[which(Egfp_target > 1 & Eyfp_target == 0 & Creert2_target == 0)]
EGFP_positive_idx = (cell_name[EGFP_positive_cell])

target_14 = data[,EGFP_positive_idx]
print(dim(target_14))

# Seurat pipeline starts
df_14 = CreateSeuratObject(counts = target_14, project = "corti")
df_14[["percent.mt"]] = PercentageFeatureSet(df_14,pattern = "mt-")
VlnPlot(df_14, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 = FeatureScatter(df_14, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(df_14, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

df_14 = subset(df_14, subset = nFeature_RNA > 200 & nFeature_RNA < 9500 & percent.mt < 10)
df_14 = NormalizeData(df_14, normalization.method = "LogNormalize", scale.factor = 10000)
df_14 = FindVariableFeatures(df_14, selection.method = "vst", nfeatures = 3000)

all.genes = rownames(df_14)
df_14 = ScaleData(df_14, features = all.genes)

# Dimension Reduction
df_14 = RunPCA(df_14, features = VariableFeatures(object = df_14),verbose = F)
DimPlot(df_14, reduction = "pca")
ElbowPlot(df_14,ndims = 50)

# Identify clusters
df_14 = FindNeighbors(df_14, dims = 1:15,k.param = 20)
df_14 = FindClusters(df_14, resolution = 0.5)

df_14 = RunUMAP(df_14, dims = 1:12)
DimPlot(df_14, reduction = "umap")
saveRDS(df_14,"~/Dropbox/Tonotopy_project/Data/Tonotopy_analysis/E14_WT_seurat.RDS")

# DE analysis
de_genes = FindAllMarkers(df_14,min.pct = 0.25,logfc.threshold = 0.25,only.pos = T, test.use = "wilcox")
de_genes = subset(de_genes,de_genes$p_val_adj < 0.05)
# write.table(de_genes,"~/Dropbox/Tonotopy_project/Data/Tonotopy_analysis/E14_WT_de_genes_all_clusters.txt",sep = "\t",col.names = T,row.names = F, quote = F)
