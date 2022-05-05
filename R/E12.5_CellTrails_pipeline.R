## E12.5 CellTrails subpopulation identification
library(CellTrails)
library(dplyr)
library(Seurat)
library(ggplot2)

# Load the data 
e12_seurat = readRDS("~/Dropbox/Tonotopy_project/Data/Tonotopy_analysis/E12_WT_seurat.RDS")
e12_seurat_cochlear = SubsetData(e12_seurat,ident.use=c(0,1,4))
E12 = GetAssayData(object = e12_seurat_cochlear, slot = 'data')

index = read.csv("~/Dropbox/Tonotopy_project/Data/data/gene_index_no_cell_cycle_with_1_starting.csv",header = FALSE)
index = index$V1
seurat_cluster = e12_seurat_cochlear@active.ident

# Generate SCE object to prepare for celltrails
E12 = as.matrix(E12)[index,]
E12 = SingleCellExperiment(assays=list(logcounts = E12))
isSpike(E12, "ERCC") = 1:80

# Removes features that are not expressed or that do not sufficiently reach the technological limit of detection
trajFeatureNames(E12) = filterTrajFeaturesByDL(E12, threshold=2, show_plot=T)
trajFeatureNames(E12) = filterTrajFeaturesByCOV(E12, threshold=0.5, show_plot=T)
E12_sub = E12[trajFeatureNames(E12), ]

# Filter using scran
var_fit = scran::trendVar(x=E12_sub, use.spikes=FALSE)
var_out = scran::decomposeVar(x=E12_sub, fit=var_fit)
tfeat = featureNames(E12_sub)[which(var_out$FDR < 0.01)]
trajFeatureNames(E12) = tfeat

# Dimension reduction 
se = embedSamples(E12)
d = findSpectrum(se$eigenvalues, frac=100)
latentSpace(E12) = se$components[, d]

# Find states
cl = findStates(E12, min_size=0.005, min_feat=2, max_pval=1e-4, min_fc=1)
states(E12) = cl
plotManifold(E12, color_by="phenoName", name="state")

E12$seurat_cluster = seurat_cluster
plotManifold(E12, color_by="phenoName", name="seurat_cluster")

# Sample ordering
set.seed(10)
E12 = connectStates(E12, l=15)
plotStateTrajectory(E12, color_by="phenoName", name="seurat_cluster", component=1, point_size=1.5, label_offset=4)

# Make trajectory
E12 = selectTrajectory(E12, component=1)
E12 = fitTrajectory(E12)
plotTrajectoryFit(E12) 

# write.ygraphml(sce=E12,file='~/Dropbox/Tonotopy_project/Data/Tonotopy_analysis/E12_trajectory3.graphml',color_by='phenoName',name='state',node_label='state')

tl = read.ygraphml("~/Dropbox/Tonotopy_project/Data/Tonotopy_analysis/E12_trajectory3.graphml")
plot(tl[,1:2], axes=FALSE, xlab="", ylab="", pch=20, cex=.25)
trajLayout(E12, adjust=TRUE) = tl
