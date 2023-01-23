## E14.5 CellTrails subpopulation identification
library(CellTrails)
library(dplyr)
library(Seurat)
library(ggplot2)

# Load the data 
e14_seurat = readRDS("~/Dropbox/Tonotopy_project/Data/Tonotopy_analysis_rerun/E14_WT_seurat_rerun.RDS")
e14_seurat_cochlear = subset(e14_seurat,idents=c(1,2,3))
E14 = GetAssayData(object = e14_seurat_cochlear, slot = 'data')

index = read.csv("~/Dropbox/Tonotopy_project/Data/data/gene_index_no_cell_cycle_with_1_starting.csv",header = T)
index = index$X0
seurat_cluster = e14_seurat_cochlear@active.ident

# Generate SCE object to prepare for celltrails
E14 = as.matrix(E14)[index,]
E14 = SingleCellExperiment(assays=list(logcounts=E14))
isSpike(E14, "ERCC") = 1:80

# Removes features that are not expressed or that do not sufficiently reach the technological limit of detection
trajFeatureNames(E14) = filterTrajFeaturesByDL(E14, threshold=2, show_plot=FALSE)
trajFeatureNames(E14) = filterTrajFeaturesByCOV(E14, threshold=0.5, show_plot=FALSE)
E14_sub = E14[trajFeatureNames(E14), ]

# Filter using scran
var_fit = scran::trendVar(x=E14_sub, use.spikes=FALSE)
var_out = scran::decomposeVar(x=E14_sub, fit=var_fit)
tfeat = featureNames(E14_sub)[which(var_out$FDR < 0.01)]
trajFeatureNames(E14) = tfeat

# Dimension reduction 
se = embedSamples(E14)
d = findSpectrum(se$eigenvalues, frac=100)
latentSpace(E14) = se$components[, d]

# Find states
cl = findStates(E14, min_size=0.01, min_feat=2, max_pval=1e-4, min_fc=1.5)

states(E14) = cl
plotManifold(E14, color_by="phenoName", name="state")

E14$Experiment = seurat_cluster
plotManifold(E14, color_by="phenoName", name="Experiment")

# Sample ordering
E14 = connectStates(E14, l=8)
plotStateTrajectory(E14, color_by="phenoName", name="Experiment", component=1, point_size=1.5, label_offset=4)
plotStateTrajectory(E14, color_by="featureName", name="Lum", component=1, point_size=5)

# Make trajectory
E14 = selectTrajectory(E14, component=1)
E14 = fitTrajectory(E14)
plotTrajectoryFit(E14) 

write.ygraphml(sce=E14,file='~/Dropbox/Tonotopy_project/Data/Tonotopy_analysis_rerun/E14_trajectory_11states_fc1_5.graphml',color_by='phenoName',name='state',node_label='state')

tl = read.ygraphml("~/Dropbox/Tonotopy_project/Data/Tonotopy_analysis_rerun/E14_trajectory_11states_fc1_5.graphml")
plot(tl[,1:2], axes=FALSE, xlab="", ylab="", pch=20, cex=.25)
trajLayout(E14, adjust=TRUE) = tl

# saveRDS(E14,"~/Dropbox/Tonotopy_project/Data/Tonotopy_analysis_rerun/E14_WT_celltrails_11_subclusters_fc1_5_rerun.RDS")


