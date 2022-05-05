## E14.5 CellTrails subpopulation identification
library(CellTrails)
library(dplyr)
library(Seurat)
library(ggplot2)

# Load the data 
E14_1 = readRDS("~/Dropbox/Tonotopy_project/Data/Tonotopy_analysis/E14_WT_1_Cochlear_Duct_Roof_seurat.RDS")
E14_1 = GetAssayData(object = E14_1, slot = 'data')

E14_2 = readRDS("~/Dropbox/Tonotopy_project/Data/Tonotopy_analysis/E14_WT_2_Cochlear_Duct_Base_seurat.RDS")
E14_2 = GetAssayData(object = E14_2, slot = 'data')

E14_3 = readRDS("~/Dropbox/Tonotopy_project/Data/Tonotopy_analysis/E14_WT_3_Cochlear_Duct_Apex_seurat.RDS")
E14_3 = GetAssayData(object = E14_3, slot = 'data')

index = read.csv("~/Dropbox/Tonotopy_project/Data/data/gene_index_no_cell_cycle_with_1_starting.csv",header = FALSE)
index = index$V1
E14_1 = as.data.frame(E14_1)
l4 = rep("Roof",dim(E14_1)[2])
E14_2 = as.data.frame(E14_2)
l5 = rep("Base",dim(E14_2)[2])
E14_3 = as.data.frame(E14_3)
l6 = rep("Apex",dim(E14_3)[2])

# Generate SCE object to prepare for celltrails
E14 = as.matrix(cbind(E14_1,E14_2,E14_3))[index,]
E14 = SingleCellExperiment(assays=list(logcounts=E14))
isSpike(E14, "ERCC") = 1:80

cluster = c(l4,l5,l6)
cluster = as.factor(cluster)

# Removes features that are not expressed or that do not sufficiently reach the technological limit of detection
trajFeatureNames(E14) = filterTrajFeaturesByDL(E14, threshold=2, show_plot=FALSE)
trajFeatureNames(E14) = filterTrajFeaturesByCOV(E14, threshold=0.5, show_plot=FALSE)
E14_sub = E14[trajFeatureNames(E14), ]

# Filter using scran
var_fit = scran::trendVar(x=E14_sub, use.spikes=FALSE)
var_out = scran::decomposeVar(x=E14_sub, fit=var_fit)
tfeat = featureNames(E14_sub)[which(var_out$FDR < 0.01)]

trajFeatureNames(E14) = tfeat
showTrajInfo(E14)

# Dimension reduction 
se = embedSamples(E14)
d = findSpectrum(se$eigenvalues, frac=100)
latentSpace(E14) = se$components[, d]

# Find states
cl <- findStates(E14, min_size=0.005, min_feat=2, max_pval=1e-4, min_fc=1)
states(E14) = cl
plotManifold(E14, color_by="phenoName", name="state")

E14$Experiment = cluster
plotManifold(E14, color_by="phenoName", name="Experiment")

# Sample ordering
E14 = connectStates(E14, l=8)
plotStateTrajectory(E14, color_by="phenoName", name="Experiment", component=1, point_size=1.5, label_offset=4)
plotStateTrajectory(E14, color_by="featureName", name="Lum", component=1, point_size=5)

# Make trajectory
E14 = selectTrajectory(E14, component=1)
E14 = fitTrajectory(E14)
plotTrajectoryFit(E14) 

# write.ygraphml(sce=E14,file='~/Dropbox/Tonotopy_project/Data/Tonotopy_analysis/trajectory4.graphml',color_by='phenoName',name='state',node_label='state')

tl = read.ygraphml("~/Dropbox/Tonotopy_project/Data/Tonotopy_analysis/trajectory4.graphml")
plot(tl[,1:2], axes=FALSE, xlab="", ylab="", pch=20, cex=.25)


