# Libraries
library(Matrix)
library(dplyr)
library(ggplot2)
require(gridExtra)
library(RColorBrewer)
library(cowplot)
library(matrixStats)
library(Seurat)
library(matrixStats)
library(SingleR) #annotate pbmc
library(celldex) # reference gene set
library(SingleCellExperiment)
library(viridis)

my.theme <- theme(axis.title = element_text(size = 12), axis.text.x = element_text(angle = 0, hjust=0.5),
                  axis.text = element_text(size=12, color='black'), plot.title = element_text(size=12, face="plain"),
                  legend.position = "none")

# Import 10x data
data_dir <- '/Users/junjie/Desktop/OneDrive - The University of Chicago (1)/python/GPS/20250521/Ctrl_filtered_feature_bc_matrix'
list.files(data_dir)
mat <- Read10X(data.dir = data_dir)

#Create seurat_obj
seurat_obj <- CreateSeuratObject(counts = mat, project = "10x_pbmc", min.cells = 3, min.features = 200)

#data preprocessing
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, subset = nCount_RNA >= 1200 & nCount_RNA <= 30000)
seurat_obj <- subset(seurat_obj, subset = percent.mt < 20)

# Normalize and scale
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 4000)
seurat_obj <- ScaleData(seurat_obj, features=rownames(seurat_obj))

# PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs = 50)
DimPlot(seurat_obj, reduction = "pca")
ElbowPlot(seurat_obj, ndims=50)

# Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, k.param = 10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.2) # resolution = 0.2

seurat_obj <- RunTSNE(seurat_obj, dims = 1:30)
DimPlot(seurat_obj, reduction = "tsne", label=TRUE, label.size=8, pt.size=0.2) +
  xlab("mRNA t-SNE 1") + ylab("mRNA t-SNE 2") +
  theme(axis.title = element_text(size = 12), axis.text.x = element_text(angle = 0, hjust=0.5),
        axis.text = element_text(size=12, color='black')) + NoLegend()

# celldex PBMC data: https://bioconductor.org/packages/3.14/data/experiment/vignettes/celldex/inst/doc/userguide.html
ref <- celldex::NovershternHematopoieticData()

# Automated annotation
singleR.results <- SingleR(as.SingleCellExperiment(seurat_obj), assay.type.test = 1,
                           ref = ref, labels = ref$label.main)
table(singleR.results$labels)
plotScoreHeatmap(singleR.results)
seurat_obj[["singleR_ann_l1"]] <- singleR.results$labels

# Automated annotation
singleR.results <- SingleR(as.SingleCellExperiment(seurat_obj), assay.type.test = 1,
                           ref = ref, labels = ref$label.fine)
table(singleR.results$labels)
plotScoreHeatmap(singleR.results)
seurat_obj[["singleR_ann_l2"]] <- singleR.results$labels

#The annotation results can then be exported for python analysis
#Perform similar analysis for ConA sample

# Import 10x data
data_dir <- '/Users/junjie/Desktop/OneDrive - The University of Chicago (1)/python/GPS/20250521/ConA_filtered_feature_bc_matrix'
list.files(data_dir)
mat <- Read10X(data.dir = data_dir)

#Create seurat_obj
seurat_obj <- CreateSeuratObject(counts = mat, project = "10x_pbmc", min.cells = 3, min.features = 200)

#data preprocessing
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, subset = nCount_RNA >= 2000 & nCount_RNA <= 30000)
seurat_obj <- subset(seurat_obj, subset = percent.mt < 20)

# Normalize and scale
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 4000)
seurat_obj <- ScaleData(seurat_obj, features=rownames(seurat_obj))

# PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs = 50)
DimPlot(seurat_obj, reduction = "pca")
ElbowPlot(seurat_obj, ndims=50)

# Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, k.param = 10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.2) # resolution = 0.2

seurat_obj <- RunTSNE(seurat_obj, dims = 1:30)
DimPlot(seurat_obj, reduction = "tsne", label=TRUE, label.size=8, pt.size=0.2) +
  xlab("mRNA t-SNE 1") + ylab("mRNA t-SNE 2") +
  theme(axis.title = element_text(size = 12), axis.text.x = element_text(angle = 0, hjust=0.5),
        axis.text = element_text(size=12, color='black')) + NoLegend()

# celldex PBMC data: https://bioconductor.org/packages/3.14/data/experiment/vignettes/celldex/inst/doc/userguide.html
ref <- celldex::NovershternHematopoieticData()

# Automated annotation
singleR.results <- SingleR(as.SingleCellExperiment(seurat_obj), assay.type.test = 1,
                           ref = ref, labels = ref$label.main)
table(singleR.results$labels)
plotScoreHeatmap(singleR.results)
seurat_obj[["singleR_ann_l1"]] <- singleR.results$labels

# Automated annotation
singleR.results <- SingleR(as.SingleCellExperiment(seurat_obj), assay.type.test = 1,
                           ref = ref, labels = ref$label.fine)
table(singleR.results$labels)
plotScoreHeatmap(singleR.results)
seurat_obj[["singleR_ann_l2"]] <- singleR.results$labels