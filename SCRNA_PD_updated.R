Project: Integrating scRNA-Seq Datasets to Correct for Batch Effects
Overview
This project demonstrates the integration and batch effect correction of scRNA-Seq datasets using the Seurat package in R. The datasets used in this project are from the study GSE184950, and the objective is to integrate multiple scRNA-Seq datasets, correct for batch effects, and perform downstream analyses such as cell clustering and cell type annotation.

Libraries Required

library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(SingleR)
library(SeuratDisk)
library(SeuratData)
library(celldex)
library(SingleCellExperiment)
library(ExperimentHub)
library(scRNAseq)
library(patchwork)
library(performance)
library(metap)
library(multtest)
Steps Involved
1. Data Loading

# Set working directory
setwd("~/Desktop/Test")

# Get data location
dirs <- list.dirs(path = 'Data project/', recursive = F, full.names = F)

# Load and create Seurat objects
for(x in dirs){
  name <- gsub('_filtered_feature_bc_matrix','', x)
  
  cts <- ReadMtx(mtx = paste0('Data project/',x,'/matrix.mtx.gz'),
                 features = paste0('Data project/',x,'/features.tsv.gz'),
                 cells = paste0('Data project/',x,'/barcodes.tsv.gz'))
  
  assign(name, CreateSeuratObject(counts = cts))
}
2. Data Merging

# Merge datasets
merged_seurat <- merge(GSM5602315_A10_Control, y = c(GSM5602316_A15_Diseased, GSM5602324_B25_Diseased, GSM5602327_B5_Control,
                                                     GSM5602332_C20_Diseased, GSM5602336_D16_Diseased, GSM5602338_D27_Diseased,
                                                     GSM5602344_p4_Control, GSM5602347_p9_Diseased, GSM6042298_D14_Diseased),
                       add.cell.ids = ls()[3:12],
                       project = 'PKD')
3. Quality Control (QC) and Filtering
r
Copy code
# Create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# Split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Dataset','patient','Conditions', 'Barcode'), sep = '_')

# Calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')

# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# QC plots
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')

# Filtering
merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 800 & nFeature_RNA > 500 & mitoPercent < 10)
4. Data Normalization and Integration

# Standard workflow steps to check for batch effects
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)

# Integration to correct for batch effects
obj.list <- SplitObject(merged_seurat_filtered, split.by = 'patient')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}

# Select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# Find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

# Integrate data
seurat.integrated <- IntegrateData(anchorset = anchors, k.weight = 50)
5. Downstream Analysis and Visualization

# Scale data, run PCA and UMAP, visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)

# Visualization
p1 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'patient')
p2 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Conditions', cols = c('red','green','blue'))
grid.arrange(p1, p2, ncol = 2, nrow = 2)

# Cell annotation using SingleR
ref <- celldex::HumanPrimaryCellAtlasData()
results <- SingleR(test = as.SingleCellExperiment(seurat.integrated), ref = ref, labels = ref$label.main)
seurat.integrated$singlr_labels <- results$labels
DimPlot(seurat.integrated, reduction = 'umap', group.by = 'singlr_labels', label = TRUE)

# Save integrated data
write.csv(seurat.integrated@meta.data, "./seurat_metadata_celltype.csv")
SaveH5Seurat(seurat.integrated, overwrite = TRUE)
saveRDS(seurat.integrated, file = "parkinson.rds")

# Finding top 10 variable features
seurat.integrated <- FindVariableFeatures(seurat.integrated, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seurat.integrated), 10)

# Plot variable features
plot1 <- VariableFeaturePlot(seurat.integrated)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2
Key Highlights
Data Integration: Merged multiple scRNA-Seq datasets and corrected for batch effects.
Quality Control: Filtered and normalized data to ensure high-quality results.
Visualization: Created UMAP plots to visualize integrated data and cell types.
Cell Annotation: Used SingleR for cell type annotation.
Files
seurat_metadata_celltype.csv: Metadata with cell type annotations.
parkinson.rds: Integrated Seurat object saved as an RDS file.

Usage

To replicate the analysis, follow these steps:

Data Loading: Load the scRNA-Seq datasets into Seurat objects.
Data Merging: Merge the datasets into a single Seurat object.
Quality Control: Perform QC and filtering to retain high-quality cells.
Normalization and Integration: Normalize the data and integrate to correct for batch effects.
Downstream Analysis: Perform PCA, UMAP, and visualize the integrated data.
Cell Annotation: Annotate cell types using SingleR and save the results.






