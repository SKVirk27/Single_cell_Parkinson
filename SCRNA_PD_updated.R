# script to integrate scRNA-Seq datasets to correct for batch effects  GSE184950
setwd("~/Desktop/Test")


# load libraries
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

# get data location
dirs <- list.dirs(path = 'Data project/', recursive = F, full.names = F)

for(x in dirs){
  name <- gsub('_filtered_feature_bc_matrix','', x)
  
  cts <- ReadMtx(mtx = paste0('Data project/',x,'/matrix.mtx.gz'),
                 features = paste0('Data project/',x,'/features.tsv.gz'),
                 cells = paste0('Data project/',x,'/barcodes.tsv.gz'))
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}




# merge datasets

ls()
merged_seurat <- merge(GSM5602315_A10_Control, y = c(GSM5602316_A15_Diseased,GSM5602324_B25_Diseased,GSM5602327_B5_Control,
                                                     GSM5602332_C20_Diseased,GSM5602336_D16_Diseased,GSM5602338_D27_Diseased,
                                                     GSM5602344_p4_Control,GSM5602347_p9_Diseased,GSM6042298_D14_Diseased),
                       add.cell.ids = ls()[3:12],
                       project = 'PKD')


merged_seurat

# QC & filtering -----------------------

View(merged_seurat@meta.data)
# create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Dataset','patient','Conditions', 'Barcode'), 
                                    sep = '_')

# calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')


# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
# QC
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')


# filtering
merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 800 &
                                   nFeature_RNA >500 &
                                   mitoPercent < 10)

merged_seurat_filtered


merged_seurat




# perform standard workflow steps to figure out if we see any batch effects --------
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)


# plot
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'patient')
p2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Conditions',
          cols = c('red','green','blue'))

grid.arrange(p1,p2, ncol = 2, nrow = 2)


# perform integration to correct for batch effects ------
obj.list <- SplitObject(merged_seurat_filtered, split.by = 'patient')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}


# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA)


anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors,k.weight = 50)
View(seurat.integrated@meta.data)

# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)


# View table

p3 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'patient')
p4 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Conditions',
              cols = c('red','green','blue'))
grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
p5 <- DimPlot(seurat.integrated, reduction ='umap', group.by ='seurat_clusters')
grid.arrange(p5 )
p6 <-DimPlot(seurat.integrated, label = TRUE, split.by = "Conditions")  + NoLegend()
grid.arrange(p6 )

# Cell annotation using SingleR library
ref <- celldex::HumanPrimaryCellAtlasData()
results <- SingleR(test=as.SingleCellExperiment(seurat.integrated),ref=ref,labels= ref$label.main)
results
seurat.integrated$singlr_labels <-results$labels
seurat.integrated[[]]
view(seurat.integrated)
# Visualize cell type in our data
DimPlot(seurat.integrated,reduction = 'umap',group.by = 'singlr_labels', label=TRUE)
# Save Seurat Integrated file
write.csv(seurat.integrated@meta.data,"./seurat_metadata_celltype.csv")
SaveH5Seurat(seurat.integrated, overwrite = TRUE)
saveRDS(seurat.integrated, file = "parkinson.rds")
# Finding top 10 variable features 
seurat.integrated<- FindVariableFeatures(seurat.integrated, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat.integrated), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat.integrated)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2
# Cell marker and cluster analysis
clusters <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
condition <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Conditions')



# let's visualize top features

FeaturePlot(seurat.integrated, features = c('GBA','LRRK2','PRKN','SNCA'),split.by = 'Conditions', min.cutoff = 'q10')
FeaturePlot(seurat.integrated, features = c('PARK7','PINK1','VPS35'),split.by = 'Conditions', min.cutoff = 'q10')
Idents(seurat.integrated) <- seurat.integrated@meta.data$singlr_labels
Idents(seurat.integrated)  
DimPlot(seurat.integrated, reduction = 'umap', label = TRUE)

filename <- file.choose()
seurat.integrated <- readRDS("parkinson.rds")
view(seurat.integrated)
seurat.integrated$Conditions <- sample(c("Control", "Diseased"), size = ncol(seurat.integrated), replace = TRUE)
features1 <- c('LRRK2','SNCA')

seurat.integrated
DotPlot(seurat.integrated, features = features1) + RotatedAxis()
# Find conserved markers
DefaultAssay(seurat_integrated) <- "RNA"
markers_cluster1 <- FindConservedMarkers(seurat.integrated,
                                         ident.1 = 1,
                                         grouping.var = 'Conditions')

head(markers_cluster1)
markers_cluster2 <- FindConservedMarkers(seurat.integrated,
                                         ident.1 = 2,
                                         grouping.var = 'Conditions')
head(markers_cluster2)
