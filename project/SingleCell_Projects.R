packageVersion("Matrix")
# Bioconductor packages
# Puis charger le package
library(celldex)
# Seurat from CRAN
#install.packages("Seurat")
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SingleR)
library(celldex)
library(RColorBrewer)
setwd("~/Single_analysis/project/data")
# set the seed  reproductibility
set.seed(1234)
# where i downloaded the file
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243665
#https://pubmed.ncbi.nlm.nih.gov/38307867/


####### 1- "Download and Prepare the Data Download raw 10X data from any dataset
#######   load the matrix and create an object for analysis."

raw_data<- Read10X(data.dir = "~/Single_analysis/project/data")
seurat_obj <- CreateSeuratObject(counts = raw_data[['Gene Expression']],
                                 project = "lung_cancer",
                                 min.cells = 3,# genes detected in >= 3 cells
                                 min.features = 200 # min 200 genes in one cells
                                  )

##### 2- Preprocess the Data
##########a- Perform quality control (filter low-quality cells).
###############calculate the percentage of  mitochondrial genes per cell
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
###############Visualize key QC Metrics per cell using violin plots
###############nFeature_RNA number of genes detected per cell
###############ncount_RNA Total count of RNA  per cell
###############percent.mt mithocondrial gene percentage per cell
VlnPlot(seurat_obj,features=c("nFeature_RNA","nCount_RNA","percent.mt"),
        pt.size = 0.1,
        ncol= 3)

####################filter low-quality cells###############
## cells expressing fewer than 500 genes (likely dead)
## cells expressing more than 6000 genes (potentially doublets)
## cells with > 15% mithocondrial gene contents (stressed or dead)
seurat_obj <- subset(seurat_obj,
                     subset=nFeature_RNA >500 & 
                       nCount_RNA < 6000 & 
                       percent.mt < 8)
############b- Normalize the data
############# using long-Normalisation
# Transform raw UMI counts into Normalized value to account for sequencing depth
seurat_obj <- NormalizeData(seurat_obj, 
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
#############c- Identify highly variable genes (2000)
seurat_obj <- FindVariableFeatures(seurat_obj,
                                   selection.method = "vst",
                                   nfeatures = 2000)
################ Visualize the variable features
# The top twenty variable genes are labeled for interpretability
var_plot <- VariableFeaturePlot(seurat_obj)
LabelPoints(plot = var_plot,
            points = head(VariableFeatures(seurat_obj),10),
            repel = TRUE)
#############d- Run PCA
############### before pca 
############## scale and center the expression values of each gene accross all cells
seurat_obj <- ScaleData(seurat_obj)
###################### RUN PCA
seurat_obj <- RunPCA(seurat_obj,npcs = 50)

########## 3- Dimensionality Reduction & Clustering
########## a- Use Elbow Plot or other criteria to select principal components.
#Visualize the standard deviation of each PC to decide how many pcs
ElbowPlot(seurat_obj,ndims = 50)

pcs <- 11 # this value depending on your datasets
###########b-identify the k-nearest neighbours for each cell based on the selected pcs
seurat_obj <- FindNeighbors(seurat_obj,dims = 1:pcs)
###########c- Perform clustering (e.g., Leiden)
install.packages('leidenbase')
seurat_obj <- FindClusters(seurat_obj, resolution = 0.6)
########### d - Visualize clusters using UMAP.
seurat_obj <- RunUMAP(seurat_obj, dims = 1:pcs)
##### Plot UMAP Clustering results
DimPlot(seurat_obj,reduction = "umap",
        label = TRUE , repel = TRUE) + ggtitle("UMAP: Cluster cells lung cancer")

########## 4-Annotate Cell Types
######### annotate clusters with potential cell types.
# load a built-in reference dataset(HUMAN RNA-seq) from celldex package
ref <- celldex::HumanPrimaryCellAtlasData()
##### using expression profile similarity
annotations <- SingleR(test = GetAssayData(seurat_obj, layer = "data")
                       # normalized data
                       ,ref = ref, 
                       labels = ref$label.main)
# Add predicted cell type label to the seurat object Metadata
seurat_obj$celltype <- annotations$labels

# Plot UMAP with predicted cell types lung cancer
DimPlot(seurat_obj, group.by = "celltype", label = TRUE, repel = TRUE) +
  ggtitle("UMAP: Annotated cell Types lung cancer Human")

######### 5- Marker Gene Identification
########### Identify top marker genes for each cluster
markers <- FindAllMarkers(seurat_obj,
                          only.pos = TRUE, #keep only positive Markers 
                          min.pct=0.25, # expressed in at least 25% of cells in a cluster
                          logfc.threshold = 0.25 # minimum log fold change
)
####### Identify top 20 Markers per cluster based on av_log2FC
library(dplyr)

top_markers <- markers %>% 
  group_by(cluster)%>%
  slice_max(n = 5 ,order_by = avg_log2FC)
###### Save Markers to cluster csv file
write.csv(top_markers,"cluster_markers.csv")
############## Visualize using violin plots, heatmaps, or dot plots
# DotPlot
DotPlot(seurat_obj, features = unique(top_markers$gene)) + RotatedAxis()

# Heatmap
DoHeatmap(seurat_obj, features = unique(top_markers$gene)) + NoLegend()

############# 6- Differential Expression Analysis
################Compare gene expression between two biologically
################relevant clusters T cells vs tumor cells).
####### Differential expression (FindMarkers)
de_genes <- FindMarkers(seurat_obj, ident.1 = 0, ident.2 = 1,
                        min.pct=0.25, logfc.threshold=0.25)
# view top de_genes
head(de_genes)
#save to csv
write.csv(de_genes, "cluster_0_vs_cluster_1_DEGs.csv")
###############Generate a volcano plot.
# Add column for significance
de_genes$gene <- rownames(de_genes)
de_genes$significant <- ifelse(de_genes$p_val_adj< 0.05 
                               & abs(de_genes$avg_log2FC)> 0.5,"Yes","No")
# Plot with ggplot
library(ggplot2)
ggplot(de_genes,aes(x=avg_log2FC, y= -log10(p_val_adj), color = significant)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()+
  labs(title = "Volcano Plots : Cluster_0 vs Cluster_1 ",
       x= "log2 Fold Change",
       y= "-Log10 adjusted P value")
