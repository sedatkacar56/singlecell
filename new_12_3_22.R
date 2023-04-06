
#script to integrate scRNA-Seq datasets to correct for batch effects
#setwd("~/Desktop/demo/single_cell_integrate/")

setwd('C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/')
library(dplyr)
library(Seurat)
library(patchwork)

library(SeuratDisk)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

library(hdf5r)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(SingleR)
browseVignettes("SingleR")
library(scater)
library(celldex)
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se
library(scRNAseq)
library(scRNAseq)
hESCs <- LaMannoBrainData('human-es')
hESCs <- hESCs[,1:100]
yes

save(hESCs, file = hESCs.RData)

library(SingleR)
pred.hesc <- SingleR(test = hESCs, ref = hpca.se, assay.type.test=1,
                     labels = hpca.se$label.main)
str(pred.hesc)

BiocManager::install("scRNAseq")
n
Read10X_h5("GSM3489182_Donor_01_filtered_gene_bc_matrices_h5.h5")
Read10X_h5("GSM3489183_IPF_01_filtered_gene_bc_matrices_h5.h5") 
Read10X_h5("GSM3489184_IPF_02_filtered_gene_bc_matrices_h5.h5")
Read10X_h5("GSM3489185_Donor_02_filtered_gene_bc_matrices_h5.h5")
Read10X_h5("GSM3489186_Cryobiopsy_01_filtered_gene_bc_matrices_h5.h5")
Read10X_h5("GSM3489187_Donor_03_filtered_gene_bc_matrices_h5.h5")
Read10X_h5("GSM3489188_IPF_03_filtered_gene_bc_matrices_h5.h5")
Read10X_h5("GSM3489189_Donor_04_filtered_gene_bc_matrices_h5.h5")
Read10X_h5("GSM3489190_IPF_04_filtered_gene_bc_matrices_h5.h5")
Read10X_h5("GSM3489191_Donor_05_filtered_gene_bc_matrices_h5.h5")
Read10X_h5("GSM3489192_HP_01_filtered_gene_bc_matrices_h5.h5")
Read10X_h5("GSM3489193_Donor_06_filtered_gene_bc_matrices_h5.h5")
Read10X_h5("GSM3489194_SSc-ILD_01_filtered_gene_bc_matrices_h5.h5")
Read10X_h5("GSM3489195_Donor_07_filtered_gene_bc_matrices_h5.h5")
Read10X_h5("GSM3489196_Myositis-ILD_01_filtered_gene_bc_matrices_h5.h5")
Read10X_h5("GSM3489197_Donor_08_filtered_gene_bc_matrices_h5.h5")
Read10X_h5("GSM3489198_SSc-ILD_02_filtered_gene_bc_matrices_h5.h5")

str(h5_seurat)
getwd()

h5_files <- list.files(pattern = ("*.h5"))

h5_read <- lapply(h5_files, Read10X_h5)
h5_seurat <- lapply(h5_read, CreateSeuratObject)
h5_seurat <- merge(h5_seurat[[1]], y = h5_seurat[2:length(h5_seurat)], add.cell.ids = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15","16", "17"), project = "project")  


# we used  https://satijalab.org/seurat/archive/v3.1/pbmc3k_tutorial.html




h5_seurat[["percent.mt"]] <- PercentageFeatureSet(h5_seurat, pattern = "^MT-")

VlnPlot(h5_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(h5_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(h5_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

h5_seurat <- subset(h5_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)



#normalizing the data

h5_seurat <- NormalizeData(h5_seurat)


#Identification of highly variable features

h5_seurat <- FindVariableFeatures(h5_seurat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(h5_seurat), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(h5_seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scaling data KALDIM

all.genes <- rownames(h5_seurat)
h5_seurat <- ScaleData(h5_seurat, features = all.genes)

#increase the memory
memory.limit(size=1000000000)
memory.size()

#garbage memory see all
gc()
gc(reset = TRUE)                 # Garbage collection
gc()
memory.limit(9999999999999)


#Renviron usage
install.packages("usethis")

library(usethis) 
usethis::edit_r_environ()

# Clear workspace
rm(list = ls())                  

#Perform linear dimensional reduction
h5_seurat <- RunPCA(h5_seurat, features = VariableFeatures(object = h5_seurat))

# Examine and visualize PCA results a few different ways
print(h5_seurat[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(h5_seurat, dims = 1:2, reduction = "pca")

DimPlot(h5_seurat, reduction = "pca")

DimHeatmap(h5_seurat, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(h5_seurat, dims = 1:20, cells = 500, balanced = TRUE)

#Determine the 'dimensionality' of the dataset
h5_seurat <- JackStraw(h5_seurat, num.replicate = 100)
h5_seurat <- ScoreJackStraw(h5_seurat, dims = 1:20)
JackStrawPlot(h5_seurat, dims = 1:20)

ElbowPlot(h5_seurat)

#CLuster the cells

h5_seurat <- FindNeighbors(h5_seurat, dims = 1:10)
h5_seurat <- FindClusters(h5_seurat, resolution = 0.1)

h5_seurat <- FindClusters(h5_seurat, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
View(h5_seurat@meta.data)


#understanding resolution

DimPlot(h5_seurat, group.by = "RNA_snn_res.0.1", label = TRUE)
DimPlot(h5_seurat, group.by = "RNA_snn_res.0.3", label = TRUE)

DimPlot(h5_seurat, group.by = "RNA_snn_res.0.5", label = TRUE)

DimPlot(h5_seurat, group.by = "RNA_snn_res.0.7", label = TRUE)

DimPlot(h5_seurat, group.by = "RNA_snn_res.1", label = TRUE)


# Look at cluster IDs of the first 5 cells
head(Idents(h5_seurat), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
h5_seurat <- RunUMAP(h5_seurat, dims = 1:0)
DimPlot(h5_seurat, reduction = "umap")


h5_seurat <- RunUMAP(h5_seurat, dims = 1:20)

DimPlot(h5_seurat, reduction = "umap")


saveRDS(h5_seurat, file ="kkk.rds")

load("kkk.rds")


#Finding differentially expressed features (cluster biomarkers)

# find all markers of cluster 1

cluster1.markers <- FindMarkers(h5_seurat, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)



# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(h5_seurat, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
 

# cluster 1 in markerlerini bul
cluster1.markers <- FindMarkers(h5_seurat, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster1.markers, n = 5)

#belli markerlerin violin plotla clusterlara gore expressionu 
VlnPlot(h5_seurat, features = c("ISG15", "KCNAB2"))

# you can plot raw counts as well
VlnPlot(h5_seurat, features = c("ISG15", "KCNAB2"), slot = "counts", log = TRUE)

#cluster umap mapi uzerinde expressyonlari
FeaturePlot(h5_seurat, features = c("ACE2", "TMPRSS2", "ISG15", "KCNAB2"))

rm(sed)

head(cluster1.markers, n = 5)
rm(raw_objects)

#heatmap
top10 <- h5_seurat.markers %>% group_by(cluster) %>% top_n(2)

top10 <- h5_seurat.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
DoHeatmap(h5_seurat, features = top10$gene) + NoLegend()

getwd()
top10 <- h5_seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(h5_seurat, features = top10$gene) + NoLegend()

#Assigning cell type identity to clusters
library(cellassign)
library(devtools)
tensorflow::tf_config()
install_tensorflow()
conda install -c "bioconda/label/cf201901" r-cellassign
tensorflow::tf_config()
reticulate::install_miniconda(force=TRUE, update=FALSE)
library(pheatmap)
library(stringr)
library(scater) # BioConductor
library(SingleCellExperiment) # BioConductor
library(DropletUtils) # BioConductor
library(DT)
install.packages("SingleCellExperiment")
install.packages("tensorflow", repos = "http://cran.us.r-project.org")
install.packages("cellassign")
remotes::install_github("rstudio/tensorflow")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(tensorflow)
install_tensorflow(version = "2.2.0-beta1")
library(tensorflow)
tf_config()
BiocManager::install("cellassign")
install.packages("reticulate")

reticulate::py_config()
```{r, eval=FALSE}
install.packages("stringr")
library(tensorflow)
tf_config()
install_tensorflow(extra_packages = "tensorflow-probability")
library(BiocManager)
library(reticulate)
path_to_python <- install_python()
virtualenv_create("r-reticulate", python = path_to_python)
library(tensorflow)
install_tensorflow(envname = "r-reticulate")
library(tensorflow)
tf$constant("Hello Tensorflow!")
BiocManager::install(c("GenomicRanges", "Organism.dplyr"))
install_tensorflow()
chooseCRANmirror()
install.packages("BiocManager")-
78
yes
BiocManager::valid()
library(devtools)
install.packages("tensorflow")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("devtools", repos = "http://cran.us.r-project.org") # If not already installed
devtools::install_github("Irrationone/cellassign")
conda install -c conda-forge -c bioconda r-cellassign
install.packages("devtools")   # unnecessary if you have it already
devtools::install_github("Bioconductor/BiocManager", ref="ghost-binary-repo")
a
install.packages("BiocManager")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("devtools") # If not already installed
devtools::install_github("Irrationone/cellassign")

BiocManager::install("BiocVersion")
R.version

.libPaths()                           # Get paths of installed packages

install.packages("installr")

library(installr)

updateR()


BiocManager::valid()
library(tidyverse) # CRAN
library(here) # CRAN
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

*1
 #cell markers are retrieved https://panglaodb.se/markers.html?cell_type=%27Airway%20goblet%20cells%27

marker_list <- list(
  `Macrophages, Alveolar Macrophages` = c("CD68", "MARCO", "FCGR3A", "LYZ", "PTPRC", "SIGLECF", "MRC1"),
  `ab T cells` = c("CD2", "CD3D", "TRAC", "IL32", "CD3E", "PTPRC"),
  `gd T cells` = c("NKG7", "FCGR3A", "HOPX", "GNLY", "KLRF1", "CMC1", "CCL3", "PTPRC"),
  `NK cells` = c("GZMK", "KLRF1", "CCL3", "CMC1", "NKG7", "PTPRC"),
  `Plasma cells` = c("CD27", "IGHG1", "CD79A", "IGHG2", "PTPRC", "IGKC"),
  `Mature B cells` = c("MS4A1", "LTB", "CD52", "IGHD", "CD79A", "PTPRC", "IGKC"),
  `Epithelial cells` = c("EPCAM", "CDH1"),
  `Endothelial cells` = c("PECAM1", "EMCN"),
  `Erythroid cells` = c("HBB", "SLC25A37", "CA1", "ALAS2"),
  `Pulmonary alveolar type I cells` = c("AGER", "CLDN18"),
  `Pulmonary alveolar type II cells` = c("SFTPC", "SFTPA1" ,"NKX2-1", "ETV5"),
  `Airway Goblet Cells` = c("MUC5B", "DUSP4","GGH")
)          
marker_list_to_mat(marker_list)





mgi <- marker_list_to_mat(Lung_marker_list, include_other = FALSE)

pheatmap(mgi)


str("LAST")




install.packages("cellassign")

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(h5_seurat)
h5_seurat <- RenameIdents(h5_seurat, new.cluster.ids)
DimPlot(h5_seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

load("kkk.rds")

harmony <- saveRDS(h5_seurat, file = "LAST.rds")
str(harmony)
