library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(SingleCellExperiment)
library(SingleR)
setwd("C:/Users/sedat/Downloads/Bismillah-to-Single_Cell")


#clear environment
rm(list = ls())

#Young data donor 6:22 and donor 3:29

female_22 <- Read10X_h5("others/GSM3489193_Donor_06_filtered_gene_bc_matrices_h5.h5")
female_29<- Read10X_h5("others/GSM3489187_Donor_03_filtered_gene_bc_matrices_h5.h5")

#OLD data donor 1:63 and donor 4:57

female_57 <- Read10X_h5("others/GSM3489189_Donor_04_filtered_gene_bc_matrices_h5.h5")
female_63 <- Read10X_h5("others/GSM3489182_Donor_01_filtered_gene_bc_matrices_h5.h5")

#Seurat object asssignemnt

F_22y = CreateSeuratObject(female_22, "22 years")
F_29y = CreateSeuratObject(female_29, "29 years")
F_57y = CreateSeuratObject(female_57, "57 years")
F_63y = CreateSeuratObject(female_63, "63 years")

#QC analysis

F_22y[["percent.mt"]] <- PercentageFeatureSet(F_22y, pattern = "^MT-")
F_29y[["percent.mt"]] <- PercentageFeatureSet(F_29y, pattern = "^MT-")
F_57y[["percent.mt"]] <- PercentageFeatureSet(F_57y, pattern = "^MT-")
F_63y[["percent.mt"]] <- PercentageFeatureSet(F_63y, pattern = "^MT-")

#Filtering the data

F_22y <- subset(F_22y, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15 & nCount_RNA < 25000 & nCount_RNA > 2000)
F_29y <- subset(F_29y, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15 & nCount_RNA < 25000 & nCount_RNA > 2000)
F_57y <- subset(F_57y, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15 & nCount_RNA < 25000 & nCount_RNA > 2000)
F_63y <- subset(F_63y, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15 & nCount_RNA < 25000 & nCount_RNA > 2000)

#Normalization -burdayiz

F_22y <- NormalizeData(F_22y, normalization.method = "LogNormalize", scale.factor = 10000)
F_29y <- NormalizeData(F_29y, normalization.method = "LogNormalize", scale.factor = 10000)
F_57y <- NormalizeData(F_57y, normalization.method = "LogNormalize", scale.factor = 10000)
F_63y <- NormalizeData(F_63y, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features

F_22y <- FindVariableFeatures(F_22y, selection.method = "vst", nfeatures = 2000)
F_29y <- FindVariableFeatures(F_29y, selection.method = "vst", nfeatures = 2000)
F_57y <- FindVariableFeatures(F_57y, selection.method = "vst", nfeatures = 2000)
F_63y <- FindVariableFeatures(F_63y, selection.method = "vst", nfeatures = 2000)


#CCA integration



anchors <- FindIntegrationAnchors(object.list = list(F_22y,F_29y,F_57y,F_63y), dims = 1:20)
object <- IntegrateData(anchorset = anchors, dims = 1:20)

saveRDS(object, file='C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/object.rds')

#QC metrics in violin plot
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +   geom_hline(yintercept = c(low_cutoff, high_cutoff), linetype = "dashed", color = "red") + ggtitle(“Title”)
#COOOLLLLLL!!!

#vinplot without dot
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, point.size.use = NA)

VlnPlot(...) +
  geom_hline(yintercept = c(low_cutoff, high_cutoff), linetype = "dashed", color = "red") + ggtitle(“Title”)



#Principal component analysis (PCA)

DefaultAssay(object) <- "integrated"
object <- ScaleData(object, verbose = TRUE)
object <- RunPCA(object = object, npcs = 50, verbose = FALSE)
object <- FindVariableFeatures(object = object, selection.method = "vst", nfeatures = 2000)
object <- JackStraw(object = object, num.replicate = 100, dims = 50)
object <- ScoreJackStraw(object = object, dims = 1:50)
JackStrawPlot(object = object, dims = 1:50)



#Young ve Old a gore ayarlama


object$condition = 1
temp = object@meta.data

x = dim(temp)[1]
for (i in 1:x)
{
  if(temp[i,1] == "22 years")
  {
    temp[i,5] = "Young"
  }
  else if(temp[i,1] == "29 years")
  {
    temp[i,5] = "Young"
  }
  else if(temp[i,1] == "57 years")
  {
    temp[i,5] = "Old"
  }    
  else if(temp[i,1] == "63 years")
  {
    temp[i,5] = "Old"
  }
}
object@meta.data = temp
Idents(object) = "condition"
DimPlot(object, split.by = "condition")


#TSeNE de resolution ekleme ve conditiona gore calistirma

object <- FindNeighbors(object = object, dims = 1:20)
object = RunTSNE(object, dims = 1:20)
object <- FindClusters(object, resolution = 0.235)
DimPlot(object, reduction = "tsne", split.by = "condition")


#Hepsini okuturken object yap

object1 <- FindClusters(object, resolution = 0.15)
object2 <- FindClusters(object, resolution = 0.2)
object3 <- FindClusters(object, resolution = 0.25)
object4 <- FindClusters(object, resolution = 0.3)
object5 <- FindClusters(object, resolution = 0.35)

getwd()


object$celltype = paste(object$condition, object$seurat_clusters, sep = "_")
Idents(object) = "celltype"
DefaultAssay(object) = 'RNA'


markers_0 = FindMarkers(object, ident.1 = "Young_0", ident.2 = "Old_0", logfc.threshold = 0.1, min.pct = 0.2, only.pos = F)
write.csv(markers_0, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/markers_0.csv")
markers_1 = FindMarkers(object, ident.1 = "Young_1", ident.2 = "Old_1", logfc.threshold = 0.1, min.pct = 0.2, only.pos = F)
write.csv(markers_1, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/markers_1.csv")
markers_2 = FindMarkers(object, ident.1 = "Young_2", ident.2 = "Old_2", logfc.threshold = 0.1, min.pct = 0.2, only.pos = F)
write.csv(markers_2, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/markers_2.csv")
markers_3 = FindMarkers(object, ident.1 = "Young_3", ident.2 = "Old_3", logfc.threshold = 0.1, min.pct = 0.2, only.pos = F)
write.csv(markers_3, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/markers_3.csv")
markers_4 = FindMarkers(object, ident.1 = "Young_4", ident.2 = "Old_4", logfc.threshold = 0.1, min.pct = 0.2, only.pos = F)
write.csv(markers_4, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/markers_4.csv")
markers_5 = FindMarkers(object, ident.1 = "Young_5", ident.2 = "Old_5", logfc.threshold = 0.1, min.pct = 0.2, only.pos = F)
write.csv(markers_5, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/markers_5.csv")
markers_6 = FindMarkers(object, ident.1 = "Young_6", ident.2 = "Old_6", logfc.threshold = 0.1, min.pct = 0.2, only.pos = F)
write.csv(markers_6, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/markers_6.csv")
markers_7 = FindMarkers(object, ident.1 = "Young_7", ident.2 = "Old_7", logfc.threshold = 0.1, min.pct = 0.2, only.pos = F)
write.csv(markers_7, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/markers_7.csv")
markers_8 = FindMarkers(object, ident.1 = "Young_8", ident.2 = "Old_8", logfc.threshold = 0.1, min.pct = 0.2, only.pos = F)
write.csv(markers_8, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/markers_8.csv")
markers_9 = FindMarkers(object, ident.1 = "Young_9", ident.2 = "Old_9", logfc.threshold = 0.1, min.pct = 0.2, only.pos = F)
write.csv(markers_9, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/markers_9.csv")
markers_10 = FindMarkers(object, ident.1 = "Young_10", ident.2 = "Old_10", logfc.threshold = 0.1, min.pct = 0.2, only.pos = F)
write.csv(markers_10, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/markers_10.csv")
markers_11 = FindMarkers(object, ident.1 = "Young_11", ident.2 = "Old_11", logfc.threshold = 0.1, min.pct = 0.2, only.pos = F)
write.csv(markers_11, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/markers_11.csv")
markers_12 = FindMarkers(object, ident.1 = "Young_12", ident.2 = "Old_12", logfc.threshold = 0.1, min.pct = 0.2, only.pos = F)
write.csv(markers_12, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/markers_12.csv")
markers_13 = FindMarkers(object, ident.1 = "Young_13", ident.2 = "Old_13", logfc.threshold = 0.1, min.pct = 0.2, only.pos = F)
write.csv(markers_13, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/markers_13.csv")

Idents(object) = 'condition'
object_markers = FindAllMarkers(object, only.pos = T, logfc.threshold = 0.1, min.pct = 0.2)
write.csv(object_markers, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/bulk_markers.csv")

Idents(object) = 'seurat_clusters'
object_markers = FindAllMarkers(object, only.pos = F, logfc.threshold = 0.1, min.pct = 0.2)
write.csv(object_markers, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/Clusters_markers.csv")

object1 <- object

object@meta.data


#Sayfa genclerin yaslilardan farkli markerlari arastiriliyor

Idents(object1) = "orig.ident"
DefaultAssay(object1) = 'RNA'


markers_y22_o57 = FindMarkers(object1, ident.1 = "22 years", ident.2 = "57 years", logfc.threshold = 0.1, min.pct = 0.2, only.pos = F)
write.csv(markers_y22_o57, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/epi_markers_y22_o57.csv")

markers_y22_o63 = FindMarkers(object1, ident.1 = "22 years", ident.2 = "63 years", logfc.threshold = 0.1, min.pct = 0.2, only.pos = F)
write.csv(markers_y22_o63, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/epi_markers_y22_o63.csv")

markers_y29_o57 = FindMarkers(object1, ident.1 = "29 years", ident.2 = "57 years", logfc.threshold = 0.1, min.pct = 0.2, only.pos = F)
write.csv(markers_y29_o57, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/epi_markers_y29_o57.csv")

markers_y29_o63 = FindMarkers(object1, ident.1 = "29 years", ident.2 = "63 years", logfc.threshold = 0.1, min.pct = 0.2, only.pos = F)
write.csv(markers_y29_o63, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/epi_markers_y29_o63.csv")


data("object")
Endothelial_markers <- c("vWF")
DotPlot(object = object, features = Endothelial_markers)


library(SingleR)
singler=SingleR(method = "cluster", sc_data=object@assays$RNA@data, clusters = object@meta.data$seurat_clusters, genes = "de", quantile.use = 0.8, p.threshold = 0.05, fine.tune = TRUE, fine.tune.thres = 0.05, sd.thres = 1, do.pvals = T, numCores = 1,ref_data=hpca$data,types=hpca$types)


DotPlot(object, features = top3$gene)+ theme(axis.text.x = element_text(angle = 45, hjust=1))








top3 <- object_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

#Error var
DotPlot(object, features = top3$gene)+ theme(axis.text.x = element_text(angle = 45, hjust=1))




View(objects@meta.data)
rm(raw_objects)
