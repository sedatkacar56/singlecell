library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

library(tidyverse)
library(patchwork)
library(scCustomize)
library(qs)
library(ggpubr)

#plots with details https://samuel-marsh.github.io/scCustomize/articles/QC_Plots.html#vlnplot-based-qc-plots

p1 <- QC_Plots_Genes(raw_objects, low_cutoff = 200, high_cutoff = 5000)
p2 <- QC_Plots_UMIs(seurat_object = raw_objects, low_cutoff = 1200, high_cutoff = 45000)
p3 <- QC_Plots_Mito(raw_objects, mito_name = "percent_mito", high_cutoff = 10)
p4 <- QC_Plots_Complexity(seurat_object = raw_objects, high_cutoff = 0.8)
rm(raw_objects)
devtools::install_github(repo = "samuel-marsh/scCustomize")

remotes::install_github(repo = "samuel-marsh/scCustomize")

#Plotting mFeature, nCount RNA and percent MT
p1 <- VlnPlot(raw_objects, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #COOOLLLLLL!!!


p1 <- VlnPlot(raw_objects, features = "nFeature_RNA", ncol = 1) + theme(legend.position = "none") + theme(axis.title.x = element_blank()) +
  geom_hline(yintercept = c(low_cutoff = 200, high_cutoff = 5000), linetype = "dashed", color = "red") # + ggtitle("nFeature_RNA") + xlab("Frequency") + ylab("Frequency")


p2 <- VlnPlot(raw_objects, features = "nCount_RNA", ncol = 1) +  theme(legend.position = "none") + theme(axis.title.x = element_blank()) +
  geom_hline(yintercept = c(low_cutoff = 2000, high_cutoff = 25000), linetype = "dashed", color = "red") # + ggtitle("nCount_RNA") + xlab("Frequency") + ylab("Frequency")

p3 <- VlnPlot(raw_objects, features = "percent.mt", ncol = 1) +  theme(legend.position = "none") + theme(axis.title.x = element_blank()) +
  geom_hline(yintercept = c(high_cutoff = 15), linetype = "dashed", color = "red") # + ggtitle("percent.mt") + xlab("Frequency") + ylab("Frequency")

p4 <- ggarrange(p1, p2, p3, ncol=3) + theme(legend.position = "none")


raw_objects[[]]

#histogram of low cwll quality

F_22y[[]]

histogram_F22 <- hist(F_22y$nCount_RNA)
  abline(v=2000, col="red")
  abline(v=25000, col='red')
  
  histogram_F29 <- hist(F_29y$nCount_RNA)
  abline(v=2000, col="red")
  abline(v=25000, col='red')
  
  histogram_F57 <- hist(F_57y$nCount_RNA)
  abline(v=2000, col="red")
  abline(v=25000, col='red')
  
  histogram_F63 <- hist(F_63y$nCount_RNA)
  abline(v=2000, col="red")
  abline(v=25000, col='red')
  
  
  mt1 <- histogram_F22 <- hist(F_22y$percent.mt)
  abline(v=15, col="red")
  
  
  mt2 <- histogram_F29 <- hist(F_29y$percent.mt)
  abline(v=15, col="red")
  
  mt3 <- histogram_F57 <- hist(F_57y$percent.mt)
  abline(v=15, col="red")
  
  mt4 <-  histogram_F63 <- hist(F_63y$percent.mt)
  abline(v=15, col="red")
  
  
mt1 + mt2 + mt3 + mt4
  
rm(list = ls())
  
#Plotting mFeature, nCount RNA and percent MT without points!!!! YAN VE ALT YAZIYI SILME

p5 <- VlnPlot(raw_objects, pt.size = FALSE, features = "nFeature_RNA", ncol = 1) + theme(legend.position = "none") + theme(axis.title.x = element_blank()) +
  geom_hline(yintercept = c(low_cutoff = 200, high_cutoff = 5000), linetype = "dashed", color = "red") # + ggtitle("nFeature_RNA") + xlab("Frequency") + ylab("Frequency")


p6 <- VlnPlot(raw_objects, pt.size = FALSE, features = "nCount_RNA", ncol = 1) +  theme(legend.position = "none") + theme(axis.title.x = element_blank()) +
  geom_hline(yintercept = c(low_cutoff = 2000, high_cutoff = 25000), linetype = "dashed", color = "red") # + ggtitle("nCount_RNA") + xlab("Frequency") + ylab("Frequency")

p7 <- VlnPlot(raw_objects, pt.size = FALSE, features = "percent.mt", ncol = 1) +  theme(legend.position = "none") + theme(axis.title.x = element_blank()) +
  geom_hline(yintercept = c(high_cutoff = 15), linetype = "dashed", color = "red") # + ggtitle("percent.mt") + xlab("Frequency") + ylab("Frequency")

p8 <- ggarrange(p5, p6, p7, ncol=3)

help(lm)


plot1 <- FeatureScatter(raw_objects, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_hline(yintercept = c(high_cutoff = 15), color = "red") +  geom_vline (xintercept = c(low_cutoff = 2000, high_cutoff = 25000), color = "red")
plot2 <- FeatureScatter(raw_objects, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +   geom_hline(yintercept = c(low_cutoff = 200, high_cutoff = 5000), color = "red") +  geom_vline (xintercept = c(low_cutoff = 2000, high_cutoff = 25000), color = "red")
plot1 + plot2


#Elbow for PCA early

DefaultAssay(object) <- "integrated"
object <- ScaleData(object, verbose = TRUE)
object <- RunPCA(object = object, npcs = 50, verbose = FALSE)
object <- FindVariableFeatures(object = object, selection.method = "vst", nfeatures = 2000)
ElbowPlot(object, ndims = 40) 

#early tsne

object = RunTSNE(object, dims = 1:15)
DimPlot(object)

saveRDS(object, file='C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/trial.rds')

object_markers[[]]

#conditiona gore etiketleme

object$condition = 1
temp = object@meta.data

x = dim(temp)[1]
for (i in 1:x)
{
  if(temp[i,1] == "22 years")
  {
    temp[i,5] = "Young"
  }
  else if(temp[i,1] == "26 years")
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



table(object@meta.data$included, object@meta.data$condition)
table(object@meta.data$included, object@meta.data$orig.ident)

length(raw_objects@meta.data$condition)




object$included = 1
temp = object@meta.data

x = dim(temp)[1]
for (i in 1:x)
{
  if(temp[i,2] > 25000)
  {
    temp[i,6] = "Excluded"
  }
  else if(temp[i,2] < 2000)
  {
    temp[i,6] = "Excluded"
  }
  else if(temp[i,3] < 200)
  {
    temp[i,6] = "Excluded"
  }
  else if(temp[i,3] > 5000)
  {
    temp[i,6] = "Excluded"
  }
  else if(temp[i,4] > 15)
  {
    temp[i,6] = "Excluded"
  }
  else
  {
    temp[i,6] = "Included"
  }
}
object@meta.data = temp
Idents(object) = "included"
DimPlot(object, cols = c("lightgreen","black"))






