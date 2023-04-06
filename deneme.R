rm(list = ls())

object[[]]
str(object)
object@meta.data

Bismillah <- object

Bismillah@meta.data

Bismillah$condition = 1

sed <-  Bismillah@meta.data

x = dim(sed)[1]

for (i in 1:x) {if(sed[i, 1] =="22 years"){sed[i, 5] = "Teen"}
  
else if(sed[i, 1] == "29 years"){sed[i, 5] = "Teen"}
  else if (sed[i, 1] == "57 years"){sed[i, 5] = "Grumpy"}
  else if (sed[i, 1] == "63 years"){sed[i, 5] = "Grumpy"}}
  
object@meta.data = sed

Idents(Bismillah) = 'condition'
DimPlot(Bismillah, split.by = "condition")


Bismillah$condition = 1


ted <-   Bismillah@meta.data


x = dim(ted)[1]

for (i in 1:x) {if(ted[i, 1] =="22 years"){ted[i, 5] = "Young"}
  
  else if(ted[i, 1] == "29 years"){ted[i, 5] = "Young"}
  else if (ted[i, 1] == "57 years"){ted[i, 5] = "Old"}
  else if (ted[i, 1] == "63 years"){ted[i, 5] = "Old"}}

Bismillah@meta.data = ted
Idents(Bismillah) = "condition"
DimPlot(Bismillah, split.by = "orig.ident", label=T)


Bismillah <- ScaleData(Bismillah, verbose = TRUE)


Bismillah <- RunPCA(object = Bismillah, npcs = 15, verbose = FALSE)

ElbowPlot(Bismillah, ndims = 15)



Bismillah = RunTSNE(Bismillah, dims = 1:15)
DimPlot(Bismillah)

Bismillah@meta.data = ted
Idents(Bismillah) = "Again"
DimPlot(Bismillah, split.by = "orig.ident", label=T)


Bismillah$included = 1
ted = Bismillah@meta.data

x = dim(ted)[1]
for (i in 1:x)
{
  if(ted[i,2] > 25000)
  {
    ted[i,7] = "Excluded"
  }
  else if(ted[i,2] < 2000)
  {
    ted[i,7] = "Excluded"
  }
  else if(ted[i,3] < 200)
  {
    ted[i,7] = "Excluded"
  }
  else if(ted[i,3] > 4000)
  {
    ted[i,7] = "Excluded"
  }
  else if(ted[i,4] > 15)
  {
    ted[i,7] = "Excluded"
  }
  
  else
  {
    ted[i,7] = "Included"
  }
}
Bismillah@meta.data = ted
Idents(Bismillah) = "included"
DimPlot(Bismillah, cols = c("lightgreen","black"))

Bismillah@meta.data = ted
Idents(Bismillah) = "included"
DimPlot(Bismillah, cols = c("lightgreen","black"), split.by = "condition")




DefaultAssay(Bismillah) <- "integrated"


FindClusters(object=Bismillah)

Bismillah <- FindNeighbors(Bismillah, dims = 1:15)

Bismillah <- FindClusters(Bismillah, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
Bismillah <- FindClusters(Bismillah, resolution = c(0.2, 0.21, 0.22, 0.23, 0.24))

View(Bismillah@meta.data)

DimPlot(Bismillah, group.by = "integrated_snn_res.0.22", label = T)

DimPlot(Bismillah, group.by = "integrated_snn_res.0.22", reduction = "tsne", label = T)
# hucrelere CLuster isimleri yerine seurat annotation isimlerini atama


Idents(Bismillah) <- Bismillah@meta.data@cell_types
#RNA assayle conserved markerlari bulma

DefaultAssay(Bismillah) <- "RNA"
nk.markers <- FindConservedMarkers(Bismillah, ident.1 = 7, grouping.var = "stim", verbose = FALSE)
head(nk.markers)

getwd()
#Find all markers
DefaultAssay(Bismillah) <- "RNA"

Idents(Bismillah) = 'integrated_snn_res.0.22'
Allmarkers = FindAllMarkers(Bismillah, only.pos = F, logfc.threshold = 0.25, min.pct = 0.2)
write.csv(Allmarkers, "/Allmarkers.csv")


Allmarkers %>%
  group_by("integrated_snn_res.0.22"") %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(Bismillah, features = top10$gene) + NoLegend()


rm(ls())
#cell labeling

ref <- celldex::HumanPrimaryCellAtlasData()
results <- SingleR(test = as.SingleCellExperiment(Bismillah), ref = ref, labels = ref$label.main)
coldata(ref) (check label main)

Bismillah$single_labels <- results$labels

Bismillah[[]]
head(Bismillah)
DimPlot(Bismillah, reduction = "umap", group.by = "single_labels", label = TRUE)
Idents(Bismillah) = "condition"
DimPlot(Bismillah, reduction = "tsne", group.by = "single_labels", label = TRUE)
FeaturePlot(alldata, reduction = "tsne", features = c("TMPRSS2"))	

Bismillah <- RunUMAP(Bismillah, dims = 1:15)

clusters <- DimPlot(Bismillah, reduction = "umap", group.by = "integrated_snn_res.0.1", label = T)

cluster1 <- DimPlot(Bismillah, reduction = "umap", group.by = "single_labels", label = TRUE) + theme(legend.position = "none")
cluster2 <- DimPlot(Bismillah, reduction = "umap", group.by = "integrated_snn_res.0.1", label = TRUE) + theme(legend.position = "none")

cluster1 + cluster2

getwd()

cluster3 <- DimPlot(Bismillah, reduction = "tsne", group.by = "condition", label = TRUE) #+ theme(legend.position = "none")
cluster4 <- DimPlot(Bismillah, reduction = "tsne", group.by = "integrated_snn_res.0.2", label = TRUE) #+ theme(legend.position = "none")
cluster5 <- FeaturePlot(Bismillah, features = c("SFTPC"), reduction = "tsne")

cluster3 + cluster4 + cluster5

cluster3 <- DimPlot(Bismillah, reduction = "umap", group.by = "condition") #+ theme(legend.position = "none")
cluster4 <- DimPlot(Bismillah, reduction = "umap", group.by = "integrated_snn_res.0.1", label = TRUE) #+ theme(legend.position = "none")
cluster5 <- FeaturePlot(Bismillah, features = c("SFTPC"))

cluster3 + cluster4 + cluster5

number_perCluster<- table(Bismillah@meta.data$condition, 
                          Bismillah@meta.data$integrated_snn_res.0.1,
                          Bismillah@meta.data$single_labels)
write.csv(number_perCluster, "C:/Users/sedat/Downloads/Bismillah-to-Single_Cell/table1.csv")

Idents(Bismillah) = "integrated_snn_res.0.1"
FeaturePlot(Bismillah, features = c("SFTPC", "SFTPB", "NAPSA", "SFTPA1", "SLPI"), reduction = "tsne"), split.by = "condition")

cluster0.markers <- FindMarkers(Bismillah, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster0.markers, n = 5)

#Find cluster
install.packages('BiocManager')
BiocManager::install('multtest')
install.packages('metap')
-------------------------
  library(metap)
install.packages('metap')
#Cluster O nun diger clusterlardan farki
DefaultAssay(Bismillah) <- "RNA"
markers_cluster_0 <- FindConservedMarkers(Bismillah, ident.1 = 0, grouping.var = "condition")
head(markers_cluster_0)

#visualize data
FeaturePlot(Bismillah, features = c("ALPL"), min.cutoff = "q10", reduction = "tsne")
FeaturePlot(Bismillah, features = c("C1QB"), min.cutoff = "q10", reduction = "tsne")
#rename cluster 3 idents
Idents(Bismillah)
Bismillah <- RenameIdents(Bismillah, "0" = "Epithelial Cells")
DimPlot(Bismillah, reduction = "tsne", label = T)

#kolayca hucreleri konditiona gore isimlendirme

Bismillah$celltype.cnd <- paste0(Bismillah$seurat_clusters, "_", Bismillah$condition)

Idents(Bismillah) <- Bismillah$celltype.cnd
DimPlot(Bismillah, reduction = "tsne", label = T)


top3 <- Bismillah_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
DoHeatmap(object, features = top3$gene)


DotPlot(Bismillah, features = top3$gene)+ theme(axis.text.x = element_text(angle = 45, hjust=1))-YAPAMADIM
bakbak


Bismillah.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


Bismillah$Celltype_condition <- paste0(Bismillah$single_labels, "_", Bismillah$condition)
Idents(Bismillah) <-  Bismillah$Celltype_condition


Bismillah$Age_Status <- paste0(Bismillah$orig.ident, "_", Bismillah$condition)

Idents(Bismillah) <- Bismillah$Age_Status

view(Bismillah@meta.data)
difference <- FindMarkers(Bismillah, ident.1 = "0_Old", ident.2 = "0_Young")

head(difference)
FeaturePlot(Bismillah, features = c("ALPL", "RAP1GAP", "HSPA1A", "RPS2Y1"), min.cutoff = "q10", split.by = "condition")

Bismillah[[]]
str(Bismillah)
View(Bismillah@meta.data)

Bismillah[["RNA"]]@data
