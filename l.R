

# set seed and put two plots in one figure
set.seed(123)
par(mfrow=c(1,2))
# original expression distribution
raw_geneExp = as.vector(Bismillah[['RNA']]@counts) %>% sample(10000)
raw_geneExp = raw_geneExp[raw_geneExp != 0]
hist(raw_geneExp)
# expression distribution after normalization
logNorm_geneExp = as.vector(Bismillah[['RNA']]@data) %>% sample(10000)
logNorm_geneExp = logNorm_geneExp[logNorm_geneExp != 0]
hist(logNorm_geneExp)


  Bismillah <- FindVariableFeatures(Bismillah, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Bismillah), 10)


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Bismillah) + 
  theme(legend.position="bottom")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + 
  theme(legend.position="none")
plot1 + plot2

trace()
savehistory('file')
myhistory <- scan('file','character')


VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
Bismillah@meta.data

set.seed(123)
par(mfrow=c(2,2))
DotPlot(Bismillah, features = c("NKG7", "MS4A1", "CD79A")) + theme(legend.position = "none" )

str(Bismillah)


DotPlot(Bismillah, features = c("NKG7"), slot = "counts", log = TRUE) + theme(legend.position = ) + theme(legend.position = "none")

Bismillah.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

Bismillah.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()