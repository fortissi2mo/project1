# load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(celldex)
library(SingleR)

### Set working directory ###
setwd("D:/YUHS/65.이영인교수님/2.ANALYSIS/QC")

### Load data ###
# modify path!
data <- Read10X(data.dir = '../../1.DATA/Study1-1.GSE181297_RAW/Ke01/')
# Should modify prefix!
prefix = "Ke02"
obj <- CreateSeuratObject(counts = data, project = prefix, min.cells = 3, min.features = 200)

### Cell QC ###
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
## YH cutoff
obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 5)

### Normalization ###
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
par(mfrow = c(1,2))
hist(Matrix::colSums(data), main="Total expression before normalization", xlab="Sum of expression")
hist(Matrix::colSums(obj), main="Total expression after normalization", xlab="Sum of expression")

### Feature Selection ###
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(obj), 10)
plot1 <- VariableFeaturePlot(obj)
plot2 <- LabelPoints(plot = plot1, points = top10)
plot1 + plot2

### Scaling Data ### 
## YH difference why not feature all genes?
obj <- ScaleData(obj, vars.to.regress = "percent.mt")

### Dimensiona Reduction (PCA) ###
obj <- RunPCA(obj, features = VariableFeatures(object = obj))
print(obj[["pca"]], dims = 1:5, nfeatures = 5)

# Visualization - DimHeatmap
DimHeatmap(obj, dims = 1:20, cells = 500, balanced = TRUE)

# Determine ‘significant’ PCs (1) - Jackstraw
## YH can't understand this analysis
obj <- JackStraw(obj, num.replicate = 100)
obj <- ScoreJackStraw(obj, dims = 1:20)
JackStrawPlot(obj, dims = 1:20)
# Determine ‘significant’ PCs (2) - Elbow
ElbowPlot(obj)

### Clustering ###
obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.5)
head(Idents(obj), 5)

### Visualization ###
# (1) t-SNE
obj <- RunTSNE(obj, dims = 1:20)
p1 <- DimPlot(obj, label = TRUE, label.size = 5, reduction = "tsne")
# (2) UMAP
obj <- RunUMAP(obj, dims = 1:20)
p2 <- DimPlot(obj, label = TRUE, label.size = 5, reduction = "umap")
# draw plots
p1+p2

### Save and Read RDS ###
outfile <- paste0("./",prefix,".rds")
saveRDS(obj, file = outfile)
obj <- readRDS(outfile)

### Detecting Doublets

library(DoubletFinder)

# Can run parameter pK optimization with paramSweep.
sweep.res <- paramSweep_v3(obj)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
bcmvn
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

# define the expected number of doublet cells

# Doublet rate as indicated from the Chromium user guide
# Multiple Rate (%)	# of Cells Loaded	# of Cells Recovered
# ~0.4%				~870				~500
# ~0.8%				~1700				~1000
# ~1.6%				~3500				~2000
# ~2.3%				~5300				~3000
# ~3.1%				~7000				~4000
# ~3.9%				~8700				~5000
# ~4.6%				~10500				~6000
# ~5.4%				~12200				~7000
# ~6.1%				~14000				~8000
# ~6.9%				~15700				~9000
# ~7.6%				~17400				~10000

nExp <- round(ncol(obj) * 0.08)  # expect 8% doublets
obj <- doubletFinder_v3(obj, pN = 0.25, pK = 0.28, nExp = nExp, PCs = 1:20)
table(obj$DF.classifications_0.25_0.28_824)

# Remove Doublets
obj <- subset(obj, subset = DF.classifications_0.25_0.28_824 == "Singlet")

### Feature Selection ###
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(obj), 10)
plot1 <- VariableFeaturePlot(obj)
plot2 <- LabelPoints(plot = plot1, points = top10)
plot1 + plot2

### Scaling Data ###
obj <- ScaleData(obj, vars.to.regress = "percent.mt")

### Dimensiona Reduction (PCA) ###
obj <- RunPCA(obj, features = VariableFeatures(object = obj))
print(obj[["pca"]], dims = 1:5, nfeatures = 5)

# Visualization - DimHeatmap
DimHeatmap(obj, dims = 1:20, cells = 500, balanced = TRUE)

# Determine ‘significant’ PCs (1) - Jackstraw
obj <- JackStraw(obj, num.replicate = 100)
obj <- ScoreJackStraw(obj, dims = 1:20)
JackStrawPlot(obj, dims = 1:20)
# Determine ‘significant’ PCs (2) - Elbow
ElbowPlot(obj)

### Clustering ###
obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.5)
head(Idents(obj), 5)

### Visualization ###
# (1) t-SNE
obj <- RunTSNE(obj, dims = 1:20)
p1 <- DimPlot(obj, label = TRUE, label.size = 5, reduction = "tsne")
# (2) UMAP
obj <- RunUMAP(obj, dims = 1:20)
p2 <- DimPlot(obj, label = TRUE, label.size = 5, reduction = "umap")
# draw plots
p1+p2

### Cell type Annotation (by SingleR)
hpca.se <- celldex::HumanPrimaryCellAtlasData()
hpca.se
head(hpca.se$label.main)
head(hpca.se$label.fine)

pred <- SingleR(test = GetAssayData(obj), ref = hpca.se, labels = hpca.se$label.main)
head(pred)
table(pred$labels)

obj <- AddMetaData(obj, pred$labels, col.name = "celltype")
p1 <- DimPlot(object = obj, reduction = "tsne", label = TRUE, label.size=5)
p2 <- DimPlot(object = obj, reduction = "tsne", group.by = "celltype", label = TRUE, label.size=5)
p1+p2

p3 <- DimPlot(object = obj, reduction = "umap", label = TRUE, label.size=5)
p4 <- DimPlot(object = obj, reduction = "umap", group.by = "celltype", label = TRUE, label.size=5)
p3+p4

p5 <- FeaturePlot(obj, reduction = "umap", features = c("NRXN1"))		# NRXN1: neuron marker
p4+p5

plotScoreHeatmap(pred)

### Save and Read RDS ###
outfile <- paste0("./",prefix,".rds")
saveRDS(obj, file = outfile)
obj <- readRDS(outfile)

### Integration multiple obj ###
# Add object to object list
obj1 <- readRDS("Ke01.rds")
obj2 <- readRDS("NS02.rds")
obj1 <- RenameCells(obj1, add.cell.id = "Ke01")
obj2 <- RenameCells(obj2, add.cell.id = "NS02")
obj.list <- c(obj1,obj2)

features <- SelectIntegrationFeatures(object.list = obj.list)
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
#anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features, reduction = "rpca")
combined <- IntegrateData(anchorset = anchors)

#DefaultAssay(combined) <- "integrated"
DefaultAssay(combined) <- "RNA"
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined, vars.to.regress = "percent.mt")
combined <- RunPCA(combined, features = VariableFeatures(object = combined))

DimHeatmap(combined, dims = 1:30, cells = 500, balanced = TRUE)
combined <- JackStraw(combined, num.replicate = 100)
combined <- ScoreJackStraw(combined, dims = 1:30)
JackStrawPlot(combined, dims = 1:30)
ElbowPlot(combined)

combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

combined <- RunTSNE(combined, dims = 1:30)
p1 <- DimPlot(combined, label = TRUE, label.size = 5, reduction = "tsne", group.by="orig.ident")
combined <- RunUMAP(combined, dims = 1:30)
p2 <- DimPlot(combined, label = TRUE, label.size = 5, reduction = "umap", group.by="orig.ident")
p1+p2

p3 <- DimPlot(object = combined, reduction = "umap", group.by = "celltype", label = TRUE, label.size=5)
p4 <- FeaturePlot(combined, reduction = "umap", features = c("NRXN1"))
p3+p4

saveRDS(combined, file = "combined.rds")

### Subset Analysis ###
# 전체적인 keloid vs normal scar sample에서 차이를 보이는 gene이 아닌, neuron cell에 specific하게 차이를 보이는 gene을 찾도록 Subset Analysis pipeline 수정 필요
sub <- subset(combined, subset = seurat_clusters == "17")
Idents(sub) <- sub$orig.ident
markers <- FindAllMarkers(sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.table(as.matrix(markers), file = "markers.txt", sep = "\t")

top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DotPlot(combined, features = top20$gene, cols = c("blue", "red"), dot.scale = 8, split.by = "orig.ident") + RotatedAxis()


