#!/usr/bin/env Rscript
#Run with R CMD BATCH SeuratSCT.R instead of Rscript to redirect std out to the log file
args = commandArgs(trailingOnly=TRUE)

# Test if there is one argument: if not, return an error
if (length(args) != 10) {
  stop("Other args must be supplied { path/to/info.csv, path/to/results, id, n. min genes per-cell, max genes per-cell, % MT, n. of PCs, cluster resolution (0.4 - 1.2), species (Human|Mouse), tissue (Blood|Brain) }", call.=FALSE)
} else if (length(args) == 10) {
  input = args[1]
  path = args[2]
  id = args[3]
  minGenes = as.numeric(args[4])
  maxGenes = as.numeric(args[5])
  percMT = as.numeric(args[6])
  npc = as.numeric(args[7])
  res = as.numeric(args[8])
  species = args[9]
  tissue = args[10]
}

# Load the libraries
suppressPackageStartupMessages({
    library(dplyr)
    library(Seurat)
    #library(celldex)
    #library(SingleCellExperiment)
    #library(SingleR)
    library(scCATCH)
    library(patchwork)
    library(ggplot2)
    library(RColorBrewer)
    library(cowplot)
    library(data.table)
})


SeuratObj <- function(sample, condition, path) {
  sobj <- Read10X(data.dir=path)
  sobj <- CreateSeuratObject(counts=sobj, project=sample, min.cells=3, min.features=200)
  sobj <- AddMetaData(sobj, metadata=c(sample, condition), col.name=c("Sample","Condition"))
  return(sobj)
}


# Output directories
out <- file.path(path,id)

if (dir.exists(out)) {
  write("Output directory already exists!", stdout())
} else {
  dir.create(out, recursive=TRUE)
}


write("Starting Seurat 10X Object...", stdout())
info <- fread(input, header=F, col.names=c("sample","condition","path"))
v <- mapply(SeuratObj, info$sample, info$condition, info$path)
data <- merge(v[[1]], v[2:length(v)], add.cell.ids=info$sample, project=id)
data <- AddMetaData(data, metadata=id, col.name="Project")
v <- NULL
#head(data@meta.data)
data
table(data$Project)
table(data$Sample)
table(data$Condition)


write("\nCompute QC metrics...", stdout())
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")

# We can define ribosomal proteins (their names begin with RPS or RPL), which often take substantial fraction of reads:
data[["percent.rb"]] <- PercentageFeatureSet(data, pattern="^RP[SL]")

# Visualize QC metrics as a violin plot
vp <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), group.by = "orig.ident", ncol=2, pt.size=0.1) & theme(plot.title = element_text(size=16))

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
fs1 <- FeatureScatter(data, feature1="nCount_RNA", feature2="percent.mt", group.by = "orig.ident")
fs2 <- FeatureScatter(data, feature1="nCount_RNA", feature2="nFeature_RNA", group.by = "orig.ident")
fs3 <- FeatureScatter(data, feature1="nCount_RNA", feature2="percent.rb", group.by = "orig.ident")
fs4 <- FeatureScatter(data, feature1="percent.rb", feature2="percent.mt", group.by = "orig.ident")
fsp <- plot_grid(fs1, fs2, fs3, fs4, labels=c("A", "B", "C", "D"), ncol=2, align=c("h","v"), label_size=20)

# Compute the relative expression of each gene per cell
# Too slow and not that useful
#C <- data@assays$RNA@counts
#C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
#most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
#bme <- boxplot(t(as.matrix(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell", col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

# Filetering step
data <- subset(data, subset = nFeature_RNA > minGenes & nFeature_RNA < maxGenes & percent.mt < percMT)
data

qc <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), group.by = "orig.ident", ncol=2, pt.size=0.1) & theme(plot.title = element_text(size=16))

pdf(file=file.path(out, "1_QC.pdf"), height=14.4, width=25.6)
vp
fsp
#print(bme)
qc
invisible(dev.off())

write("Integrating Seurat Objects based on 'Condition'...", stdout())
data.list <- SplitObject(data, split.by = "Condition")
data.list <- lapply(X = data.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = features)

# Variable features Venn
#hvgs_per_dataset <- lapply(data.list, function(x) {
#    x@assays$SCT@var.features
#})
#
#pdf(file=file.path(out, "2_Venn.pdf"), height=14.4, width=25.6)
#venn::venn(hvgs_per_dataset, opacity = 0.4, zcolor = (scales::hue_pal())(3), cexsn = 1, cexil = 1, lwd = 1, col = "white", frame = F, borders = NA)
#invisible(dev.off())

anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = features)
combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

combined.sct <- RunPCA(combined.sct, verbose = FALSE)
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:npc)

combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:npc)
combined.sct <- FindClusters(combined.sct, resolution = res)

head(combined.sct@meta.data)

# We are skipping the VariableFeatures plot because
# based on #2778 and @satijalab 's reply you can't run this on the integrated slot of an integrated dataset.....which makes sense because that would mean that
# 1) from a coding perspective the meta.features slot would need to be populated with values for things like sct.detection_rate sct.gmean etc etc, and 
# 2) conceptually, the integrated slot only exists to allow for optimal clustering. You should definitely not use it for anything else (such as gene expression, most variable genes etc).
# https://github.com/satijalab/seurat/issues/2172

pdf(file=file.path(out, "3_PCA.pdf"), height=14.4, width=25.6)
#VizDimLoadings(data, dims = 1:2, reduction="pca")
VizDimLoadings(combined.sct, dims = 1:15, reduction = "pca") & theme(axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold"))
DimPlot(combined.sct, reduction="pca")
DimHeatmap(combined.sct, dims=1, cells=500, balanced=TRUE)
DimHeatmap(combined.sct, dims = 1:15, cells=500, balanced=TRUE)
ElbowPlot(combined.sct, ndims = 30)
#JackStrawPlot(data, dims = 1:15)
invisible(dev.off())

pdf(file=file.path(out, "4_UMAP.pdf"), height=14.4, width=25.6)
DimPlot(combined.sct, label = TRUE, reduction="umap") + DimPlot(combined.sct, group.by="Condition", reduction="umap")
DimPlot(combined.sct, label = TRUE, split.by="Condition", ncol=4)
DimPlot(combined.sct, label = TRUE, split.by="Sample", ncol=4)
#DimPlot(data, reduction="tsne")
#DimPlot(data, reduction="pca")
# Plot cluster id over cell population:
#DimPlot(data,label.size=4, repel=T, label=T)
invisible(dev.off())

#table(combined.sct@meta.data$seurat_clusters)


write("\nFinding BioMarkers...", stdout())
DefaultAssay(combined.sct) <- "RNA"
data.markers <- FindAllMarkers(combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#data.markers <- FindAllMarkers(combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
#data.markers <- FindAllMarkers(combined.sct, logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1, min.diff.pct = 0, only.pos = TRUE, max.cells.per.ident = 20, return.thresh = 1, assay = "RNA")

#data.markers %>%
#    group_by(cluster) %>%
#    top_n(n = 2, wt = avg_log2FC)

#VlnPlot(data, features = c("IL7R", "LTB"))
#VlnPlot(data, features = c("IL7R", "LTB"), slot = "counts", log = TRUE)

data.markers %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = avg_log2FC) -> top1

data.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

DefaultAssay(combined.sct) <- "integrated"

pdf(file=file.path(out, "5_Biomarkers.pdf"), height=14.4, width=25.6)
FeaturePlot(combined.sct, features = c(head(top1$gene, 12)))
DoHeatmap(combined.sct, features = top10$gene) + NoLegend()
invisible(dev.off())


write("\nPerforming cell type annotation...", stdout())
#clu_ann <- scCATCH(object = clu_markers$clu_markers,
clu_ann <- scCATCH(object = data.markers,
    species = species,
    cancer = NULL,
    tissue = tissue)
#print(clu_markers)
#print(clu_ann)

convertSeurat <- function(seurat_object, scCATCH_anno) {
  tmp1 <- data.frame(cluster = levels(Idents(seurat_object)))
  tmp <- merge(tmp1, scCATCH_anno, by = 'cluster', all = T)
  tmp$cell_type[which(is.na(tmp$cell_type))] <- "Unclassified"
  
  new.cluster.ids <- tmp$cell_type
  names(new.cluster.ids) <- levels(seurat_object)
  seurat_object <- RenameIdents(seurat_object, new.cluster.ids)
  
  return(seurat_object)
}

# data is the Seurat object and clu_ann the scCATCH result
combined.sct <- convertSeurat(combined.sct, clu_ann)

pdf(file=file.path(out, "6_AnnoUMAP.pdf"), height=14.4, width=25.6)
DimPlot(combined.sct, label = T , repel = T, label.size = 5) #+ NoLegend()
invisible(dev.off())
