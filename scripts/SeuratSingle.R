#!/usr/bin/env Rscript
#Run with R CMD BATCH Seurat_1sample.R instead of Rscript to redirect std out to the log file
args = commandArgs(trailingOnly=TRUE)

# Test if there is one argument: if not, return an error
if (length(args) != 10) {
  stop("Other arguments must be supplied { path/to/filtered_feature_bc_matrix, path/to/results, id, n. min genes per-cell, max genes per-cell, % MT, n. of PCs, cluster resolution (0.4 - 1.2), species (Human|Mouse), tissue (Blood|Brain) }", call.=FALSE)
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
})

# Output directories
out <- file.path(path,id)

if (dir.exists(out)) {
  write("Output directory already exists!", stdout())
} else {
  dir.create(out, recursive=TRUE)
}

write("Starting Seurat 10X Object...", stdout())
# Load the data dataset
cell.ranger <- Read10X(data.dir=input)

#adj.matrix[1:10, 1:5]
#str(adj.matrix)

#dense.size <- object.size(as.matrix(adj.matrix))
#dense.size

#sparse.size <- object.size(adj.matrix)
#sparse.size

# Initialize the Seurat object with the raw (non-normalized data).
data <- CreateSeuratObject(counts=cell.ranger, project=id, min.cells=3, min.features=200)
data
cell.ranger <- NULL
#str(data)

# Show QC metrics for the first 5 cells
#head(data@meta.data, 5)
# nCount_RNA -> number of UMI reads detected per cell
# nFeature_RNA -> number of expressed (detected) genes per same cell

#meta <- data@meta.data
#summary(meta$nCount_RNA)
#summary(meta$nFeature_RNA)


write("\nCompute QC metrics...", stdout())
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")

# we can define ribosomal proteins (their names begin with RPS or RPL), which often take substantial fraction of reads:
data[["percent.rb"]] <- PercentageFeatureSet(data, pattern="^RP[SL]")

# Visualize QC metrics as a violin plot
vp <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol=4, pt.size=0.1) & theme(plot.title = element_text(size=16))

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

fs1 <- FeatureScatter(data, feature1="nCount_RNA", feature2="percent.mt")
fs2 <- FeatureScatter(data, feature1="nCount_RNA", feature2="nFeature_RNA")
fs3 <- FeatureScatter(data, feature1="nCount_RNA", feature2="percent.rb")
fs4 <- FeatureScatter(data, feature1="percent.rb", feature2="percent.mt")

fsp <- plot_grid(fs1, fs2, fs3, fs4, labels=c("A", "B", "C", "D"), ncol=2, align=c("h","v"), label_size=20)

data <- subset(data, subset = nFeature_RNA > minGenes & nFeature_RNA < maxGenes & percent.mt < percMT)
data

qc <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), group.by = "orig.ident", ncol=2, pt.size=0.1) & theme(plot.title = element_text(size=16))

#pdf(file="./Seurat_results/QC.pdf", height=14.4, width=25.6)
pdf(file=file.path(out, "1_QC.pdf"), height=14.4, width=25.6)
vp
fsp
qc
invisible(dev.off())


write("Normalizing data...\n", stdout())
data <- NormalizeData(data, normalization.method="LogNormalize", scale.factor=10000)
#data[["RNA"]]@data

# Next step discovers the most variable features (genes) - these are usually most interesting for downstream analysis.
data <- FindVariableFeatures(data, selection.method="vst", nfeatures=2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data), 10)

# plot variable features with and without labels
vf1 <- VariableFeaturePlot(data)
vf2 <- LabelPoints(plot=vf1, points=top10, repel=TRUE, xnudge = 0, ynudge = 0)

pdf(file=file.path(out, "2_VariableFeatures.pdf"), height=14.4, width=25.6)
#vf1
suppressWarnings(print(vf2))
invisible(dev.off())

write("\n", stdout())
# Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA.
all.genes <- rownames(data)
data <- ScaleData(data, features=all.genes)


write("\nPerforming linear dimensional reduction...\n", stdout())
data <- RunPCA(data, features=VariableFeatures(object=data))

# Examine and visualize PCA results a few different ways
#print(data[["pca"]], dims = 1:5, nfeatures = 5)


# To overcome the extensive technical noise in any single feature for scRNA-seq data, Seurat clusters cells based on their PCA scores.
# The top principal components therefore represent a robust compression of the dataset. We randomly permute a subset of the data (1% by default) and rerun PCA
# We identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.

# The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). 
# Significant PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line).
#data <- JackStraw(data, num.replicate = 100)
#data <- ScoreJackStraw(data, dims = 1:20)

pdf(file=file.path(out, "3_PCA.pdf"), height=14.4, width=25.6)
#VizDimLoadings(data, dims = 1:2, reduction="pca")
VizDimLoadings(data, dims = 1:15, reduction = "pca") & theme(axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold"))
DimPlot(data, reduction="pca")
DimHeatmap(data, dims=1, cells=500, balanced=TRUE)
DimHeatmap(data, dims = 1:15, cells=500, balanced=TRUE)
ElbowPlot(data)
#JackStrawPlot(data, dims = 1:15)
invisible(dev.off())


write("\nClustering cells...", stdout())
data <- FindNeighbors(data, dims = 1:npc)
data <- FindClusters(data, resolution=res)
data <- RunUMAP(data, dims = 1:npc)
#data <- RunTSNE(data, dims = 1:npc)

pdf(file=file.path(out, "4_UMAP.pdf"), height=14.4, width=25.6)
DimPlot(data, reduction="umap")
#DimPlot(data, reduction="tsne")
#DimPlot(data, reduction="pca")
# Plot cluster id over cell population:
#DimPlot(data,label.size=4, repel=T, label=T)
invisible(dev.off())

#saveRDS(data, file = "../output/pbmc_tutorial.rds")


write("\nFinding differentially expressed features...", stdout())
data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
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

pdf(file=file.path(out, "5_Biomarkers.pdf"), height=14.4, width=25.6)
FeaturePlot(data, features = c(head(top1$gene, 12)))
DoHeatmap(data, features = top10$gene) + NoLegend()
invisible(dev.off())


write("\nPerforming cell type annotation...", stdout())
#monaco.ref <- celldex::MonacoImmuneData()

#hpca.ref <- celldex::HumanPrimaryCellAtlasData()
#dice.ref <- celldex::DatabaseImmuneCellExpressionData()

# Convert to Single Cell Experiment object
#sce <- as.SingleCellExperiment(DietSeurat(data))
#monaco.main <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.main, clusters = Idents(data))
#monaco.fine <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.fine, clusters = Idents(data))

# hpca.main <- SingleR(test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.main)
# hpca.fine <- SingleR(test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.fine)
# dice.main <- SingleR(test = sce, assay.type.test = 1, ref = dice.ref, labels = dice.ref$label.main)
# dice.fine <- SingleR(test = sce, assay.type.test = 1, ref = dice.ref, labels = dice.ref$label.fine)

#table(monaco.main$pruned.labels)

#table(hpca.main$pruned.labels)
#table(dice.main$pruned.labels)

#data@meta.data$monaco.main <- monaco.main$pruned.labels
#data@meta.data$monaco.fine <- monaco.fine$pruned.labels

# data@meta.data$hpca.main   <- hpca.main$pruned.labels
# data@meta.data$dice.main   <- dice.main$pruned.labels
# data@meta.data$hpca.fine   <- hpca.fine$pruned.labels
# data@meta.data$dice.fine   <- dice.fine$pruned.labels

#data <- SetIdent(data, value = "monaco.fine")
#pdf(file=file.path(out, "AnnoUMAP.pdf"), height=14.4, width=25.6)
#DimPlot(data, label = T , repel = T, label.size = 3) + NoLegend()
#invisible(dev.off())


#clu_markers <- findmarkergenes(object = data,
#    species = 'Human',
#    cluster = 'All',
#    match_CellMatch = FALSE,
#    cancer = NULL,
#    tissue = NULL,
#    cell_min_pct = 0.25,
#    logfc = 0.25,
#    pvalue = 0.05)


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
data <- convertSeurat(data, clu_ann)

pdf(file=file.path(out, "6_AnnoUMAP.pdf"), height=14.4, width=25.6)
DimPlot(data, label = T , repel = T, label.size = 5) #+ NoLegend()
invisible(dev.off())
