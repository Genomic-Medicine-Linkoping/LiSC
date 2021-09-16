#!/usr/bin/env Rscript

test_input <- function(csv) {
    if (file.exists(csv)) {
    fread(csv, header=F, col.names=c("sample","condition","path"))
    } else {
    stop(sprintf("File %s does not exists!", csv), call.=FALSE)
    }
}


test_output <- function(out) {
    if (dir.exists(out)) {
    write("Output directory already exists!", stdout())
    } else {
    dir.create(out, recursive=TRUE)
    }
}


SeuratSingle <- function(info, id) {
    cell.ranger <- Read10X(data.dir=info$path)
    data <- CreateSeuratObject(counts=cell.ranger, project=id, min.cells=3, min.features=200)
    data <- AddMetaData(data, metadata=c(info$sample, info$condition, id), col.name=c("Sample","Condition","Project"))
    cell.ranger <- NULL
    return(data)
}

SeuratObj <- function(sample, condition, path) {
    sobj <- Read10X(data.dir=path)
    sobj <- CreateSeuratObject(counts=sobj, project=sample, min.cells=3, min.features=200)
    sobj <- AddMetaData(sobj, metadata=c(sample, condition), col.name=c("Sample","Condition"))
    return(sobj)
}

SeuratIntegrate <- function(SeuratObj, info, id) {
    v <- mapply(SeuratObj, info$sample, info$condition, info$path)
    data <- merge(v[[1]], v[2:length(v)], add.cell.ids=info$sample, project=id)
    data <- AddMetaData(data, metadata=id, col.name="Project")
    v <- NULL
    return(data)
}

SeuratQC <- function(data, minGenes, maxGenes, percMT, n, out) {
    write("\nCompute QC metrics...", stdout())
    data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")
    data[["percent.rb"]] <- PercentageFeatureSet(data, pattern="^RP[SL]")

    vp1 <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), group.by = "orig.ident", ncol=n, pt.size=0.1) & theme(plot.title = element_text(size=16))
    fs1 <- FeatureScatter(data, feature1="nCount_RNA", feature2="percent.mt", group.by = "orig.ident")
    fs2 <- FeatureScatter(data, feature1="nCount_RNA", feature2="nFeature_RNA", group.by = "orig.ident")
    fs3 <- FeatureScatter(data, feature1="nCount_RNA", feature2="percent.rb", group.by = "orig.ident")
    fs4 <- FeatureScatter(data, feature1="percent.rb", feature2="percent.mt", group.by = "orig.ident")
    fsp <- plot_grid(fs1, fs2, fs3, fs4, labels=c("A", "B", "C", "D"), ncol=2, align=c("h","v"), label_size=20)

    data <- subset(data, subset = nFeature_RNA > minGenes & nFeature_RNA < maxGenes & percent.mt < percMT)

    vp2 <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), group.by = "orig.ident", ncol=n, pt.size=0.1) & theme(plot.title = element_text(size=16))
    #vp3 <- plot_grid(vp1, vp2, labels=c("A", "B"), ncol=2, align=c("h","v"), label_size=20)
    #qcl <- c(vp1,fsp,vp2,vp3)
    pdf(file=file.path(out, "1_QC.pdf"), height=14.4, width=25.6)
    print(vp1)
    print(fsp)
    print(vp2)
    #print(vp3)
    invisible(dev.off())
    return(data)
}

SingleNorm <- function(data) {
    write("Normalizing data...\n", stdout())
    data <- NormalizeData(data, normalization.method="LogNormalize", scale.factor=10000)
    data <- FindVariableFeatures(data, selection.method="vst", nfeatures=2000)
    return(data)
}

IntegrateNorm <- function(data) {
    write("Normalizing data...\n", stdout())
    data.list <- SplitObject(data, split.by = "Condition")
    data <- lapply(X = data.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })
    return(data)
}

IntegrateList <- function(data.list) {
    features <- SelectIntegrationFeatures(object.list = data.list)
    anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features)
    combined <- IntegrateData(anchorset = anchors)
    return(combined)
}

SeuratPCA <- function(combined, npc, res) { 
# Run the standard workflow for visualization and clustering
    write("\nPerforming linear dimensional reduction and clustering...", stdout())
    combined <- ScaleData(combined, verbose = FALSE)
    combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
    combined <- FindNeighbors(combined, reduction = "pca", dims = 1:npc)
    combined <- FindClusters(combined, resolution = res)
    combined <- RunUMAP(combined, reduction = "pca", dims = 1:npc)
}

SeuratTop <- function(data.markers, n) { 
    data.markers %>%
    group_by(cluster) %>%
    top_n(n = n, wt = avg_log2FC)
}

SeuratAnno <- function(data.markers, species, tissue) { 
    write("\nPerforming cell type annotation...", stdout())
    scCATCH(object = data.markers,
        species = species,
        cancer = NULL,
        tissue = tissue)
}

convertSeurat <- function(seurat_object, scCATCH_anno) {
  tmp1 <- data.frame(cluster = levels(Idents(seurat_object)))
  tmp <- merge(tmp1, scCATCH_anno, by = 'cluster', all = T)
  tmp$cell_type[which(is.na(tmp$cell_type))] <- "Unclassified"
  
  new.cluster.ids <- tmp$cell_type
  names(new.cluster.ids) <- levels(seurat_object)
  seurat_object <- RenameIdents(seurat_object, new.cluster.ids)
  
  return(seurat_object)
}

SeuratPlot <- function(out, data, top1, top10) {
    #DefaultAssay(data) <- "RNA"
    #top10 <- head(VariableFeatures(data), 10)
    #vf1 <- VariableFeaturePlot(data)
    #vf2 <- LabelPoints(plot=vf1, points=top10, repel=TRUE, xnudge = 0, ynudge = 0)
    #pdf(file=file.path(out, "2_VariableFeatures.pdf"), height=14.4, width=25.6)
    #vf1
    #suppressWarnings(print(vf2))
    #invisible(dev.off())

    #DefaultAssay(combined) <- "data"
    pdf(file=file.path(out, "2_PCA.pdf"), height=14.4, width=25.6)
    #VizDimLoadings(data, dims = 1:2, reduction="pca")
    print(VizDimLoadings(data, dims = 1:15, reduction = "pca") & theme(axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold")))
    print(DimPlot(data, reduction="pca", group.by = "Condition"))
    print(DimHeatmap(data, dims=1, cells=500, balanced=TRUE))
    print(DimHeatmap(data, dims = 1:15, cells=500, balanced=TRUE))
    print(ElbowPlot(data))
    #JackStrawPlot(data, dims = 1:15)
    invisible(dev.off())

    pdf(file=file.path(out, "4_Biomarkers.pdf"), height=14.4, width=25.6)
    print(FeaturePlot(data, features = c(head(top1$gene, 12))))
    print(DoHeatmap(data, features = top10$gene) + NoLegend())
    invisible(dev.off())
}

PlotIntegratedUMAP <- function(out, fname, data, ncol) {
    pdf(file=file.path(out, fname), height=14.4, width=25.6)
    print(DimPlot(data, label = TRUE, reduction="umap") + DimPlot(data, group.by="Condition", reduction="umap"))
    print(DimPlot(data, label = TRUE, split.by="Condition", ncol=ncol))
    print(DimPlot(data, label = FALSE, split.by="Sample", ncol=ncol))
    invisible(dev.off())
}








SeuratPlotty <- function(out, list) {
    pdf(file=file.path(out, "1_QC.pdf"), height=14.4, width=25.6, onefile=TRUE)
    #do.call(list[[2]])
    plot.new()
    list
    invisible(dev.off())
}

SeuratPlotQC <- function(data, n, out) {
    vp1 <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), group.by = "orig.ident", ncol=n, pt.size=0.1) & theme(plot.title = element_text(size=16))
    fs1 <- FeatureScatter(data, feature1="nCount_RNA", feature2="percent.mt", group.by = "orig.ident")
    fs2 <- FeatureScatter(data, feature1="nCount_RNA", feature2="nFeature_RNA", group.by = "orig.ident")
    fs3 <- FeatureScatter(data, feature1="nCount_RNA", feature2="percent.rb", group.by = "orig.ident")
    fs4 <- FeatureScatter(data, feature1="percent.rb", feature2="percent.mt", group.by = "orig.ident")
    fsp <- plot_grid(fs1, fs2, fs3, fs4, labels=c("A", "B", "C", "D"), ncol=2, align=c("h","v"), label_size=20)
    vp2 <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), group.by = "orig.ident", ncol=n, pt.size=0.1) & theme(plot.title = element_text(size=16))
    vp3 <- plot_grid(vp1, vp2, labels=c("A", "B"), ncol=2, align=c("h","v"), label_size=20)
    qcl <- c(vp1,fsp,vp2,vp3)
    #pdf(file=file.path(out, "1_QC.pdf"), height=14.4, width=25.6)
    #print(vp1)
    #print(fsp)
    #print(vp2)
    #print(vp3)
    #invisible(dev.off())
    return(qcl)
}