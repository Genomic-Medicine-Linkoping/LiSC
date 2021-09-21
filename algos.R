#!/usr/bin/env Rscript

SeuratSingleAlgo <- function(info, args, out) {
    data <- SeuratSingle(info[1,], args$id)
    print(data)
    data <- SeuratQC(data, as.numeric(args$minGPC), as.numeric(args$maxGPC), as.numeric(args$percMT), 4, out, "1_preQC.pdf")
    PlotQC(data, 4, out, "2_postQC.pdf")
    print(data)
    if (args$SCT == FALSE) {
        write("\nSCT is disabled!\n", stdout())
        data <- SingleNorm(data)
    } else {
        write("\nSCT is enabled!\n", stdout())
        data <- SCTransform(data, vars.to.regress = "percent.mt", verbose = TRUE)
    }
    data <- SeuratPCA(data, as.numeric(args$npc), as.numeric(args$res))
    PlotPCA(out, data, "3_PCA.pdf")
    PlotUMAP(out, data, "4_UMAP.pdf")
    markers <- SeuratAllMarkers(data)
    PlotMarkers(out, data, "5_Markers.pdf", markers)
    anno <- SeuratAnno(markers[[1]], args$species, args$tissue)
    data <- convertSeurat(data, anno)
    PlotUMAP(out, data, "6_annoUMAP.pdf")
    saveRDS(data, file.path(out, paste(args$id, "rds", sep=".")))
    print(data)
}

SeuratIntegrateAlgo <- function(info, args, out) {
    data <- SeuratIntegrate(SeuratObj, info, args$id)
    print(data)
    data <- SeuratQC(data, as.numeric(args$minGPC), as.numeric(args$maxGPC), as.numeric(args$percMT), 2, out, "1_preQC.pdf")
    PlotQC(data, 2, out, "2_postQC.pdf")
    print(data)
    if (args$SCT == FALSE) {
        write("\nSCT is disabled!\n", stdout())
        data.list <- SeuratNorm(data)
        combined <- IntegrateNorm(data.list)
        combined <- ScaleData(combined, verbose = FALSE)
    } else {
        write("\nSCT is enabled!\n", stdout())
        data.list <- SeuratSCT(data)
        combined <- IntegrateSCT(data.list)
    }
    combined <- SeuratPCA(combined, as.numeric(args$npc), as.numeric(args$res))
    PlotPCA(out, combined, "3_PCA.pdf")
    PlotIntegratedUMAP(out, combined, "4_UMAP.pdf", 4)
    if (args$SCT == FALSE) {
        DefaultAssay(combined) <- "RNA"
        all.genes <- rownames(combined)
        combined <- ScaleData(combined, features=all.genes, verbose = TRUE)
        markers <- SeuratAllMarkers(combined)
        PlotMarkers(out, combined, "5_Markers.pdf", markers)
    } else {
        DefaultAssay(combined) <- "SCT"
        markers <- SeuratAllMarkers(combined)
        PlotMarkers(out, combined, "5_Markers.pdf", markers)
    }    
    DefaultAssay(combined) <- "integrated"
    anno <- SeuratAnno(markers[[1]], args$species, args$tissue)
    combined <- convertSeurat(combined, anno)
    PlotIntegratedUMAP(out, combined, "6_annoUMAP.pdf", 4)
    saveRDS(combined, file.path(out, paste(args$id, "rds", sep=".")))
    print(combined)
}
