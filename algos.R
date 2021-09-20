#!/usr/bin/env Rscript

SeuratIntegrateAlgo <- function(SeuratObj, info, args, out) {
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
    markers <- SeuratAllMarkers(combined)
    PlotPCA(out, combined, "3_PCA.pdf")
    PlotIntegratedUMAP(out, combined, "4_UMAP.pdf", 4)
    PlotMarkers(out, combined, "5_Markers.pdf", markers)
    anno<- SeuratAnno(markers[[1]], args$species, args$tissue)
    combined <- convertSeurat(combined, anno)
    PlotIntegratedUMAP(out, combined, "6_annoUMAP.pdf", 4)
    saveRDS(combined, file.path(out, paste(args$id, "rds", sep=".")))
    print(combined)
}

#tail(data@meta.data)
#table(data@meta.data$Project)
#table(data@meta.data$Sample)
#table(data@meta.data$Condition)
