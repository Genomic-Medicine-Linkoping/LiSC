#!/usr/bin/env Rscript

'LiSC is a Seurat (v4.0.4) wrapper and uses scCATCH (v2.1) for cell-type annotation.

Usage:
  LiSC.R (single | integrate) <id> <csv> <out> [--minGPC=<n>, --maxGPC=<n>, -npc=<n>, --res=<n>, --species=<c>, --tissue=<c>]
  LiSC.R (-h | --help)
  LiSC.R --version

Options:
  -h --help         Show this screen.
  --version         Show version.
  --minGPC=<n>      min number of genes per cell [default: 200].
  --maxGPC=<n>      max number of genes per cell [default: 2500].
  --percMT=<n>      max percentage of MT [default: 5].
  -n --npc=<n>      max number of Principal Components [default: 30].
  -r --res=<n>      cluster resolution [default: 0.5].
  -s --species=<c>  species for scCATCH [default: Human].
  -t --tissue=<c>   tissue for scCATCH [default: Blood].
' -> doc

# Load the libraries
suppressPackageStartupMessages({
    library(docopt)
    library(dplyr)
    library(Seurat)
    library(data.table)
    library(scCATCH)
    library(patchwork)
    library(ggplot2)
    library(RColorBrewer)
    library(cowplot)
    library(gridExtra)
})

args <- docopt(doc, version = 'LiSC 0.1.0\n')
print(args)
source("fun.R")

# Input/Output checks
info <- test_input(args$csv)
out <- file.path(args$out, args$id)
test_output(out)
#print(info)

write("Starting Seurat 10X Object...", stdout())
# Make Seurat object for single or integrated data analysis
if (args$single == TRUE) {
    data <- SeuratSingle(info[1,], args$id)
    print(data)
    data <- SeuratQC(data, as.numeric(args$minGPC), as.numeric(args$maxGPC), as.numeric(args$percMT), 4, out)
    print(data)
    data <- SingleNorm(data)
    print(data)
} else if (args$integrate == TRUE) {
    data <- SeuratIntegrate(SeuratObj, info, args$id)
    print(data)
    data <- SeuratQC(data, as.numeric(args$minGPC), as.numeric(args$maxGPC), as.numeric(args$percMT), 2, out)
    print(data)
    data <- IntegrateNorm(data)
    combined <- IntegrateList(data)
    DefaultAssay(combined) <- "integrated"
    combined <- SeuratPCA(combined, as.numeric(args$npc), as.numeric(args$res))
    PlotIntegratedUMAP(out, "3_UMAP.pdf", combined, 4)    
    write("\nFinding differentially expressed features...", stdout())
    #DefaultAssay(combined) <- "RNA"
    data.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    top1 <- SeuratTop(data.markers, 1)
    top10 <- SeuratTop(data.markers, 10)
    SeuratPlot(out, combined, top1, top10)
    clu_ann <- SeuratAnno(data.markers, args$species, args$tissue)
    combined <- convertSeurat(combined, clu_ann)
    PlotIntegratedUMAP(out, "5_AnnoUMAP.pdf", combined, 4)    
}




#tail(data@meta.data)
#data
#table(data$Project)
#table(data$Sample)
#table(data$Condition)





