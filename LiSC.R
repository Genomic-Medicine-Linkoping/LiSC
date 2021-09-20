#!/usr/bin/env Rscript

'LiSC is a Seurat (v4.0.4) wrapper and uses scCATCH (v2.1) for cell-type annotation.

Usage:
  LiSC.R (single | integrate) <id> <csv> <out> [--minGPC=<n>, --maxGPC=<n>, -npc=<n>, --res=<n>, --species=<c>, --tissue=<c>, --SCT]
  LiSC.R (-h | --help)
  LiSC.R --version

Options:
  -h --help         Show this screen.
  --version         Show version.
  --SCT             Use SCTtransform algorithm for normalization.
  --minGPC=<n>      Min number of genes per cell [default: 200].
  --maxGPC=<n>      Max number of genes per cell [default: 2500].
  --percMT=<n>      Max percentage of MT [default: 5].
  -n --npc=<n>      Max number of Principal Components [default: 30].
  -r --res=<n>      Cluster resolution [default: 0.5].
  -s --species=<c>  Species for scCATCH [default: Human].
  -t --tissue=<c>   Tissue for scCATCH [default: Blood].
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
#print(args)
source("fun.R")
source("algos.R")

# Input/Output checks
info <- test_input(args$csv)
out <- file.path(args$out, args$id)
test_output(out)
#print(info)

write("Starting Seurat 10X Object...", stdout())
# Make Seurat object for single or integrated data analysis
if (args$single == TRUE) {
  SeuratSingleAlgo(info, args, out)
} else if (args$integrate == TRUE) {
    SeuratIntegrateAlgo(info, args, out)   
}
