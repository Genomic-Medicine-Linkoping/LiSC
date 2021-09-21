#!/usr/bin/env Rscript

mySeurat <- combined
mySeurat$celltype.condition <- paste(Idents(mySeurat), mySeurat$Condition, sep="_")
head(mySeurat@meta.data)
mySeurat$celltype <- Idents(mySeurat)
head(mySeurat@meta.data)
Idents(mySeurat) <- "celltype.condition"
head(Idents(mySeurat))
FindMarkers(mySeurat, ident.1 = "Angiogenic T Cell_T1", ident.2="CD14+ Monocyte_T1", min.pct=0.25, logfc.threshold=0.25)
FindAllMarkers(mySeurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
allvsall <- FindAllMarkers(mySeurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

prova <- FindMarkers(mySeurat, ident.1 = "Angiogenic T Cell_T1", ident.2="CD14+ Monocyte_T1", min.pct=0.25, logfc.threshold=0.25)
min(prova$avg_log2FC)
p[p$avg_log2FC < 0,]
prova[prova$avg_log2FC < 0,]
?FindMarkers()
prova <- FindMarkers(mySeurat, ident.1 = "Angiogenic T Cell_T1", ident.2="CD14+ Monocyte_T1", min.pct=0.25, logfc.threshold=0.25, group.by="Condition")
mySeurat@meta.data
head(mySeurat@meta.data)
prova <- FindMarkers(mySeurat, ident.1 = "Angiogenic T Cell_T1", ident.2="CD14+ Monocyte_T1", min.pct=0.25, logfc.threshold=0.25, group.by="celltype.condition")

#savehistory("history.txt")

#tail(data@meta.data)
#table(data@meta.data$Project)
#table(data@meta.data$Sample)
#table(data@meta.data$Condition)
