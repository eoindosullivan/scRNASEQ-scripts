## SCRIPT: classify a renal single cell dataset ->  basic walkthrough

## 19.11.20 MEDEOS @eoinrenal , https://www.youtube.com/channel/UCfsIprbFgwES4i0DG2kzDhA

## Data: Conway, O'Sullivan- "Single Cell Transcriptomics reveals renal myeloid dynamics post reversable unilateral ureteric injury"

## source:
## "Conway BR, O'Sullivan ED, Cairns C, O'Sullivan J et al.
## Kidney Single-Cell Atlas Reveals Myeloid Heterogeneity in Progression and Regression of Kidney Disease. 
## J Am Soc Nephrol 2020 Sep 25. PMID: 32978267")

## Heres a GEO example: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140023 ( my paper )

## make a new project folder always for each project!!

## the library is open ----

library(Seurat) ## your main single cell workflow toolkit
library(ggplot2) ## ubiquitous - better plotting functions than default, seurat plots are ggplot objects
library(tidyverse) ## ubiquitous - data manipulatation toolkit

library(reticulate) ## allows python functions to be called in R
library(Matrix) ## I need to read in a .mtx file

library(AnnotationDbi) ## to convert ensembl ID to symbol
library(org.Mm.eg.db)
library(viridis) ## beautiful and scientifically better heatmap colour schemes

## data ---- 

## start with public data in a matrix format initially which is easies to work with
## later you can get the raw data (.fastq files etc) and align yourself if youd like, which you need to learn to do anyway if its your own data
## get your data from GEO/ supplemental files in papers/  zenodo repository etc

## if files are are .tar.gz so unzip with 7zip twice to reveal a folder with .mtx, barcodes and genes.

## read in dataset ----

# usually convienient function ## change the names to "barcodes,features,matrix" if required, as in this case

dat.sham<-Read10X(
  data.dir = "GEO_example_Sham_renamed/",
  gene.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
) 

dat.sham
##change ENSEMBL ID to genesymbol----
## ensembl better for anlysis but less readable by humans
## ideally analyse as ENS and convert later, but...

## biomart takes ages and is cumbersome imo

#Mapping the Ensembl Gene ID back to the first symbol

x<- select(x = org.Mm.eg.db, 
       keys = rownames(dat.sham), 
       column = "SYMBOL", 
       keytype = "ENSEMBL",
       multiVals = "first",asNA=F)

head(x)
x$FINAL <- ifelse(is.na(x$SYMBOL), x$ENSEMBL, x$SYMBOL) # issues with excess matches so collapse into a single col

for (name in rownames(dat.sham)) {
  rownames(dat.sham)[match(name, rownames(dat.sham))]<-x$FINAL[match(name, x$ENSEMBL)]
}
head(rownames(dat.sham)) #check

##  create seurat object ----

## if this doesnt work convieniently for the data I want to look at though as the data its a little more than filenames to change, 
# as not in the optimal format required
## so easier to manually load in each 


dat <- CreateSeuratObject(
  dat.sham,
  project = "exampleset",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = NULL
)

## review the object

dat


## your best friend : Cheat sheet----
# https://satijalab.org/seurat/essential_commands.html

## QC

## now this is probably already done but nice to check anyway
## following the offical vignettes is the key from here
## https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html

dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^Mt") ## if human, "MT"

rownames(dat)
range(dat[["percent.mt"]])## !

VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dat <- subset(dat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # if Mt genes were present

## normalise ----

dat <- NormalizeData(dat)

## Variable gene selection ----
dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)

## what are our variable genes ?
top10 <- head(VariableFeatures(dat), 20)

# plot variable features with and without labels for the craic / sanity test
plot1 <- VariableFeaturePlot(dat)
plot1 <- LabelPoints(plot = plot1, points = top10, repel = T)
plot1

## scale----

all.genes <- rownames(dat)
dat <- ScaleData(dat, features = all.genes)

## dim reduction
dat <- RunPCA(dat, features = VariableFeatures(object = dat))
ElbowPlot(dat) # fastest, easiest - can dive in deeper with jackstraw, looking at the individiual loadings etc

## cluster ----

dat <- FindNeighbors(dat, dims = 1:20) ## iterativly decide the dims tbh
dat <- FindClusters(dat, resolution = 0.8) ## semi - supervise resolution, experiment! youll return through this cycle alot
dat <- RunUMAP(dat, dims = 1:20, n.neighbors = 20,min.dist = .01,spread = 6) ## twiddle the dials to get umaps of face validity

#dat <- RunUMAP(dat, dims = 1:15, n.neighbors = 4,min.dist = .001,spread = 5) ## twiddle the dials to get umaps of face validity

## plot 1 : dimension reduction plots ----

DimPlot(dat, reduction = "umap",pt.size = 2)

png(filename = "dat_umap_spread6.nn20_8cluster.png", width = 1000, height = 1000, pointsize = 12,family = "sans", bg = "white")
DimPlot(dat, reduction = "umap",pt.size = 6,)
dev.off()

## marker genes ----

dat.markers <- FindAllMarkers(dat, only.pos = F, min.pct = 0.2, logfc.threshold = 0.2)
dat.markers.0v7 <- FindMarkers(dat, ident.1= "0", ident.2="7",only.pos = F, min.pct = 0.2, logfc.threshold = 0.2)
#?FindAllMarkers#tweak
top10.dat.markers <- dat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

png(filename = "dat_marker_heatmap_10cluster.png", width = 1000, height = 1600, pointsize = 12,family = "sans", bg = "white")
DoHeatmap(dat, features = top10.dat.markers$gene,) + NoLegend() + scale_fill_viridis(option = "B")+ 
  theme(text = element_text(size = 20))
dev.off()

## rank and save the marker genes

dat.markers$rank <-(dat.markers$avg_logFC)*(-log(dat.markers$p_val_adj)) ## rank for GSEA analysis later
write.csv(dat.markers,"dat.DEG_genes.csv")

#webgestalt.org

## plot 2: key genes

png(filename = "epithelial cells feature.png", width = 10000, height = 10000, pointsize = 12,family = "sans", bg = "white",res=900)
FeaturePlot(dat, features = c("Slc34a1","Slc12a3","Umod","Aqp2","Lcn2","Slc12a1","Foxi1","Havcr1","Mki67"),order=T, pt.size = 1)
dev.off()

png(filename = "cells feature.png", width = 10000, height = 10000, pointsize = 12,family = "sans", bg = "white",res=900)
FeaturePlot(dat, features = c("Slc34a1","Umod","Tfcp2l1","Emcn","Cd3g","Gzma","Cx3cr1","Itgax","S100a9","Pdgfrb","Cd79a","Mki67"),order=T, pt.size = 1)
dev.off()


png(filename = "cells feature dot.png", width = 1000, height = 1000, pointsize = 12,family = "sans", bg = "white")
DotPlot(object = dat,features = c("Slc34a1","Umod","Tfcp2l1","Emcn","Cd3g","Gzma","Cx3cr1","Itgax","S100a9","Pdgfrb","Cd79a","Mki67"),dot.scale = 15) + coord_flip()
dev.off()

## rename ----


dat <- RenameIdents(object = dat, `0` = "PT1",`1` = "PT2",`2` = "PT3",`3` = "Loh/DCT",`4` = "collecting duct",
                     `5` = "Endothelial",`6` = "Macrophage",`7` = "Mesenchymal")

png(filename = "final_named_cluster.png", width = 1000, height = 1000, pointsize = 12,family = "sans", bg = "white")
DimPlot(dat, reduction = "umap",pt.size = 6, label = T,repel=T)
dev.off()

