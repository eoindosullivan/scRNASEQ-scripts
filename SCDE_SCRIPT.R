## __________________________________________
## Project: Multiple
## Script purpose: Differential expression
## Date: 14.10.18
## Author: Eoin O'Sullivan
_______________
## Load libraries--------------------------------------------------------------------

library(dplyr)  
library(scde)

## Old flexmix and scde needed to avoid errors--------------------------------------------------------------------
## not resolved for multiple people on forums
  
install.packages("rlang",type="win.binary") 
install.packages("devtools")
require(devtools)
install_version("flexmix", version = "2.3-13", repos = "http://cran.us.r-project.org")

## run on desktop rather than laptop for speed(single core only in current stable platform for windows machines as no forking of processes outside linux)

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("scde")
require(devtools)
devtools::install_github('hms-dbmi/scde', build_vignettes = FALSE)

## older stable version of scde here

install.packages("scde-1.99.2.tar.gz", repos = NULL,
                 type="source")

## Raw Object prep--------------------------------------------------------------------
  
rawdata<-readRDS("senepiirirawdat") ## this way is much more convienient to import than the raw matrix if running lots of iterations of the clusters
as.data.frame(rawdata)->rawdata

## when relevant make the cell names the sorting variable, i.e. the injury
#rownames(rawdata)<-rawdata[,1]
#dim(rawdata)
#rawdata <- Octmacall16@raw.data[1:576]
#colnames(rawdata) <- make.unique(as.character(Octmacall16@meta.data$injury), sep = "_")

## Dataprep--------------------------------------------------------------------

## extract the groups you want
## in this case cells prior indentified as "Senescent"

df<-NULL
df <- rawdata %>% dplyr:: select(starts_with(c("Senescent")))
df<-cbind(df, rawdata %>% dplyr:: select(starts_with(c("Non_Senescent"))))
#df[1,] ## just ensure the rows are genes, sometimes gets bumped to col 1
row.names(df)

# clean up the dataset
cd <- clean.counts(df, min.lib.size=1000, min.reads = 1, min.detected = 1)
sg <- factor(gsub("(Senescent|Non_Senescent).*", "\\1", colnames(cd)), levels = c("Senescent", "Non_Senescent"))

# clean up the dataset
# the group factor should be named accordingly
names(sg) <- colnames(cd)
names(sg)==colnames(cd) ## needs to be TRUE
table(sg) ## final checks
rownames(cd) ## dont fuck me 
## if false then remove 
# use this line to create integers rather than characters in the counts matrix
cd <-apply(cd,2,function(x) {storage.mode(x) <- 'integer'; x}) 

## Error models--------------------------------------------------------------------

## calculate the error models 
## this takes along time with SMARTSEQ data on laptop ( this took me about 70 second per cell) unless use multiple cores ( ld version 1.99 doesnt support >1 core on windows)
## Way faster with 10x data.
## would benefit from a formal script and job submission to a linux cluster to allow forking
## fit the models

o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)

## save object
o.ifm->o.ifmshamd7

# filter out cells that don't show positive correlation with
# the expected expression magnitudes (very poor fits)

valid.cells <- o.ifm$corr.a > 0
table(valid.cells)

o.ifm <- o.ifm[valid.cells, ]

# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)

## Differential expression --------------------------------------------------------------------

## actually check for differential expression now

# define two groups of cells
groups <- factor(gsub("(Senescent|Non_Senescent).*", "\\1", rownames(o.ifm)), levels  =  c("Senescent", "Non_Senescent"))
names(groups) <- row.names(o.ifm)

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)


# top upregulated genes (tail would show top downregulated ones)
head(ediff[order(ediff$Z, decreasing  =  TRUE), ])


# write out a table with all the results, showing most significantly different genes (in both directions) on top
write.table(ediff[order(abs(ediff$Z), decreasing = TRUE), ], file = "senepi_vs_normalepi_results.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
