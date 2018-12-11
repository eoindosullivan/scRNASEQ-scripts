
## Load librarys--------------------------------------------------------------------
  
library(dplyr)  
library(scde)


## Old flexmix and scde for avoid errors--------------------------------------------------------------------
  
install.packages("rlang",type="win.binary") 
install.packages("devtools")
require(devtools)
install_version("flexmix", version = "2.3-13", repos = "http://cran.us.r-project.org")

## run on desktop  (single core only in current stable platform)

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("scde")
require(devtools)
devtools::install_github('hms-dbmi/scde', build_vignettes = FALSE)

## older stable version of scde here

install.packages("scde-1.99.2.tar.gz", repos = NULL,
                 type="source")

## Raw Object prep--------------------------------------------------------------------
  
## make the cell names the sorting variable, i.e. the injury
testdat<-read.csv("SCDERAWDATA.csv",header=T)

rawdata<-readRDS("senepiirirawdat") ## this way is much more convienient
as.data.frame(rawdata)->rawdata
rownames(rawdata)
colnames(rawdata)
#rownames(rawdata)<-rawdata[,1]
#dim(rawdata)
#rawdata <- Octmacall16@raw.data[1:576]
#colnames(rawdata) <- make.unique(as.character(Octmacall16@meta.data$injury), sep = "_")
## try loop

#a<-c("D2UUO","D2UUO","D7UUO","R-UUO")
#testdataframe<-data.frame(a,a)


 
## Dataprep--------------------------------------------------------------------

## extract the groups you want

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
## this takes absofuckinglutlyages ( this took me about 70 second per cell) unless use multiple cores ( ld version 1.99 doesnt support >1 core on windows, this one is one to run in the evenings total pain.)
## Way faster with 10x data actuially.
## fit the models

o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)

head(o.ifm)

## save object, takes ages
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


## Dataprep--------------------------------------------------------------------

## extract the groups you want

df<-NULL
df <- rawdata %>% dplyr:: select(starts_with(c("D2UUO")))
df<-cbind(df, rawdata %>% dplyr:: select(starts_with(c("D7UUO"))))

row.names(df)
# clean up the dataset
cd <- clean.counts(df, min.lib.size=1000, min.reads = 1, min.detected = 1)
sg <- factor(gsub("(D2UUO|D7UUO).*", "\\1", colnames(cd)), levels = c("D2UUO", "D7UUO"))
# clean up the dataset
# the group factor should be named accordingly
names(sg) <- colnames(cd)
names(sg)==colnames(cd)
table(sg)
rownames(cd)
## if false then remove 
# use this line to create integers rather than characters in the counts matrix
cd <-apply(cd,2,function(x) {storage.mode(x) <- 'integer'; x}) 

## Error models--------------------------------------------------------------------

## calculate the error models 
## this takes absofuckinglutlyages ( this took me about 70 second per cell) unless use multiple cores ( ld version 1.99 doesnt support >1 core on windows, this one is one to run in the evenings total pain.)

## fit the models

o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)

head(o.ifm)

## save object, takes ages
o.ifm->o.ifmd2d7


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
groups <- factor(gsub("(D2UUO|D7UUO).*", "\\1", rownames(o.ifm)), levels  =  c("D2UUO", "D7UUO"))
names(groups) <- row.names(o.ifm)

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)


# top upregulated genes (tail would show top downregulated ones)
head(ediff[order(ediff$Z, decreasing  =  TRUE), ])


# write out a table with all the results, showing most significantly different genes (in both directions) on top
write.table(ediff[order(abs(ediff$Z), decreasing = TRUE), ], file = "SCDE_D2_D7_results.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

## Dataprep--------------------------------------------------------------------

## extract the groups you want

df<-NULL
df <- rawdata %>% dplyr:: select(starts_with(c("D2UUO")))
df<-cbind(df, rawdata %>% dplyr:: select(starts_with(c("R.UUO"))))

row.names(df)
# clean up the dataset
cd <- clean.counts(df, min.lib.size=1000, min.reads = 1, min.detected = 1)
sg <- factor(gsub("(D2UUO|R.UUO).*", "\\1", colnames(cd)), levels = c("D2UUO", "R.UUO"))
# clean up the dataset
# the group factor should be named accordingly
names(sg) <- colnames(cd)
names(sg)==colnames(cd)
table(sg)
rownames(cd)
## if false then remove 
# use this line to create integers rather than characters in the counts matrix
cd <-apply(cd,2,function(x) {storage.mode(x) <- 'integer'; x}) 

## Error models--------------------------------------------------------------------

## calculate the error models 
## this takes absofuckinglutlyages ( this took me about 70 second per cell) unless use multiple cores ( ld version 1.99 doesnt support >1 core on windows, this one is one to run in the evenings total pain.)

## fit the models

o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)

head(o.ifm)

## save object, takes ages
o.ifm->o.ifmd2ruuo


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
groups <- factor(gsub("(D2UUO|R.UUO).*", "\\1", rownames(o.ifm)), levels  =  c("D2UUO", "R.UUO"))
names(groups) <- row.names(o.ifm)

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)


# top upregulated genes (tail would show top downregulated ones)
head(ediff[order(ediff$Z, decreasing  =  TRUE), ])


# write out a table with all the results, showing most significantly different genes (in both directions) on top
write.table(ediff[order(abs(ediff$Z), decreasing = TRUE), ], file = "SCDE_D2_RUUO_results.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

## Dataprep--------------------------------------------------------------------

## extract the groups you want

df<-NULL
df <- rawdata %>% dplyr:: select(starts_with(c("D7UUO")))
df<-cbind(df, rawdata %>% dplyr:: select(starts_with(c("R.UUO"))))

row.names(df)
# clean up the dataset
cd <- clean.counts(df, min.lib.size=1000, min.reads = 1, min.detected = 1)
sg <- factor(gsub("(D7UUO|R.UUO).*", "\\1", colnames(cd)), levels = c("D7UUO", "R.UUO"))
# clean up the dataset
# the group factor should be named accordingly
names(sg) <- colnames(cd)
names(sg)==colnames(cd)
table(sg)
rownames(cd)
## if false then remove 
# use this line to create integers rather than characters in the counts matrix
cd <-apply(cd,2,function(x) {storage.mode(x) <- 'integer'; x}) 

## Error models--------------------------------------------------------------------

## calculate the error models 
## this takes absofuckinglutlyages ( this took me about 70 second per cell) unless use multiple cores ( ld version 1.99 doesnt support >1 core on windows, this one is one to run in the evenings total pain.)

## fit the models

o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)

head(o.ifm)

## save object, takes ages
o.ifm->o.ifmd7ruuo


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
groups <- factor(gsub("(D7UUO|R.UUO).*", "\\1", rownames(o.ifm)), levels  =  c("D7UUO", "R.UUO"))
names(groups) <- row.names(o.ifm)

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)


# top upregulated genes (tail would show top downregulated ones)
head(ediff[order(ediff$Z, decreasing  =  TRUE), ])


# write out a table with all the results, showing most significantly different genes (in both directions) on top
write.table(ediff[order(abs(ediff$Z), decreasing = TRUE), ], file = "SCDE_D7_RUUO_results.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)


