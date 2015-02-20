#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tree))
suppressPackageStartupMessages(library(rpart))
#suppressPackageStartupMessages(library(mvpart))
suppressPackageStartupMessages(library(partykit))
suppressPackageStartupMessages(library(diptest))
source('~nikolas/bin/FCS/fcs.R',chdir=T)

option_list <- list( 
make_option(c("--in.file"), default=NULL, help = ".RData file"),
make_option(c('--channels'), default='FSCW,SSCW,FSCH,SSCH,FSCA,SSCA,CD4,CD25,CD45RA,FOXP3', help=''),
make_option(c("--cell.subset"), default=NULL, help = "any of the subsets defined in the CLR file"),
make_option(c("--clr"), default='magnetic-manual-gates', help = "any of the subsets defined in the CLR file"),
make_option(c("--out.dir"), default='~/Plots-tree-pstat5/', help = "Output directory"),
make_option(c("--splits"), default=NULL, help = "number of splits")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

if (!is.null(opt$channels)) {
    channels <- unlist(strsplit(opt$channels, ","))
} else  {
    channels <- NULL
}

base.dir <- '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All'
file.name <- 'CB01477E_2012-09-26.RData'
file.name <- opt$in.file
file.name <- file.path(base.dir,file.name)
print(basename(file.name))
clr <- 'magnetic-manual-gates2'
clr <- opt$clr
gate.dir <- sprintf('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/%s/CLR/',clr)
load(file.path(gate.dir,basename(file.name)))

#splits <- 3
splits <- opt$splits
#cell.subset <- 'CD4'
cell.subset <- opt$cell.subset

print(load(file.name)) 
load('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/transform-w1.RData')
fcs.data <- applyTransforms(fcs.data,transforms)
fcs.data <- baseline.relative.pstat5(fcs.data) 
#t <- tree::tree('diff.PSTAT5.4 ~ FSCA + SSCA + CD45RA + CD25 + FOXP3', data=data.frame(fcs.data))
t <- tree::tree('diff.PSTAT5.4 ~ FSCA + SSCA + CD4', data=data.frame(fcs.data))
#t <- prune.tree(t, best=3)
w <- as.factor(t$where)
levels(w) <- 1:length(unique(w))
table(w) 
plotClusters(fcs.data[,c('diff.PSTAT5.4','FSCA','SSCA','CD4')],classification=w,ellipses=FALSE,chull.lwd=2)

cd4<-which(as.logical(CLR[,'CD4'])) 
plotClusters( fcs.data[cd4,c('CD45RA','CD25','FOXP3')], classification=w[cd4], ellipses=FALSE,chull.lwd=2)

memory<-which(as.logical(CLR[,'Memory']))
plotClusters( fcs.data[memory,c('diff.PSTAT5.4','CD45RA','CD25','FOXP3')], classification=w[memory], ellipses=FALSE,chull.lwd=2,outliers=TRUE)

naive<-which(as.logical(CLR[,'Naive']))

plotClusters( fcs.data[,c('diff.PSTAT5.4','SSCA','FSCA')], classification=as.numeric(w==9), ellipses=FALSE,chull.lwd=2)


#lymphocytes
d2 <- data.frame(fcs.data[which(w==2),])

#for some reason CD4+ are split in 2!
#t <- tree::tree('diff.PSTAT5.4 ~ CD4 + CD25', data=d2)
t <- tree::tree('diff.PSTAT5.3 ~ CD4 + SSCA', data=d2)
t <- prune.tree(t, best=2)
w <- as.factor(t$where)
levels(w) <- 1:length(unique(w))
table(w) 
#plotClusters(d2[,c('diff.PSTAT5.4','CD4','SSCA','CD25')],classification=w,ellipses=FALSE,chull.lwd=2)
plotClusters(d2[,c('diff.PSTAT5.4','CD4','SSCA','CD25')],classification=w,chull.lwd=2)

#t <- tree::tree('diff.PSTAT5.4 ~ CD45RA + CD25 + CD4', data=data.frame(fcs.data[which(w==1),]))
#w <- as.factor(t$where)
#levels(w) <- 1:length(unique(w))
#table(w) 
#plotClusters(fcs.data[which(w==1),c('CD45RA','CD25','CD4')],classification=w,ellipses=FALSE,chull.lwd=2)

d22 <- data.frame(d2[which(w==2),])

t <- tree::tree('diff.PSTAT5.3 ~ CD45RA + SSCA', data=d22)
t <- prune.tree(t, best=2)
w <- as.factor(t$where)
levels(w) <- 1:length(unique(w))
table(w) 
plotClusters(d22[,c('diff.PSTAT5.4','SSCA','CD45RA')],classification=w,ellipses=FALSE,chull.lwd=2)

d221 <- d22[which(w==1),]
d222 <- d22[which(w==2),]

#memory
t <- tree::tree('diff.PSTAT5.2 ~ CD25 + FOXP3', data=d221)
t <- prune.tree(t, best=3)
w <- as.factor(t$where)
levels(w) <- 1:length(unique(w))
table(w) 
plotClusters(d221[,c('diff.PSTAT5.4','CD25','FOXP3')],classification=w,ellipses=FALSE,chull.lwd=2)


#naive
t <- tree::tree('diff.PSTAT5.2 ~ CD25 + FOXP3', data=d222)
t <- prune.tree(t, best=3)
w <- as.factor(t$where)
levels(w) <- 1:length(unique(w))
table(w) 
plotClusters(d222[,c('diff.PSTAT5.4','CD25','FOXP3')],classification=w,ellipses=FALSE,chull.lwd=2)




