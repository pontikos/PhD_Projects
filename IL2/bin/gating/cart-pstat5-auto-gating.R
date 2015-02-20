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

file.name <- opt$in.file

base.dir <- '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All'
file.name <- 'CB01477E_2012-09-26.RData'
#file.name <- opt$in.file
file.name <- file.path(base.dir,file.name)
print(basename(file.name))
load(file.name)

load('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/transform-w1.RData')
fcs.data <- applyTransforms(fcs.data,transforms)
fcs.data <- baseline.relative.pstat5(fcs.data) 
clr <- 'magnetic-manual-gates2'
gate.dir <- sprintf('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/%s/CLR/',clr)
load(file.path(gate.dir,basename(file.name)))

fcs.data2 <- fcs.data 
fcs.data <- fcs.data[!as.logical(CLR[,'Lymphocytes']),]
t <- tree::tree( diff.PSTAT5.4 ~ FSCA + SSCA + CD4 + CD45RA + FOXP3 + CD25, data=data.frame(fcs.data) )
w <- as.factor(t$where)
levels(w) <- 1:length(unique(w))
#plotClusters( fcs.data[,c('SSCA','FSCA','CD4','CD45RA','FOXP3','CD25')], classification=w, ellipses=TRUE, chulls=FALSE)
plotClusters( fcs.data2[,c('SSCA','FSCA','CD4','CD45RA','FOXP3','CD25')], plot.points=fcs.data[which(w==5),], ellipses=TRUE, chulls=FALSE)
smoothPlot( fcs.data[,c('SSCA','FSCA')], classification=as.numeric(w==5), ellipses=TRUE, chulls=FALSE)
classification.to.ellipse( fcs.data2[,c('SSCA','FSCA')], w==5)

par(mfrow=c(2,5))
sapply(sort(unique(w)), function(i) {
    X <- fcs.data[which(w==i),]
    xrange <- range(X[,grep('^diff.PSTAT5',colnames(X))])
    plot(normalised.density(X[,'diff.PSTAT5.2']),xlim=xrange)
    lines(normalised.density(X[,'diff.PSTAT5.3']))
    lines(normalised.density(X[,'diff.PSTAT5.4']))
})


sapply(sort(unique(w)), function(i) 100*prop.table(table({as.numeric(w==i)+CLR[,'Lymphocytes']})))


title(nextElem(figure.labels), adj=0) 
lymphocytes <- points.in.ellipse(fcs.data[,c('SSCA','FSCA')], w==2) 
fcs.data <- fcs.data[as.logical(lymphocytes),]
#single cells
t <- tree::tree( diff.PSTAT5.4 ~ SSCH + SSCW, data=data.frame(fcs.data) )
t <- prune.tree( t, best=2 )
w <- as.factor(t$where)
levels(w) <- 1:length(unique(w))
p1 <- points.in.ellipse(fcs.data[,c('SSCH','SSCW')], w==1)
p2 <- points.in.ellipse(fcs.data[,c('SSCH','SSCW')], w==2)
single.cells <- as.logical(p1) | as.logical(p2) 
smoothPlot( fcs.data[,c('SSCH','SSCW')], classification=single.cells, ellipses=TRUE, chulls=FALSE)
title(nextElem(figure.labels), adj=0) 
fcs.data <- fcs.data[as.logical(single.cells),]
# CD4+ lymphocytes
t <- tree::tree( diff.PSTAT5.3 ~ CD4 + SSCA, data=data.frame(fcs.data) )
if (!'CD4' %in% t$frame$var)
t <- tree::tree( diff.PSTAT5.4 ~ CD4 + SSCA, data=data.frame(fcs.data) )
print(t)
cat('CD4\n')
cd4.positive <- unique(as.numeric(gsub('>|<','',t$frame[which(t$frame$var=='CD4'),'splits'])))
cd4 <- as.numeric(fcs.data[,'CD4'] > cd4.positive)
smoothPlot( fcs.data[,c('CD4','SSCA')], classification=cd4, ellipses=TRUE, chulls=FALSE)
title(nextElem(figure.labels), adj=0) 
fcs.data <- fcs.data[as.logical(cd4),]
# memory / naive
t <- tree::tree( diff.PSTAT5.4 ~ CD45RA + SSCA, data=data.frame(fcs.data) )
print(t)
cat('CD45RA\n')
print(cd45ra.positive <- (min(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='CD45RA'),'splits'])))))
memory <- as.numeric(fcs.data[,'CD45RA']<cd45ra.positive)
print(cd45ra.positive <- (max(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='CD45RA'),'splits'])))))
naive <- as.numeric(fcs.data[,'CD45RA']>cd45ra.positive)
smoothPlot( fcs.data[,c('CD45RA','SSCA')], classification=memory+naive*2, ellipses=TRUE, chulls=FALSE)
title(nextElem(figure.labels), adj=0)
naive <- fcs.data[as.logical(naive),]
memory <- fcs.data[as.logical(memory),]
# naive effector / treg
t <- tree::tree( diff.PSTAT5.2 ~ CD25 + FOXP3, data=data.frame(naive) )
if (!'CD25' %in% t$frame$var)
t <- tree::tree( diff.PSTAT5.3 ~ CD25 + FOXP3, data=data.frame(naive) )
if (!'CD25' %in% t$frame$var)
t <- tree::tree( diff.PSTAT5.4 ~ CD25 + FOXP3, data=data.frame(naive) )
cat('Naive\n')
print(t)
# naive eff
print(foxp3.positive <- min(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='FOXP3'),'splits'])))))
print(cd25.positive <- min(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='CD25'),'splits'])))))
naive.eff <- as.numeric(naive[,'FOXP3'] < foxp3.positive & naive[,'CD25'] < cd25.positive)
# naive treg
print(foxp3.positive <- max(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='FOXP3'),'splits'])))))
print(cd25.positive <- max(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='CD25'),'splits'])))))
naive.tregs <- as.numeric(naive[,'FOXP3'] > foxp3.positive & naive[,'CD25'] > cd25.positive)
smoothPlot( naive[,c('CD25','FOXP3')], classification=naive.eff+naive.tregs*2, ellipses=TRUE, chulls=FALSE, outliers=TRUE )
title(nextElem(figure.labels), adj=0)
# memory effector / treg
t <- tree::tree( diff.PSTAT5.2 ~ CD25 + FOXP3, data=data.frame(memory) )
cat('Memory\n')
print(t)
# memory eff
print(foxp3.positive <- min(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='FOXP3'),'splits'])))))
print(cd25.positive <- min(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='CD25'),'splits'])))))
memory.eff <- as.numeric(memory[,'FOXP3'] < foxp3.positive & memory[,'CD25'] < cd25.positive)
# memory treg
print(foxp3.positive <- max(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='FOXP3'),'splits'])))))
print(cd25.positive <- max(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='CD25'),'splits'])))))
memory.tregs <- as.numeric(memory[,'FOXP3'] > foxp3.positive & memory[,'CD25'] > cd25.positive)
smoothPlot( memory[,c('CD25','FOXP3')], classification=memory.eff+memory.tregs*2, ellipses=TRUE, chulls=FALSE, outliers=TRUE )
title(nextElem(figure.labels), adj=0)
dev.off()






