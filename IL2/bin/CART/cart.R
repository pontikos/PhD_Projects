library(tree)

library(spade)
library(iterators)
library(flowCore)
source('~nikolas/bin/FCS/fcs.R')
source('~nikolas/Projects/IL2/bin/common.R')
source('~nikolas/Projects/IL2/bin/SPADE/spade-functions.R')

# this one looks ok
individual <- 'CB00086S'
day <- '2012-09-18'

BASE.DIR <- '~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5'

setwd(BASE.DIR)

load(file.path(BASE.DIR,'RData','pstat5-join', sprintf('%s_%s.RData',individual,day)))
load('~/dunwich/Projects/IL2/transforms.RData')
fcs.data <- applyTransforms(fcs.data,transforms)
fcs.data <- baseline.relative.pstat5(fcs.data)
print(load(file.path(BASE.DIR,'CLR',sprintf('%s_0U_%s.RData',individual,day))))
# Lymphocyte, single cells which are not memory or naive
print(dim(fcs.data <- fcs.data[which(apply(CLR[,CELL.TYPES],1,function(x) sum(2**(0:8)[as.logical(x)]))==3),]))
MARKERS <- c('CD25','CD3','CD4','CD45RA','CD56','CD8','FOXP3')
pdf('~/Thesis/figures/cart-lymphocytes-trees-pstat5.pdf')
figure.labels <- iter(paste(letters,')',sep=''))
par(mfrow=c(2,2),width=10,height=10)
for (i in 2:4) {
    fm <- as.formula(paste(sprintf('diff.PSTAT5.%d',i), paste(MARKERS, collapse='+'), sep='~'))
    tr <- tree::tree(fm, data=as.data.frame(fcs.data))
    #tr <- rpart(fm, data=as.data.frame(fcs.data))
    plot(tr)
    text(tr)
    title(paste(nextElem(figure.labels),DOSES[[i]], sep='\t'), adj=0)
}
dev.off()


par(mfrow=c(2,2))
for (i in 2:4) {
    diff.pstat5 <- sprintf('diff.PSTAT5.%d',i)
    fm <- as.formula(paste(diff.pstat5, paste(MARKERS, collapse='+'), sep='~'))
    tr <- tree::tree(fm, data=as.data.frame(fcs.data))
    w <- as.factor(tr$where)
    levels(w) <- 1:length(unique(w))
    plot(normalised.density(fcs.data[,diff.pstat5]),xlim=c(-.5,2))
    for (k in sort(unique(w))) lines(normalised.density(fcs.data[which(w==k),diff.pstat5]),col=k)
}


figure.labels <- iter(paste(letters,')',sep=''))
for (i in 2:4) {
    diff.pstat5 <- sprintf('diff.PSTAT5.%d',i)
    fm <- as.formula(paste(diff.pstat5, paste(MARKERS, collapse='+'), sep='~'))
    tr <- tree::tree(fm, data=as.data.frame(fcs.data))
    w <- as.factor(tr$where)
    levels(w) <- 1:length(unique(w))
    plotClusters(fcs.data[,c(diff.pstat5,grep('leaf',unique(tr$frame$var),value=TRUE,invert=TRUE))],classification=w,ellipses=TRUE,chulls=FALSE,outliers=TRUE)
    title(paste(nextElem(figure.labels),dose, sep='\t'), adj=0)
}

# non-lymphocytes
load(file.path(BASE.DIR,'RData','pstat5-join', sprintf('%s_%s.RData',individual,day)))
load('~/dunwich/Projects/IL2/transforms.RData')
fcs.data <- transform.scatter(applyTransforms(fcs.data,transforms))
fcs.data <- baseline.relative.pstat5(fcs.data)
print(load(file.path(BASE.DIR,'CLR',sprintf('%s_0U_%s.RData',individual,day))))
# Not lymphocytes
print(dim(fcs.data <- fcs.data[!as.logical(CLR[,'Lymphocytes']),])) 
MARKERS <- c('SSCA','FSCA','CD25','CD3','CD4','CD45RA','CD56','CD8','FOXP3') 
figure.labels <- iter(paste(letters,')',sep=''))
pdf('~/Thesis/figures/cart-nonlymphocytes-trees-pstat5.pdf')
par(mfrow=c(2,2))
for (i in 2:4) {
    fm <- as.formula(paste(sprintf('diff.PSTAT5.%d',i), paste(MARKERS, collapse='+'), sep='~'))
    tr <- tree::tree(fm, data=as.data.frame(fcs.data))
    w <- as.factor(tr$where)
    levels(w) <- 1:length(unique(w))
    #tr <- rpart(fm, data=as.data.frame(fcs.data))
    plot(tr)
    text(tr)
    title(paste(nextElem(figure.labels),DOSES[[i]], sep='\t'), adj=0)
}
dev.off()

smoothPlot(fcs.data[,c('SSCA','FSCA')],classification=w,ellipses=FALSE,chull.lwd=2)

par(mfrow=c(2,2),width=10,height=10)
for (i in 2:4) {
    diff.pstat5 <- sprintf('diff.PSTAT5.%d',i)
    fm <- as.formula(paste(diff.pstat5, paste(MARKERS, collapse='+'), sep='~'))
    tr <- tree::tree(fm, data=as.data.frame(fcs.data))
    w <- as.factor(tr$where)
    levels(w) <- 1:length(unique(w))
    plot(normalised.density(fcs.data[,diff.pstat5]),xlim=c(-.5,2))
    for (k in sort(unique(w))) lines(normalised.density(fcs.data[which(w==k),diff.pstat5]),col=k)
}







