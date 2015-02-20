source('~nikolas/Projects/IL2/bin/common.R')
source('~nikolas/Projects/IL2/bin/RPART/rpart-functions.R')

################################################## common

DOSES <- DOSES_

#setwd('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/SPADE/CB00086S_2012-09-18/')
### Recursive partionning on side/forward scatter
MARKERS <- c('FSCA','SSCA')
print(load(file.path(base.dir,'RData',FILES[[1]])))
fcs.data <- applyTransforms(fcs.data,transforms)
cs <- create.splits(fcs.data[,MARKERS],fun=var,branch.length=7)
# cluster stats
clusters <- c()
for (i in 1:4) {
    print(load(file.path(base.dir,'RData',FILES[[i]])))
    fcs.data <- applyTransforms(fcs.data,transforms)
    clusters <- rbind(clusters,cbind(cluster.stats(fcs.data[,c(MARKERS,'PSTAT5')],cs),dose=DOSES[[i]]))
} 
print(clusters)

# Overlay splits on top of smooth scatter
pdf('~/Thesis/figures/scatter-rpart-128bin-var.pdf')
figure.labels <- iter(paste(letters,')',sep=''))
#BLTR
par(mfrow=c(2,2), mai=c(.75,.75,.5,.5))
for (dose in DOSES) {
    print(f <- file.path(base.dir,'RData',sprintf('CB00086S_%s_2012-09-18.RData', dose)))
    print(load(f))
    fcs.data <- applyTransforms(fcs.data,transforms)
    cs2 <- apply.splits(fcs.data[,MARKERS], cs$splits)
    K <- cs2$X$cluster
    smoothPlot(fcs.data[,c('FSCA','SSCA')],outliers=TRUE)
    for (k in unique(sort(K))) {
        X <- fcs.data[which(K==k),c('FSCA','SSCA')]
        x <- range(X[,1])
        y <- range(X[,2])
        rect(x[1],y[1],x[2],y[2])
    }
    title(paste(nextElem(figure.labels),dose, sep='\t'), adj=0)
}
dev.off()

# bins coloured by pSTAT5 response
palette <- colorRampPalette(c("blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red"))
colorscale <- palette(100)
boundaries <- quantile(tapply(clusters$PSTAT5,clusters$cluster,function(x)max(x)-min(x)),probs=c(.02,.98))
grad <- seq(boundaries[[1]],boundaries[[2]],length.out=length(colorscale)) 
pdf('~/Thesis/figures/scatter-rpart-128bin-var-pstat5.pdf')
figure.labels <- iter(paste(letters,')',sep=''))
#BLTR
layout(matrix(c(1,2,3,4,5,5),nrow=3,ncol=2,byrow=T), heights=c(5,5,1),widths=c(5,5))
#par(mai=c(.65,.65,.5,.5))
par(mai=c(.5,.6,.5,.5))
#par(mfrow=c(2,2), mai=c(.75,.75,.5,.5))
for (i in 1:4) {
    dose <- DOSES_[[i]]
    f <- file.path(base.dir,'RData',FILES[[i]])
    print(load(f))
    fcs.data <- applyTransforms(fcs.data,transforms)
    cs2 <- apply.splits(fcs.data[,MARKERS], cs$splits)
    K <- cs2$X$cluster
    cluster <- as.data.frame(do.call('rbind',by(fcs.data,K,function(x)colMedians(x[,c(MARKERS,'PSTAT5')]))))
    cluster$n <- table(K)
    cluster$prop <- prop.table(table(K))
    cluster$cluster <- as.numeric(rownames(cluster))
    cluster$dose <- dose
    #smoothPlot(fcs.data[,c('FSCA','SSCA')], outliers=FALSE)
    plot(NULL, xlim=range(fcs.data[,'FSCA']), ylim=range(fcs.data[,'SSCA']), xlab='FSCA', ylab='SSCA', cex.lab=2, cex.main=2)
    for (k in unique(sort(K))) {
        X <- fcs.data[which(K==k),c('FSCA','SSCA')]
        x <- range(X[,1])
        y <- range(X[,2])
        pstat5 <- cluster[which(cluster$cluster==k),'PSTAT5']-clusters[which(clusters$cluster==k&clusters$dose=='0U'),'PSTAT5']
        col <- colorscale[findInterval(pstat5, grad, all.inside=TRUE)]
        rect(x[1],y[1],x[2],y[2],col=col)
    }
    title(paste(nextElem(figure.labels),dose, sep='\t'), adj=0)
}
par(mai=c(.4,2,.15,2))
image( grad, 1, matrix(1:length(colorscale),ncol=1), col=colorscale, xlab='', ylab="", yaxt="n", main='pSTAT5 MFI')
dev.off()


# cs from first sample: bins coloured by increase/decrease in prop between samples
palette <- colorRampPalette(c('red','white','green'))
colorscale <- palette(10)
prop <- as.numeric(do.call('rbind',tapply(clusters$prop,clusters$cluster,function(x)x-x[1])))
boundaries <- quantile(prop,probs=c(0,.99))
grad <- seq(boundaries[[1]],boundaries[[2]],length.out=length(colorscale)) 
pdf('~/Thesis/figures/scatter-rpart-128bin-var-prop.pdf')
figure.labels <- iter(paste(letters,')',sep=''))
#BLTR
#par(mfrow=c(2,2), mai=c(.75,.75,.5,.5))
layout(matrix(c(1,2,3,4,5,5),nrow=3,ncol=2,byrow=T), heights=c(5,5,1),widths=c(5,5))
par(mai=c(.5,.5,.5,.5))
for (dose in DOSES_) {
    print(f <- file.path(base.dir,'RData',sprintf('CB00086S_%s_2012-09-18.RData', dose)))
    print(load(f))
    fcs.data <- applyTransforms(fcs.data,transforms)
    cs2 <- apply.splits(fcs.data[,c('FSCA','SSCA')], cs$splits)
    K <- cs2$X$cluster
    cluster <- as.data.frame(do.call('rbind',by(fcs.data,K,function(x)colMedians(x[,c(MARKERS,'PSTAT5')]))))
    cluster$n <- table(K)
    cluster$prop <- prop.table(table(K))
    cluster$cluster <- as.numeric(rownames(cluster))
    cluster$dose <- dose
    #smoothPlot(fcs.data[,c('FSCA','SSCA')],outliers=TRUE)
    plot(NULL, xlim=range(fcs.data[,'FSCA']), ylim=range(fcs.data[,'SSCA']), xlab='FSCA', ylab='SSCA')
    for (k in unique(sort(K))) {
        X <- fcs.data[which(K==k),c('FSCA','SSCA')]
        x <- range(X[,1])
        y <- range(X[,2])
        pstat5 <- cluster[which(cluster$cluster==k),'prop']-clusters[which(clusters$cluster==k&clusters$dose=='0U'),'prop']
        col <- colorscale[findInterval(pstat5, grad, all.inside=TRUE)]
        rect(x[1],y[1],x[2],y[2],col=col)
    }
    title(paste(nextElem(figure.labels),dose, sep='\t'), adj=0)
}
par(mai=c(.4,2,.15,2))
image( grad, 1, matrix(1:length(colorscale),ncol=1), col=colorscale, xlab='', ylab="", yaxt="n", main='change in proportion')
dev.off()

# cs.all: bins coloured by increase/decrease in prop between samples
pool.fcs.data <- c()
for (dose in DOSES_) {
    print(load(file.path(base.dir,'CLR',sprintf('CB00086S_%s_2012-09-18.RData',dose))))
    print(load(file.path(base.dir, 'RData', sprintf('CB00086S_%s_2012-09-18.RData',dose))))
    trans.fcs.data <- applyTransforms(fcs.data,transforms)
    pool.fcs.data <- rbind(pool.fcs.data, trans.fcs.data)
}
cs.all <- create.splits(pool.fcs.data[,MARKERS],fun=var,branch.length=7)
# cluster stats
clusters <- c()
for (dose in DOSES_) {
    print(load(file.path(base.dir,'RData',sprintf('CB00086S_%s_2012-09-18.RData',dose))))
    fcs.data <- applyTransforms(fcs.data,transforms)
    clusters <- rbind(clusters,cbind(cluster.stats(fcs.data[,c(MARKERS,'PSTAT5')],cs.all),dose=dose))
} 
print(clusters) 

#
palette <- colorRampPalette(c('red','white','green'))
colorscale <- palette(10)
prop <- as.numeric(do.call('rbind',tapply(clusters$prop,clusters$cluster,function(x)x-mean(x))))
boundaries <- quantile(prop,probs=c(0,.99))
grad <- seq(boundaries[[1]],boundaries[[2]],length.out=length(colorscale)) 
pdf('~/Thesis/figures/scatter-rpartall-128bin-var-prop.pdf')
figure.labels <- iter(paste(letters,')',sep=''))
#BLTR
#par(mfrow=c(2,2), mai=c(.75,.75,.5,.5))
layout(matrix(c(1,2,3,4,5,5),nrow=3,ncol=2,byrow=T), heights=c(5,5,1),widths=c(5,5))
par(mar=c(.5,1,1,1),mai=c(.1,.5,.5,.5))
for (dose in DOSES) {
    print(f <- file.path(base.dir, 'RData', sprintf('CB00086S_%s_2012-09-18.RData', dose)))
    print(load(f))
    fcs.data <- applyTransforms(fcs.data,transforms)
    cs2 <- apply.splits(fcs.data[,c('FSCA','SSCA')], cs.all$splits)
    K <- cs2$X$cluster
    cluster <- as.data.frame(do.call('rbind',by(fcs.data,K,function(x)colMedians(x[,c(MARKERS,'PSTAT5')]))))
    cluster$n <- table(K)
    cluster$prop <- prop.table(table(K))
    cluster$cluster <- as.numeric(rownames(cluster))
    cluster$dose <- dose
    #smoothPlot(fcs.data[,c('FSCA','SSCA')],outliers=TRUE)
    plot(NULL, xlim=range(fcs.data[,'FSCA']), ylim=range(fcs.data[,'SSCA']), xlab='FSCA', ylab='SSCA')
    for (k in unique(sort(K))) {
        X <- fcs.data[which(K==k),c('FSCA','SSCA')]
        x <- range(X[,1])
        y <- range(X[,2])
        pstat5 <- cluster[which(cluster$cluster==k),'prop']-clusters[which(clusters$cluster==k&clusters$dose=='0U'),'prop']
        col <- colorscale[findInterval(pstat5, grad, all.inside=TRUE)]
        rect(x[1],y[1],x[2],y[2],col=col)
    }
    title(paste(nextElem(figure.labels),dose, sep='\t'), adj=0)
}
par(mai=c(.4,2,.15,2))
image( grad, 1, matrix(1:length(colorscale),ncol=1), col=colorscale, xlab='', ylab="", yaxt="n", main='change in proportion')
dev.off()

#plot max pSTAT5 fold increase against average prop diff from resting sample
zoom <- 5
pdf('~/Thesis/figures/rpart-1024bin-var.pdf',width=2*zoom,height=zoom)
par(mfrow=c(1,2), mai=c(1,1,.5,.5))
figure.labels <- iter(paste(letters,')',sep=''))
x <- tapply(clusters$PSTAT5,clusters$cluster,function(x)max(x)-min(x))
y <- tapply(clusters$prop,clusters$cluster,function(x)max(abs(x-x[1])))
plot(x,y,xlab='pSTAT5 response',ylab='change in prop')
X<-cbind(x,y)
K<-as.numeric(names(which(X[,2]<.0006&X[,1]>1)))
points(X[as.character(K),],col='red',pch=20)
#G <- locator(type='l')
#K <- clusters$cluster[which(as.logical(point.in.polygon(x,y,G$x,G$y)))]
title(nextElem(figure.labels), adj=0) 
print(load(file.path(base.dir,'RData','CB00086S_0U_2012-09-18.RData')))
fcs.data <- applyTransforms(fcs.data,transforms)
cs2 <- apply.splits(fcs.data[,MARKERS], cs$splits)
print(load(file.path(base.dir,'RData','CB00086S_0U_2012-09-18.RData')))
#plot(NULL, xlim=range(fcs.data[,'FSCA']), ylim=range(fcs.data[,'SSCA']), xlab='FSCA', ylab='SSCA')
smoothPlot(fcs.data[,c('FSCA','SSCA')],outliers=TRUE)
title(nextElem(figure.labels), adj=0) 
for (k in K) {
    X <- fcs.data[which(cs2$X$cluster==k),c('FSCA','SSCA')]
    x <- range(X[,1])
    y <- range(X[,2])
    rect(x[1],y[1],x[2],y[2])
}
dev.off()


################################################## Lymphocytes
source('~nikolas/Projects/IL2/bin/RPART/thesis-lymphocytes.R')


################################################## Non.Lymphocytes
source('~nikolas/Projects/IL2/bin/RPART/thesis-nonlymphocytes.R')


