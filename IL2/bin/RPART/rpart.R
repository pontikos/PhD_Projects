library(iterators)
source('~nikolas/Projects/IL2/bin/RPART/rpart-functions.R')
load('~/dunwich/Projects/IL2/transforms.RData')

DOSES <- c('0U','01U','10U','1000U')

# 
individual <- 'CB00086S'
date <- '2012-09-18'

# shit data, CD56 stain didn't work, need plot proving this.
#individual <- 'CB00406Q'
#date <- '2012-06-12'


#setwd('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/SPADE/CB00086S_2012-09-18/')
### Recursive partionning on side/forward scatter
MARKERS <- c('FSCA','SSCA')
print(load(sprintf('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/RData/CB00086S_%s_2012-09-18.RData','0U')))
fcs.data <- applyTransforms(fcs.data,transforms)
cs <- create.splits(fcs.data[,MARKERS],fun=var,branch.length=7)
# cluster stats
clusters <- c()
for (dose in DOSES) {
    print(load(sprintf('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/RData/CB00086S_%s_2012-09-18.RData',dose)))
    fcs.data <- applyTransforms(fcs.data,transforms)
    clusters <- rbind(clusters,cbind(cluster.stats(fcs.data[,c(MARKERS,'PSTAT5')],cs),dose=dose))
} 
print(clusters)


# Overlay splits on top of smooth scatter
pdf('~/Thesis/figures/scatter-rpart-128bin-var.pdf')
figure.labels <- iter(paste(letters,')',sep=''))
#BLTR
par(mfrow=c(2,2), mai=c(.75,.75,.5,.5))
for (dose in DOSES) {
    print(f <- sprintf('CB00086S_%s_2012-09-18.RData', dose))
    print(load(f))
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
par(mfrow=c(2,2), mai=c(.75,.75,.5,.5))
for (dose in DOSES) {
    print(f <- sprintf('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/RData/CB00086S_%s_2012-09-18.RData', dose))
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
    plot(NULL, xlim=range(fcs.data[,'FSCA']), ylim=range(fcs.data[,'SSCA']), xlab='FSCA', ylab='SSCA')
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
par(mfrow=c(2,2), mai=c(.75,.75,.5,.5))
for (dose in DOSES) {
    print(f <- sprintf('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/RData/CB00086S_%s_2012-09-18.RData', dose))
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
dev.off()

# cs.all: bins coloured by increase/decrease in prop between samples
pool.fcs.data <- c()
for (dose in DOSES) {
    print(load(sprintf('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/CLR/CB00086S_%s_2012-09-18.RData',dose)))
    print(load(sprintf('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/RData/CB00086S_%s_2012-09-18.RData',dose)))
    trans.fcs.data <- applyTransforms(fcs.data,transforms)
    pool.fcs.data <- rbind(pool.fcs.data, trans.fcs.data)
}
cs.all <- create.splits(pool.fcs.data[,MARKERS],fun=var,branch.length=7)
# cluster stats
clusters <- c()
for (dose in DOSES) {
    print(load(sprintf('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/RData/CB00086S_%s_2012-09-18.RData',dose)))
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
par(mfrow=c(2,2), mai=c(.75,.75,.5,.5))
for (dose in DOSES) {
    print(f <- sprintf('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/RData/CB00086S_%s_2012-09-18.RData', dose))
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
dev.off()



tapply(clusters$prop,clusters$cluster,function(x)x-x[1])
range(clusters$prop)

# how the prop distribution changes between samples
par(mfrow=c(1,1))
plot(density(clusters[which(clusters$dose=='10U'),'prop']),xlim=range(clusters$prop))
for (dose in DOSES[2:4]) lines(density(clusters[which(clusters$dose==dose),'prop']))
abline(v=clusters[which(clusters$dose=='0U'),'prop'])

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
print(load('CB00086S_0U_2012-09-18.RData'))
fcs.data <- transform.scatter(fcs.data)
cs2 <- apply.splits(fcs.data[,MARKERS], cs$splits)
print(load('CB00086S_0U_2012-09-18.RData'))
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



plot(tapply(clusters$PSTAT5,clusters$cluster,function(x)max(x)-min(x)),tapply(clusters$CD25,clusters$cluster,function(x)mean(x)))

library(ash)
b <- bin2(as.matrix(fcs.data[,c('FSCA','SSCA')]),nbin=c(512,512))
b <- ash2(b)

hist(clusters[which(clusters$dose=='0U'),'prop'])

smoothPlot(Y1[,1:2], classification=Y1$cluster, chull.lwd=2, ellipses=FALSE)
smoothPlot(Y2[,1:2], classification=Y2$cluster, chull.lwd=2, ellipses=FALSE)

setwd('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/SPADE/CB00086S_2012-09-18/')

print(load(list.files(pattern='.*_0U_.*.RData')))
X2 <- fcs.data[,c('FSCA','SSCA')]
Y2 <- apply.splits(X2)

chisq.test(table(Y1$cluster),table(Y2$cluster))

#smoothPlot(X) 
#which.max(apply(X,2,function(x)max(x)-min(x))) 
#apply(X,2,median)

setwd('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/SPADE/CB00086S_2012-09-18/')

x <- c()
for (dose in DOSES) {
    #sapply(DOSES, function(dose) sprintf('CB00086S_%s_2012-09-18.RData', dose))) {
    print(f <- sprintf('CB00086S_%s_2012-09-18.RData', dose))
    load(f)
    x <- rbind(x,cbind(fcs.data,dose=dose))
}

ggplot(x,aes(x=FSCA,y=SSCA,col=dose))+geom_density2d()

ggplot(X1,aes(x=FSCA,y=SSCA))+geom_density2d(n=10000)+geom_density2d(data=X2,col='red',n=10000)

smoothPlot(Y[,1:2])
Y <- rpart(X,fun=var,branch.length=4)
#smoothPlot(Y[,1:2],classification=1+Y[,3],chull.lwd=2,ellipses=FALSE)
smoothPlot(Y[,1:2])
Y <- rpart(X,fun=range,branch.length=4)
smoothPlot(Y[,1:2],classification=1+Y[,3],chull.lwd=2,ellipses=FALSE)

by(Y[,1:2], Y[,3], function(x) rect(min(x[,1]),min(x[,2]),max(x[,1]),max(x[,2])))
smoothPlot(Y[,1:2])

## Lymphocytes/single cells which are not CD4+
MARKERS <- c('CD25','CD3','CD4','CD45RA','CD56','CD8','FOXP3')
CELL.TYPES <- c("Lymphocytes", "Single cells", "CD4", "Memory", "Naive", "Naive Eff", "Naive Treg", "Memory Eff", "Memory Treg")
pool.fcs.data <- c()
for (dose in DOSES) {
    print(load(sprintf('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/CLR/CB00086S_%s_2012-09-18.RData',dose)))
    print(load(sprintf('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/RData/CB00086S_%s_2012-09-18.RData',dose)))
    #fcs.data <- fcs.data[as.logical(CLR[,'Single cells']),]
    #print(dim(fcs.data <- fcs.data[which(apply(CLR[,CELL.TYPES],1,function(x) strtoi(paste(x,collapse=''),base=2))==strtoi('110000000',base=2)),]))
    print(dim(fcs.data <- fcs.data[which(apply(CLR[,CELL.TYPES],1,function(x) sum(2**(0:8)[as.logical(x)]))==3),]))
    trans.fcs.data <- applyTransforms(fcs.data,transforms)
    pool.fcs.data <- rbind(pool.fcs.data, trans.fcs.data)
}
cs <- create.splits(pool.fcs.data[,MARKERS], branch.length=10) 
clusters <- c()
for (dose in DOSES) {
    print(load(sprintf('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/CLR/CB00086S_%s_2012-09-18.RData',dose)))
    print(load(sprintf('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/RData/CB00086S_%s_2012-09-18.RData',dose)))
    #print(dim(fcs.data <- fcs.data[which(apply(CLR[,CELL.TYPES],1,function(x) strtoi(paste(x,collapse=''),base=2))==strtoi('110000000',base=2)),]))
    print(dim(fcs.data <- fcs.data[which(apply(CLR[,CELL.TYPES],1,function(x) sum(2**(0:8)[as.logical(x)]))==3),]))
    trans.fcs.data <- applyTransforms(fcs.data,transforms)
    clusters <- rbind(clusters,cbind(cluster.stats(trans.fcs.data[,c(MARKERS,'PSTAT5')],cs),dose=dose))
}


# bin medians 
X <- do.call( 'rbind',by( clusters[,MARKERS] ,clusters$cluster, colMedians ) )
PSTAT5 <- do.call( 'rbind',by( clusters, clusters$cluster, function(X) X[,'PSTAT5']-X[which(X$dose=='0U'),'PSTAT5']) )
colnames(PSTAT5) <- DOSES
X <- cbind(X,PSTAT5)

# MST
mst.layout <- compute.mst.layout(X[,MARKERS])
palette <- colorRampPalette(c("blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red"))
colorscale <- palette(100)
boundaries <- quantile(X[,'1000U'],probs=c(.02,.98))
grad <- seq(boundaries[[1]],boundaries[[2]],length.out=length(colorscale)) 

#plot mst from X
pdf('~/Thesis/figures/rpart-lymphocytes-mst-1024bin.pdf',width=10,height=10)
par(mfrow=c(2,2))
figure.labels <- iter(paste(letters,')',sep=''))
for (dose in DOSES) {
    plot.mst(mst.layout, col=colorscale[findInterval(X[,dose], grad, all.inside=TRUE)]) 
    title(paste(nextElem(figure.labels),dose, sep='\t'), adj=0)
}
#dput(G1 <- locator(type='l',lwd=2, col='pink'))
G1 <- structure(list(x = c(-4503.23657353247, -4123.86218348976, -3768.96356054659, -3499.73012245176, -3389.58917050388, -3585.39530730011, -3891.34239604423, -4368.61985448506, -4552.18810773153), y = c(-3686.44502978582, -3677.95890775179, -3991.94542301115, -4305.93193827051, -4475.65437895125, -4611.43233149584, -4755.69640607447, -4475.65437895125, -3686.44502978582)), .Names = c("x", "y"))
lines(G1$x,G1$y,lwd=2, col='pink')
#dput(G2 <- locator(type='l',lwd=2, col='purple'))
G2 <- structure(list(x = c(2802.7799056771, 3120.96487797098, 3390.19831606581, 3439.14985026487, 3280.05736411793, 3157.67852862028, 3182.15429571981, 3316.77101476722, 3439.14985026487, 3365.72254896628, 3194.39217926957, 2839.49355632639), y = c(-5748.57268405678, -5544.90575523989, -5468.53065693356, -5646.73921964834, -5723.11431795467, -5740.08656202274, -5807.97553829504, -6003.15634507788, -6232.28163999688, -6334.11510440532, -6334.11510440532, -5765.54492812485)), .Names = c("x", "y"))
lines(G2$x,G2$y,lwd=2, col='purple')
dev.off()

k1 <- as.numeric(rownames(X)[which(as.logical(point.in.polygon(mst.layout[,1],mst.layout[,2],G1$x,G1$y)))])
k2 <- as.numeric(rownames(X)[which(as.logical(point.in.polygon(mst.layout[,1],mst.layout[,2],G2$x,G2$y)))])
#k2 <- names(which(X[as.character(k),'10U'] >= quantile(X[as.character(k),'10U'],.99)))
#k2 <- c(k2, '4089')


base.dir <- '~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5'
cell.types <- c("Memory Eff", "Memory Treg", "Naive Eff", "Naive Treg")

#dose response
PSTAT5 <- c()
for (dose in c('0U','01U','10U','1000U')) {
print(load(file.path(base.dir,'RData',f <- sprintf("%s.RData",paste(individual,dose,date,sep='_')))))
print(load(file.path(base.dir,'CLR',f)))
fcs.data <- applyTransforms(fcs.data[which(as.logical(CLR[,'Single cells'])),],transforms)
CLR <- CLR[which(as.logical(CLR[,'Single cells'])),]
csx <- apply.splits(as.data.frame(fcs.data[,MARKERS]),cs$splits)
CLR <- cbind(CLR,purple=0)
CLR[,'purple'] <- as.numeric(csx$X$cluster %in% k2 & (rowSums(CLR[,cell.types])==0))
CLR <- cbind(CLR,pink=0)
CLR[,'pink'] <- as.numeric(csx$X$cluster %in% k1 & (rowSums(CLR[,cell.types])==0))
PSTAT5 <- cbind(PSTAT5,sapply( c(CELL.TYPES,'purple','pink'), function(cell.type) median(fcs.data[as.logical(CLR[,cell.type]),'PSTAT5']) ) )
} 

pdf('~/Thesis/figures/rpart-lymphocytes-cd56bright.pdf',height=10,width=5)
par(mfrow=c(4,2))
cols <- c('black','red','darkgreen','blue','purple','pink')
for (marker in MARKERS)
smoothPlot1D(fcs.data[,marker], outliers=TRUE, posteriors=CLR[,c(cell.types,'purple','pink')], clusters.col=cols, main=marker, ellipse.lwd=3) 
ylim <- range(sapply(c(cell.types,'purple'), function(cell.type) PSTAT5[cell.type,]))
plot(NULL, xlim=c(1,4), ylim=ylim, xaxt='n', xlab='dose', ylab='pSTAT5 MFI', main='dose response')
#title(nextElem(figure.labels), adj=0)
axis(1, at=1:4, labels=DOSES)
i <- 1
for (cell.type in c(cell.types,'purple','pink')) {
    lines(1:4,  PSTAT5[cell.type,] , col=cols[[i]], lwd=3)
    i <- i+1
} 
dev.off()



pdf('~/Thesis/figures/rpart-lymphocytes-cd56bright.pdf',width=10,height=10)
plotClusters(all.trans.fcs.data[as.logical(CLR[,'Single cells']),MARKERS],
             #plot.points.col='purple',plot.points=trans.fcs.data[which(cs0U$X$cluster %in% k2),],
             outliers=TRUE, posteriors=CLR[as.logical(CLR[,'Single cells']),c('Memory Eff','Memory Treg','Naive Eff','Naive Treg','purple','pink')],chulls=FALSE,clusters.col=c('black','red','green','blue','purple','pink'))
dev.off()




## Non-lymphocytes
MARKERS <- c('SSCA','FSCA','CD25','CD3','CD4','CD45RA','CD56','CD8','FOXP3')
pool.fcs.data <- c()
for (dose in DOSES) {
    print(load(sprintf('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/CLR/CB00086S_%s_2012-09-18.RData',dose)))
    print(load(sprintf('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/RData/CB00086S_%s_2012-09-18.RData',dose)))
    fcs.data <- fcs.data[!as.logical(CLR[,'Lymphocytes']),]
    trans.fcs.data <- transform.scatter(applyTransforms(fcs.data,transforms)) 
    pool.fcs.data <- rbind(pool.fcs.data, trans.fcs.data)
}
cs <- create.splits(pool.fcs.data[,MARKERS], branch.length=10)
clusters <- c()
for (dose in DOSES) {
    print(load(sprintf('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/CLR/CB00086S_%s_2012-09-18.RData',dose)))
    print(load(sprintf('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/RData/CB00086S_%s_2012-09-18.RData',dose)))
    fcs.data <- fcs.data[!as.logical(CLR[,'Lymphocytes']),]
    trans.fcs.data <- transform.scatter(applyTransforms(fcs.data,transforms)) 
    clusters <- rbind(clusters,cbind(cluster.stats(trans.fcs.data[,c(MARKERS,'PSTAT5')],cs),dose=dose))
}

# bin medians 
X <- do.call( 'rbind',by( clusters[,MARKERS] ,clusters$cluster, colMedians ) )
PSTAT5 <- do.call( 'rbind',by( clusters, clusters$cluster, function(X) X[,'PSTAT5']-X[which(X$dose=='0U'),'PSTAT5']) )
colnames(PSTAT5) <- DOSES
X <- cbind(X,PSTAT5)

which(X[,'1000U'] > 1.2)
plot(density(X[,'1000U']))
lines(density(X[,'01U']))

plot(density(X[,'1000U']),xlim=c(boundaries[[1]],boundaries[[2]]))
gradient.rect(boundaries[[1]],0,boundaries[[2]],20,col=colorscale)
lines(density(X[,'1000U']),lwd=2)


pdf('~/Thesis/figures/rpart-nonlymphocytes-mst-1024bin.pdf')
par(mfrow=c(2,2))
figure.labels <- iter(paste(letters,')',sep=''))
for (dose in DOSES) {
    plot.mst(mst.layout, col=colorscale[findInterval(X[,dose], grad, all.inside=TRUE)]) 
    title(paste(nextElem(figure.labels),dose, sep='\t'), adj=0)
}
#purple cells in non-lymphocytes
#G2 <- locator(col='purple',lwd=2,type='l') 
lines(G2$x,G2$y,lwd=2,col='purple')
k <- as.numeric(rownames(X)[which(as.logical(point.in.polygon(mst.layout[,1],mst.layout[,2],G2$x,G2$y)))])
k2 <- names(which(X[as.character(k),'10U'] >= quantile(X[as.character(k),'10U'],.99)))
k2 <- c(k2, '4089')

#purple cells in lymphocytes
#G <- locator(col='purple',lwd=2,type='l') 
k <- as.numeric(rownames(X)[which(as.logical(point.in.polygon(mst.layout[,1],mst.layout[,2],G$x,G$y)))])
#lines(G$x,G$y,lwd=2,col='purple')
k2 <- names(which(X[as.character(k),'10U'] >= quantile(X[as.character(k),'10U'],.5)))
points(mst.layout[which(rownames(X) %in% k2),],col='purple',pch=20)
#dev.off()
identify(mst.layout)


#\myfigure{scale=.7}{rpart-lymphocytes-cd56bright}
pdf('~/Thesis/figures/rpart-nonlymphocytes-response.pdf')
plotClusters(all.trans.fcs.data[,c('SSCA','FSCA','CD4')],plot.points.pch='.', plot.points.col='purple',plot.points=trans.fcs.data[which(cs0U$X$cluster %in% k),],outliers=FALSE,chulls=FALSE,classification=CLR[,'Lymphocytes'])
dev.off()

### different MST layout?
adjacency  <- as.matrix(dist(X, method='manhattan'))
full_graph <- graph.adjacency(adjacency,mode="undirected",weighted=TRUE)
mst_graph  <- minimum.spanning.tree(full_graph)
m <- mst_graph
V(m)
get.diameter(m)
l <- layout.svd(m)


# PCA
PSTAT5 <- do.call( 'rbind',by( clusters, clusters$cluster, function(X) X[,'PSTAT5']-X[which(X$dose=='0U'),'PSTAT5']) )
colnames(PSTAT5) <- DOSES
X <- cbind(X,PSTAT5) 
palette <- colorRampPalette(c("blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red"))
colorscale <- palette(100)
boundaries <- quantile(X[,'1000U'],probs=c(.02,.98))
grad <- seq(boundaries[[1]],boundaries[[2]],length.out=length(colorscale)) 
pca <- princomp(X[,MARKERS])
par(mfrow=c(2,2))
figure.labels <- iter(paste(letters,')',sep=''))
for (dose in DOSES) {
    plot(pca$scores[,1:2], col=colorscale[findInterval(X[,dose], grad, all.inside=TRUE)], pch=20)
    title(paste(nextElem(figure.labels),dose, sep='\t'), adj=0)
}

# PLS
library(pls)
X <- do.call( 'rbind',by( clusters[,MARKERS] ,clusters$cluster, colMedians ) )
colnames(PSTAT5) <- paste('PSTAT5',1:4,sep='.')
X <- as.data.frame(cbind(X,PSTAT5))
par(mfrow=c(2,2))
figure.labels <- iter(paste(letters,')',sep=''))
for (pstat5 in paste('PSTAT5',2:4,sep='.')) {
    fm <- as.formula(sprintf('%s ~ %s', pstat5, paste(MARKERS, collapse='+')))
    pls <- plsr(fm , data=X)
    #plot(pls)
    #plot(pls,plottype='scores',comps=1:4,pch=20,col=colorscale[findInterval(X[,'PSTAT5.4'], grad, all.inside=TRUE)])
    plot(pls$scores[,1:2],pch=20,col=colorscale[findInterval(X[,pstat5], grad, all.inside=TRUE)])
    #title(paste(nextElem(figure.labels),dose, sep='\t'), adj=0)
}

palette <- colorRampPalette(c("blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red"))
colorscale <- palette(100)
boundaries <- quantile(X[,'PSTAT5.4'],probs=c(.02,.98))
grad <- seq(boundaries[[1]],boundaries[[2]],length.out=length(colorscale)) 

plot(pls,plottype='scores',comps=1:4,pch=20,col=colorscale[findInterval(X[,'PSTAT5.4'], grad, all.inside=TRUE)])

X[1,MARKERS] * pls$loadings[,1]
X[1,MARKERS] * pls$loading.weights[,1]

pls$scores[1,]









