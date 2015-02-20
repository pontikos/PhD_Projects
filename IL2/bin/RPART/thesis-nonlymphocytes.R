################################################## Non.Lymphocytes

## Non-lymphocytes

#core fluorescent markers (no scatter)
MARKERS <- CORE.MARKERS
pool.fcs.data <- c()
for (i in 1:4) {
    print(load(file.path(base.dir,'RData',FILES[[i]])))
    print(load(file.path(base.dir,'CLR',FILES[[i]])))
    fcs.data <- fcs.data[!as.logical(CLR[,'Lymphocytes']),]
    trans.fcs.data <- transform.scatter(applyTransforms(fcs.data,transforms)) 
    pool.fcs.data <- rbind(pool.fcs.data, trans.fcs.data)
}
cs <- create.splits(pool.fcs.data[,MARKERS], branch.length=10) 

clusters <- c()
for (i in 1:4) {
    print(load(file.path(base.dir,'RData',FILES[[i]])))
    print(load(file.path(base.dir,'CLR',FILES[[i]])))
    fcs.data <- fcs.data[!as.logical(CLR[,'Lymphocytes']),]
    trans.fcs.data <- transform.scatter(applyTransforms(fcs.data,transforms)) 
    clusters <- rbind(clusters,cbind(cluster.stats(trans.fcs.data[,c(MARKERS,'PSTAT5')],cs),dose=DOSES[[i]]))
}

# bin medians 
# necessary for the colouring of the tree
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

pdf('~/Thesis/figures/rpart-nonlymphocytes-mst-1024bin.pdf')
layout(matrix(c(1,2,3,4,5,5),nrow=3,ncol=2,byrow=T), heights=c(5,5,1),widths=c(5,5))
par(mar=c(.5,1,1,1),mai=c(.1,.5,.5,.5)) 
figure.labels <- iter(paste(letters,')',sep=''))
for (dose in DOSES) {
    plot.mst(mst.layout, col=colorscale[findInterval(X[,dose], grad, all.inside=TRUE)]) 
    title(paste(nextElem(figure.labels),dose, sep='\t'), adj=0)
}
#purple cells in non-lymphocytes
#dput(G1 <- locator(type='l',lwd=2, col='purple'))
G1 <- structure(list(x = c(837.438131018855, 1290.38416181167, 1250.99755043838, 896.518048078787, 837.438131018855), y = c(-8253.88312243777, -8067.69171038357, -8702.93299856847, -8757.6951785844, -8231.97825043139)), .Names = c("x", "y"))
lines(G1$x,G1$y,lwd=2,col='purple') 
par(mai=c(.4,2,.15,2))
image( grad, 1, matrix(1:length(colorscale),ncol=1), col=colorscale, xlab='', ylab="", yaxt="n", main='pSTAT5 MFI')
dev.off()

k <- as.numeric(rownames(X)[which(as.logical(point.in.polygon(mst.layout[,1],mst.layout[,2],G1$x,G1$y)))])

# CLR
print(load(file.path(base.dir,'RData',FILES[[1]])))
dim(fcs.data <- transform.scatter(applyTransforms(fcs.data,transforms)))
print(load(file.path(base.dir,'CLR',FILES[[1]])))
csx <- apply.splits(as.data.frame(fcs.data[which(!as.logical(CLR[,'Lymphocytes'])),MARKERS]),cs$splits)
CLR <- cbind(CLR,RPART.1=0)
CLR[which(!as.logical(CLR[,'Lymphocytes'])),'RPART.1'] <- as.numeric(csx$X$cluster %in% k)
#now remove all those which have already been included as part of the lymphocytes
CLR[which(CLR[,'Lymphocytes']>1),grep('RPART',colnames(CLR))]<-0 
#################### univariate clusters
pdf('~/Thesis/figures/rpart-nonlymphocytes-clusters.pdf')
#par(mfrow=c(3,3))
r <- c(1,1,2,2,3,3)
layout(rbind(r,3+r,6+r,rep(10,length(r))),heights=c(4,4,4,1)) 
par(mai=rep(0.5, 4))
for (marker in CORE.MARKERS)
smoothPlot1D(fcs.data[,marker], outliers=TRUE, posteriors=CLR[,c('Lymphocytes','RPART.1')], clusters.col=c('black','purple'), main=marker, col='white')
par(mai=c(0,0,0,0))
plot.new()
legend('center',ncol=2,legend=c('Lymphocytes','RPART.1'),fill=c('black','purple'), bty='n',cex=2)
dev.off()

#################### multivariate clusters 
pdf('~/Thesis/figures/rpart-nonlymphocytes-scatter-clusters.pdf')
plotClusters(fcs.data[,c('SSCA','FSCA','CD4')],plot.points.col='purple',plot.points.pch='.',plot.points=fcs.data[which(as.logical(CLR[,'RPART.1'])),],outliers=FALSE,chulls=FALSE,classification=CLR[,'Lymphocytes'])
dev.off()


#################### dose response
pdf('~/Thesis/figures/rpart-nonlymphocytes-dose-response.pdf') 
PSTAT5 <- c()
for (i in 1:4) {
print(load(file.path(base.dir,'RData',FILES[[i]])))
print(load(file.path(base.dir,'CLR',FILES[[i]])))
fcs.data <- transform.scatter(applyTransforms(fcs.data,transforms))
csx <- apply.splits(as.data.frame(fcs.data[which(!as.logical(CLR[,'Lymphocytes'])),MARKERS]),cs$splits)
CLR <- cbind(CLR,RPART.1=0)
CLR[which(!as.logical(CLR[,'Lymphocytes'])),'RPART.1'] <- as.numeric(csx$X$cluster %in% k)
PSTAT5 <- cbind(PSTAT5,sapply( c('Lymphocytes','RPART.1'), function(cell.type) median(fcs.data[as.logical(CLR[,cell.type]),'PSTAT5']) ) )
} 
ylim <- range(sapply(c('Lymphocytes','RPART.1'), function(cell.type) PSTAT5[cell.type,]))
par(mfrow=c(1,1))
plot(NULL, xlim=c(1,4), ylim=ylim, xaxt='n', xlab='dose', ylab='pSTAT5 MFI', main='dose response')
axis(1, at=1:4, labels=DOSES)
i <- 1
for (cell.type in c('Lymphocytes','RPART.1')) {
    lines(1:4,  PSTAT5[cell.type,] , col=c('black','purple')[[i]], lwd=3)
    i <- i+1
} 
legend('topleft', c('Lymphocytes','RPART.1'), fill=c('black','purple'), bty='n')
dev.off()


#################### newcells
print(load(file.path(base.dir,'RData',FILES[[1]])))
print(load(file.path(base.dir,'CLR',FILES[[1]])))
fcs.data <- transform.scatter(applyTransforms(fcs.data,transforms))
csx <- apply.splits(as.data.frame(fcs.data[which(!as.logical(CLR[,'Lymphocytes'])),MARKERS]),cs$splits)
CLR <- cbind(CLR,RPART.1=0)
CLR[which(!as.logical(CLR[,'Lymphocytes'])),'RPART.1'] <- as.numeric(csx$X$cluster %in% k) 
save(CLR, file=file.path(base.dir,'newcells',sprintf('RPART-nonlymphocytes-%s',FILES[[1]])))


