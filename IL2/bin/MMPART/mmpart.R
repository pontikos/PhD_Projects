source('~nikolas/Projects/IL2/bin/MMPART/mmpart-functions.R',chdir=T)

MARKERS <- CORE.FMARKERS

################################################## Lymphocytes
print(load('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/RData/pstat5-join/CB00086S_2012-09-18.RData'))
print(load('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/CLR/CB00086S_2012-09-18.RData'))
# data is already transformed
#print(load('~/dunwich/Projects/IL2/transforms.RData'))
#fcs.data <- as.data.frame(baseline.relative.pstat5(applyTransforms(fcs.data,transforms)))
fcs.data <- as.data.frame(baseline.relative.pstat5(fcs.data))
fcs.data <- fcs.data[as.logical(CLR[,'Single cells']),]
X <- split.response(fcs.data, 4)

#
pdf('~/Thesis/figures/mmpart-lymphocytes-tree.pdf')
plot.tree(X)
dev.off()

dim(X$d2$d1$d) 
print(load('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/CLR/CB00086S_2012-09-18.RData'))
CLR <- cbind(CLR,pink=0) 
CLR <- cbind(CLR,purple=0) 
#print(dim(y <- X$d1$d1$d2$d))
#pink
print(dim(y <- X$d1$d2$d))
#rownames here are from the original fcs.data (not the subset)
CLR[as.numeric(rownames(y)),'pink'] <- 1
#purple
print(dim(y <- X$d2$d2$d))
CLR[as.numeric(rownames(y)),'purple'] <- 1
#plot.points=y, plot.points.col='purple',
CLR <- CLR[as.logical(CLR[,'Single cells']),]
colSums(CLR)

#################### univariate cluster plots
pdf('~/Thesis/figures/mmpart-lymphocytes-clusters.pdf')
par(mfrow=c(3,3))
for (marker in CORE.MARKERS)
smoothPlot1D(fcs.data[,marker], posteriors=CLR[as.logical(CLR[,'Single cells']),c(CELL.TYPES,'purple','pink')],outliers=TRUE, clusters.col=c(CELL.TYPES.COL,'purple','pink'), main=marker)
dev.off()

#################### dose response
pdf('~/Thesis/figures/mmpart-lymphocytes-dose-response.pdf')
#dose response
ylim <- range( sapply(c(CELL.TYPES,'pink','purple'), function(cell.type) colMedians(fcs.data[as.logical(CLR[,cell.type]),paste('diff.PSTAT5',1:4,sep='.')])) )
plot(NULL, xlim=c(0,3), ylim=ylim, xaxt='n', xlab='dose', ylab='pSTAT5 MFI', main='dose response')
#title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
cols <- c(CELL.TYPES.COL,'pink','purple')
for (cell.type in c(CELL.TYPES, 'pink', 'purple')) {
    mfi <- fcs.data[which(as.logical(CLR[,cell.type])),paste('diff.PSTAT5',1:4,sep='.')]
    lines(0:3, colMedians( mfi ), col=cols[[i]], lwd=3)
    i <- i+1
}
dev.off()

#plotClusters(fcs.data[,MARKERS], posteriors=CLR[as.logical(CLR[,'Single cells']),c('Memory Eff','Memory Treg','Naive Eff','Naive Treg','purple')],chulls=FALSE,outliers=TRUE, clusters.col=c('black','red','green','blue','purple'))


################################################## Non-lymphocytes
MARKERS <- CORE.FMARKERS
print(load('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/RData/pstat5-join/CB00086S_2012-09-18.RData'))
print(load('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/CLR/CB00086S_2012-09-18.RData'))
fcs.data <- baseline.relative.pstat5(fcs.data)
X <- split.response(fcs.data, 4)

pdf('~/Thesis/figures/mmpart-nonlymphocytes-tree.pdf')
plot.tree(X)
dev.off()

dim(X$d2$d2$d2$d) 
print(load('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/CLR/CB00086S_2012-09-18.RData'))
CLR <- cbind(CLR,pink=0) 
CLR <- cbind(CLR,purple=0) 
#print(dim(y <- X$d1$d1$d2$d))
#pink
#print(dim(y <- X$d1$d2$d))
print(dim(y <- X$d2$d2$d2$d))
#rownames here are from the original fcs.data (not the subset)
CLR[as.numeric(rownames(y)),'pink'] <- 1
#purple
print(dim(y <- X$d2$d2$d))
CLR[as.numeric(rownames(y)),'purple'] <- 1
#plot.points=y, plot.points.col='purple',
colSums(CLR)


pdf('~/Thesis/figures/mmpart-nonlymphocytes-clusters.pdf')
par(mfrow=c(3,3))
for (marker in CORE.MARKERS)
smoothPlot1D(fcs.data[,marker], posteriors=CLR[,c('Lymphocytes','purple','pink')],outliers=TRUE, clusters.col=c('black','purple','pink'), main=marker)
dev.off()

pdf('~/Thesis/figures/mmpart-nonlymphocytes-dose-response.pdf')
par(mfrow=c(1,1))
#dose response
print( ylim <- range( sapply(c('Lymphocytes','pink','purple'), function(cell.type) colMedians(fcs.data[as.logical(CLR[,cell.type]),paste('diff.PSTAT5',1:4,sep='.')]))) )
plot(NULL, xlim=c(0,3), ylim=ylim, xaxt='n', xlab='dose', ylab='pSTAT5 MFI', main='dose response')
#title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
cols <- c('black','pink','purple')
for (cell.type in c('Lymphocytes', 'pink', 'purple')) {
    mfi <- fcs.data[which(as.logical(CLR[,cell.type])),paste('diff.PSTAT5',1:4,sep='.')]
    lines(0:3, colMedians( mfi ), col=cols[[i]], lwd=3)
    i <- i+1
}
dev.off()



y <- X$d2$d2$d2$d
MARKERS <- c('SSCA','FSCA','CD25','CD3','CD4','CD45RA','CD56','CD8','FOXP3')
plotClusters(fcs.data[,c('FSCA','SSCA')], plot.points=X$d2$d2$d2$d, plot.points.col='purple', plot.points.pch='.', chulls=FALSE,outliers=TRUE) 
i <- which(apply(t(apply(y[,paste('diff','PSTAT5',1:4,sep='.')], 1, diff)),1,function(x) all(x>0))) 
print(load('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/RData/CB00086S_01U_2012-09-18.RData'))
print(load('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/RData/CB00086S_10U_2012-09-18.RData'))
print(load('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/RData/CB00086S_1000U_2012-09-18.RData'))
plotClusters(fcs.data[,c('FSCA','SSCA')], plot.points=y[i,], plot.points.col='purple', plot.points.pch=20, chulls=FALSE,outliers=TRUE)

