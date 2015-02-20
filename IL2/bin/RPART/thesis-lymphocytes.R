################################################## Lymphocytes
## Lymphocytes/single cells which are not CD4+

#core fluorescent markers (no scatter)
MARKERS <- CORE.FMARKERS
ALL.CELL.TYPES <- CLR.CELL.TYPES
pool.fcs.data <- c()
for (i in 1:4) {
    print(load(file.path(base.dir,'RData',FILES[[i]])))
    print(load(file.path(base.dir,'CLR',FILES[[i]])))
    #only keep cells which are CD4+ lymphocytes but not memory, naive etc
    print(dim(fcs.data <- fcs.data[which(apply(CLR[,ALL.CELL.TYPES],1,function(x) sum(2**(0:8)[as.logical(x)]))==3),]))
    trans.fcs.data <- applyTransforms(fcs.data,transforms)
    pool.fcs.data <- rbind(pool.fcs.data, trans.fcs.data)
}
#create partitioning from pooled data
cs <- create.splits(pool.fcs.data[,MARKERS], branch.length=10) 
clusters <- c()
for (i in 1:4) {
    print(load(file.path(base.dir,'RData',FILES[[i]])))
    print(load(file.path(base.dir,'CLR',FILES[[i]])))
    print(dim(fcs.data <- fcs.data[which(apply(CLR[,ALL.CELL.TYPES],1,function(x) sum(2**(0:8)[as.logical(x)]))==3),]))
    trans.fcs.data <- applyTransforms(fcs.data,transforms)
    clusters <- rbind(clusters,cbind(cluster.stats(trans.fcs.data[,c(MARKERS,'PSTAT5')],cs),dose=DOSES[[i]]))
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
pdf('~/Thesis/figures/rpart-lymphocytes-mst-1024bin.pdf')
layout(matrix(c(1,2,3,4,5,5),nrow=3,ncol=2,byrow=T), heights=c(5,5,1),widths=c(5,5))
par(mai=c(.5,.5,.5,.5))
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
par(mai=c(.4,2,.15,2))
image( grad, 1, matrix(1:length(colorscale),ncol=1), col=colorscale, xlab='', ylab="", yaxt="n", main='pSTAT5 MFI')
dev.off()

k1 <- as.numeric(rownames(X)[which(as.logical(point.in.polygon(mst.layout[,1],mst.layout[,2],G1$x,G1$y)))])
k2 <- as.numeric(rownames(X)[which(as.logical(point.in.polygon(mst.layout[,1],mst.layout[,2],G2$x,G2$y)))])


#################### univariate clusters
print(load(file.path(base.dir,'RData',FILES[[1]])))
print(load(file.path(base.dir,'CLR',FILES[[1]])))
dim(fcs.data <- applyTransforms(fcs.data[which(as.logical(CLR[,'Single cells'])),],transforms))
dim(CLR <- CLR[which(as.logical(CLR[,'Single cells'])),])
csx <- apply.splits(as.data.frame(fcs.data[,MARKERS]),cs$splits)
CLR <- cbind(CLR,RPART.1=0)
CLR[,'RPART.1'] <- as.numeric(csx$X$cluster %in% k2 & (rowSums(CLR[,CELL.TYPES])==0))
CLR <- cbind(CLR,RPART.2=0)
CLR[,'RPART.2'] <- as.numeric(csx$X$cluster %in% k1 & (rowSums(CLR[,CELL.TYPES])==0)) 

plotLymphocytesClusters(file.name='~/Thesis/figures/rpart-lymphocytes-clusters.pdf',CLR=CLR,cells=c(CELL.TYPES,'RPART.1','RPART.2'),fcs.data=fcs.data)

#################### dose response
pdf('~/Thesis/figures/rpart-lymphocytes-dose-response.pdf')
PSTAT5 <- c()
for (i in 1:4) {
print(load(file.path(base.dir,'RData',FILES[[i]])))
print(load(file.path(base.dir,'CLR',FILES[[i]])))
fcs.data <- applyTransforms(fcs.data[which(as.logical(CLR[,'Single cells'])),],transforms)
CLR <- CLR[which(as.logical(CLR[,'Single cells'])),]
csx <- apply.splits(as.data.frame(fcs.data[,MARKERS]),cs$splits)
CLR <- cbind(CLR,RPART.1=0)
CLR[,'RPART.1'] <- as.numeric(csx$X$cluster %in% k2 & (rowSums(CLR[,CELL.TYPES])==0))
CLR <- cbind(CLR,RPART.2=0)
CLR[,'RPART.2'] <- as.numeric(csx$X$cluster %in% k1 & (rowSums(CLR[,CELL.TYPES])==0))
PSTAT5 <- cbind(PSTAT5,sapply( c(CELL.TYPES,'RPART.1','RPART.2'), function(cell.type) median(fcs.data[as.logical(CLR[,cell.type]),'PSTAT5']) ) )
} 
ylim <- range(sapply(c(CELL.TYPES,'RPART.1','RPART.2'), function(cell.type) PSTAT5[cell.type,]))
cols <- c(CELL.TYPES.COL,'purple','pink')
par(mfrow=c(1,1))
plot(NULL, xlim=c(1,4), ylim=ylim, xaxt='n', xlab='dose', ylab='pSTAT5 MFI', main='dose response')
axis(1, at=1:4, labels=DOSES)
i <- 1
for (cell.type in c(CELL.TYPES,'RPART.1','RPART.2')) {
    lines(1:4,  PSTAT5[cell.type,] , col=cols[[i]], lwd=3)
    i <- i+1
} 
legend('topleft', c(CELL.TYPES,'RPART.1','RPART.2'), fill=cols, bty='n')
dev.off()


#################### newcells
print(load(file.path(base.dir,'RData',FILES[[1]])))
print(load(file.path(base.dir,'CLR',FILES[[1]])))
dim(fcs.data <- applyTransforms(fcs.data[which(as.logical(CLR[,'Single cells'])),],transforms))
dim(CLR <- CLR[which(as.logical(CLR[,'Single cells'])),])
csx <- apply.splits(as.data.frame(fcs.data[,MARKERS]),cs$splits)
CLR <- cbind(CLR,RPART.1=0)
CLR[,'RPART.1'] <- as.numeric(csx$X$cluster %in% k2 & (rowSums(CLR[,CELL.TYPES])==0))
CLR <- cbind(CLR,RPART.2=0)
CLR[,'RPART.2'] <- as.numeric(csx$X$cluster %in% k1 & (rowSums(CLR[,CELL.TYPES])==0))
#now remove all those which have already been included
CLR[which(rowSums(CLR[,CELL.TYPES])>0),grep('RPART',colnames(CLR))]<-0 
save(CLR, file=file.path(base.dir,'newcells',sprintf('RPART-lymphocytes-%s',FILES[[1]])))


