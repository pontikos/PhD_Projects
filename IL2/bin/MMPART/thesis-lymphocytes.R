################################################## Lymphocytes
print(load(file.path(base.dir,'RData','pstat5-join','CB00086S_2012-09-18.RData')))
print(load(file.path(base.dir,'CLR',FILES[[1]])))
# data is already transformed
#print(load('~/dunwich/Projects/IL2/transforms.RData'))
#fcs.data <- as.data.frame(baseline.relative.pstat5(applyTransforms(fcs.data,transforms)))
fcs.data <- as.data.frame(baseline.relative.pstat5(fcs.data))
fcs.data <- fcs.data[as.logical(CLR[,'Single cells']),]
X <- split.response(fcs.data, 4)

print(load(file.path(base.dir,'CLR',FILES[[1]])))
CLR <- cbind(CLR,MMPART.1=0) 
CLR <- cbind(CLR,MMPART.2=0) 
CLR <- cbind(CLR,MMPART.3=0) 
#rownames here are from the original fcs.data (not the subset)
# high, high, high
print(dim(y <- X$d2$d2$d2$d))
CLR[as.numeric(rownames(y)),'MMPART.1'] <- 1
# low, high, high
print(dim(y <- X$d2$d2$d1$d))
CLR[as.numeric(rownames(y)),'MMPART.2'] <- 1
# low, low, high
print(dim(y <- X$d2$d1$d1$d))
CLR[as.numeric(rownames(y)),'MMPART.3'] <- 1
#now remove all those which have already been included
CLR[which(rowSums(CLR[,CELL.TYPES])>0),grep('MMPART',colnames(CLR))]<-0
CLR <- CLR[as.logical(CLR[,'Single cells']),]
colSums(CLR)

#
pdf('~/Thesis/figures/mmpart-lymphocytes-tree.pdf')
plot.tree(X,fcs.data,CLR,cells=c('MMPART.1','MMPART.2','MMPART.3'))
dev.off()


#################### univariate cluster plots
plotLymphocytesClusters(fcs.data=fcs.data, file.name='~/Thesis/figures/mmpart-lymphocytes-clusters.pdf', CLR=CLR, cells=c(CELL.TYPES,'MMPART.1','MMPART.2','MMPART.3'))

#################### dose response
plotLymphocytesDoseResponse(fcs.data=fcs.data, file.name='~/Thesis/figures/mmpart-lymphocytes-dose-response.pdf', CLR=CLR, cells=c(CELL.TYPES,'MMPART.1','MMPART.2','MMPART.3'))

#################### newcells
save(CLR, file=file.path(base.dir,'newcells',sprintf('MMPART-lymphocytes-%s',FILES[[1]])))



