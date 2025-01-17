METHODS <- c('SPADE','RPART','PLSR')

for (m in METHODS) {
load(file.path(base.dir,'newcells',sprintf('%s-lymphocytes-%s',m,FILES[[1]])))
print(dim(CLR))
}

newcells <- do.call('cbind', lapply(METHODS, function(m) {
    f <- file.path(base.dir,'newcells',sprintf('%s-lymphocytes-%s',m,FILES[[1]]))
    if (file.exists(f)) {
        print(load(f))
        i <- which(!colnames(CLR) %in%  CLR.CELL.TYPES)
        colnames(CLR)[i] <- paste(m,1:length(i),sep='.')
        CLR
    }
}))

print(head(newcells <- newcells[,c(!colnames(newcells) %in% CLR.CELL.TYPES)]))
CLR.newcells <- newcells

load(file.path(base.dir,'RData','pstat5-join',sprintf('%s_%s.RData',individual,date)))
load(file.path(base.dir,'CLR',FILES[[1]]))
single.cells <- as.logical(CLR[,'Single cells'])
#fcs.data <- baseline.relative.pstat5(applyTransforms(fcs.data[single.cells,],transforms)
fcs.data <- baseline.relative.pstat5(fcs.data[single.cells,])
CLR <- CLR[single.cells,]


#remove cells which also belong to the known subsets
CLR.newcells[which(rowSums(CLR[,CELL.TYPES])>0),] <- 0

CLR.consensus <- clr.consensus <- cbind(CLR[,CELL.TYPES],consensus=as.numeric(rowSums(CLR.newcells)>1))

################### write out summary table
summary.table <- round( t(apply(cbind(CLR[,CELL.TYPES],CLR.newcells,consensus=as.numeric(rowSums(CLR.newcells)>1)), 2, function(x) c(colMedians(fcs.data[as.logical(x),]),freq=100*sum(x)/nrow(CLR)))), digit=2 )
write.csv( summary.table[,-grep('PSTAT5',colnames(summary.table))], quote=FALSE, file='' )

#################### univariate clusters
plotLymphocytesClusters(file.name='~/Thesis/figures/final-lymphocytes-clusters.pdf', CLR=cbind('Single cells'=1, CLR.newcells), cells=colnames(CLR.newcells), fcs.data=fcs.data, cols=NEW.CELL.TYPES.COL)

# cell to more than one group
pdf('~/Thesis/figures/final-lymphocytes-clusters-consensus.pdf',width=10,height=10)
par(mfrow=c(3,3))
for (marker in CORE.MARKERS)
smoothPlot1D(fcs.data[,marker], outliers=TRUE, main=marker, clusters.col=c(CELL.TYPES.COL,'purple'), posteriors=clr.consensus, col='white')
dev.off()
 
#################### dose response
plotLymphocytesDoseResponse(file.name='~/Thesis/figures/final-lymphocytes-dose-response.pdf', CLR=CLR.newcells, cells=colnames(CLR.newcells), fcs.data=fcs.data, cols=NEW.CELL.TYPES.COL) 

pdf('~/Thesis/figures/final-lymphocytes-dose-response-consensus.pdf',width=5,height=5)
ylim <- range(sapply(colnames(CLR.consensus), function(cell.type) colMedians(fcs.data[as.logical(CLR.consensus[,cell.type]),paste('diff.PSTAT5',1:4,sep='.')]) ))
par(mfrow=c(1,1))
plot(NULL, xlim=c(0,3), ylim=ylim, xaxt='n', xlab='dose', ylab='pSTAT5 MFI')
axis(1, at=0:3, labels=DOSES)
i <- 1
cols <- c(CELL.TYPES.COL,'purple')
for (cell.type in colnames(CLR.consensus)) {
    mfi <- fcs.data[which(as.logical(CLR.consensus[,cell.type])),paste('diff.PSTAT5',1:4,sep='.')]
    lines(0:3, colMedians( mfi ), col=cols[[i]], lwd=3)
    #lines(0:3, colMedians( mfi ), lwd=3)
    i <- i+1
} 
dev.off()


