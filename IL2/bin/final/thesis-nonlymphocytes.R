
METHODS <- c('SPADE','RPART','PLSR','MMPART')
  
for (m in METHODS) {
load(file.path(base.dir,'newcells',sprintf('%s-nonlymphocytes-%s',m,FILES[[1]])))
print(dim(CLR))
}


newcells <- do.call('cbind', lapply(c('SPADE','RPART','PLSR','MMPART'), function(m) {
    f <- file.path(base.dir,'newcells',sprintf('%s-nonlymphocytes-%s',m,FILES[[1]]))
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
#fcs.data <- baseline.relative.pstat5(applyTransforms(fcs.data[single.cells,],transforms)
fcs.data <- baseline.relative.pstat5(fcs.data)

#remove cells which also belong to the lymphocytes cluster
CLR.newcells[which(as.logical(CLR[,'Lymphocytes'])),] <- 0

clr.consensus <- cbind(CLR[,CLR.CELL.TYPES],as.numeric(rowSums(CLR.newcells)>1))

################### write out summary table
summary.table <- round( t(apply(cbind(Lymphocytes=CLR[,'Lymphocytes'],CLR.newcells,consensus=as.numeric(rowSums(CLR.newcells)>1)), 2, function(x) c(colMedians(fcs.data[as.logical(x),]),freq=100*sum(x)/nrow(CLR)))), digit=2 )
write.csv( summary.table[,-grep('PSTAT5',colnames(summary.table))], quote=FALSE, file='' )

#################### univariate clusters
plotNonLymphocytesClusters(file.name='~/Thesis/figures/final-nonlymphocytes-clusters.pdf',CLR=CLR.newcells,cells=colnames(CLR.newcells),fcs.data=fcs.data)

# cell to more than one group
pdf('~/Thesis/figures/final-nonlymphocytes-clusters-consensus.pdf')
par(mfrow=c(3,3))
for (marker in CORE.MARKERS)
smoothPlot1D(fcs.data[,marker], outliers=TRUE, main=marker, clusters.col=c(CELL.TYPES.COL,'purple'), posteriors=clr.consensus, col='white')
dev.off()


#################### dose response
plotNonLymphocytesDoseResponse(file.name='~/Thesis/figures/final-nonlymphocytes-dose-response.pdf',CLR=CLR.newcells,cell=colnames(CLR.newcells),fcs.data=fcs.data)


