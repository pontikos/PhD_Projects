################################################## Non-lymphocytes
print(load(file.path(base.dir,'RData','pstat5-join','CB00086S_2012-09-18.RData')))
print(load(file.path(base.dir,'CLR',FILES[[1]])))
fcs.data <- baseline.relative.pstat5(fcs.data)
X <- split.response(fcs.data, 4)

dim(X$d2$d2$d2$d) 
print(load(file.path(base.dir,'CLR',FILES[[1]])))
CLR <- cbind(CLR,MMPART.1=0) 
CLR <- cbind(CLR,MMPART.2=0) 
print(dim(y <- X$d2$d2$d2$d))
#rownames here are from the original fcs.data (not the subset)
CLR[as.numeric(rownames(y)),'MMPART.1'] <- 1
print(dim(y <- X$d2$d2$d1$d))
CLR[as.numeric(rownames(y)),'MMPART.2'] <- 1
#now remove all those which have already been included
CLR[which(CLR[,'Lymphocytes']>1),grep('MMPART',colnames(CLR))]<-0
colSums(CLR)

pdf('~/Thesis/figures/mmpart-nonlymphocytes-tree.pdf')
plot.tree(X,fcs.data,CLR,cells=c('MMPART.1','MMPART.2'))
dev.off()



#################### univariate cluster plots
plotNonLymphocytesClusters(file.name='~/Thesis/figures/mmpart-nonlymphocytes-clusters.pdf',CLR=CLR,cells=c('Lymphocytes','MMPART.1','MMPART.2'),fcs.data=fcs.data)


#################### dose response
plotNonLymphocytesDoseResponse(file.name='~/Thesis/figures/mmpart-nonlymphocytes-dose-response.pdf',CLR=CLR,cells=c('Lymphocytes','MMPART.1','MMPART.2'),fcs.data=fcs.data)


#################### newcells
save(CLR, file=file.path(base.dir,'newcells',sprintf('MMPART-nonlymphocytes-%s',FILES[[1]])))



