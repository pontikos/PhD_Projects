#### All lymphocytes
print(f <- file.path(base.dir,'RData','pstat5-join', sprintf('%s_%s.RData',individual,date)))
print(load(f))
fcs.data <- baseline.relative.pstat5(fcs.data)
#fcs.data <- applyTransforms(fcs.data,transforms)
print(load(file.path(BASE.DIR,'CLR',FILES[[1]])))
print(dim(fcs.data <- fcs.data[as.logical(CLR[,'Single cells']),]))
MARKERS <- CORE.FMARKERS

#
fm <- as.formula(paste(sprintf('diff.PSTAT5.%d',4), paste(MARKERS, collapse='+'), sep='~'))
plsr(fm, data=as.data.frame(fcs.data),scale=FALSE)->p
Y <- p$scores[,1:2]
smoothPlot(Y, outliers=TRUE)
g <- structure(list(x = c(0.549386221184538, 0.299526136195612, 0.103689853366454, 0.0834309275565415, 0.21173745768599, 0.468350517944886, 0.542633245914567), y = c(-1.4566307230717, -1.58756893468617, -1.92800828488379, -2.16369706578983, -2.22916617159706, -1.92800828488379, -1.4566307230717)), .Names = c("x", "y"))
lines(g$x,g$y,col='lightblue')
j <- as.logical(point.in.polygon(Y[,1],Y[,2],g$x,g$y))
g <- structure(list(x = c(-0.444233868802116, -0.92602314405666, -1.55234920188757, -1.27933527924333, -0.861784574022721, -0.492412796327571), y = c(-0.163647369673267, 0.787426504230537, -0.269322244551468, -1.07949628528434, -1.1499462018698, -0.128422411380534)), .Names = c("x", "y"))
lines(g$x,g$y,col='pink')
k <- as.logical(point.in.polygon(Y[,1],Y[,2],g$x,g$y))
g <- structure(list(x = c(-0.701188148937873, -0.379995298768177, 0.439046469164547, 0.053615048960912, -0.701188148937873), y = c(0.857876420816004, 0.0829273383758669, 1.03400121227967, 1.70327541984161, 0.998776253986938)), .Names = c("x", "y"))
lines(g$x,g$y,col='purple')
l <- as.logical(point.in.polygon(Y[,1],Y[,2],g$x,g$y))

#identification of pink and blue high density clusters
print(load(file.path(BASE.DIR,'CLR',FILES[[1]])))
print(dim(CLR <- CLR[as.logical(CLR[,'Single cells']),]))
CLR <- cbind(CLR,PLSR.1=as.numeric(l))
CLR <- cbind(CLR,PLSR.2=as.numeric(k))
CLR <- cbind(CLR,PLSR.3=as.numeric(j))

pdf('~/Thesis/figures/plsr-lymphocytes.pdf')
figure.labels <- iter(paste(letters,')',sep=''))
par(mfrow=c(2,2),mai=c(.75,.75,.5,.1))
cols <- c(CELL.TYPES.COL,'purple','pink','lightblue')
for (i in 2:4) {
    fm <- as.formula(paste(sprintf('diff.PSTAT5.%d',i), paste(MARKERS, collapse='+'), sep='~'))
    plsr(fm, data=as.data.frame(fcs.data),scale=FALSE)->p
    Y <- p$scores[,1:2]
    xlab <- 'Comp 1'
    #xlab <- ''
    ylab <- 'Comp 2'
    #ylab <- ''
    smoothPlot(Y,posteriors=CLR[,c(CELL.TYPES,'PLSR.1','PLSR.2','PLSR.3')],chulls=FALSE,clusters.col=cols, outliers=TRUE,ellipse.lwd=3,ylab=ylab,xlab=xlab)
    xlab <- paste(round(p$coefficients[,,1],2),'*',rownames(p$coefficients),collapse=' + ',sep='')
    ylab <- paste(round(p$coefficients[,,2],2),'*',rownames(p$coefficients),collapse=' + ',sep='')
    print(p$loadings)
    print(p$coefficients[,,1])
    print(p$coefficients[,,2])
    #plot(NULL,xlim=range(p$loadings[,1]),ylim=range(p$loadings[,2]))
    #for (m in MARKERS) {
        #segments(0,0,p$loadings[m,1],p$loadings[m,2])
        #text(p$loadings[m,1],p$loadings[m,2], label=m)
    #}
    #segments(0,0,p$loadings[m,1],p$loadings[m,2])
    #segments(0,0,p$Yloadings[,1],p$Yloadings[,2])
    #text(p$Yloadings[,1],p$Yloadings[,2],label=m)
    title(paste(nextElem(figure.labels),DOSES[[i]], sep='\t'), adj=0)
}
plot.new()
legend('center',legend=c(CELL.TYPES,paste('PLSR',1:3,sep='.')),fill=cols,cex=2,bty='n')
dev.off()

#################### univariate clusters
plotLymphocytesClusters(file.name='~/Thesis/figures/plsr-lymphocytes-clusters.pdf', CLR=CLR, cells=c(CELL.TYPES,'PLSR.1'), fcs.data=fcs.data)


#################### dose response
plotLymphocytesDoseResponse(file.name='~/Thesis/figures/plsr-lymphocytes-dose-response.pdf', CLR=CLR, cells=c(CELL.TYPES,'PLSR.1','PLSR.2','PLSR.3'), fcs.data=fcs.data)

#################### newcells
# only keep the purple population not the other two
CLR <- CLR[,c(!colnames(CLR) %in% c('PLSR.2','PLSR.3'))] 
#now remove all those which have already been included
CLR[which(rowSums(CLR[,CELL.TYPES])>0),grep('PLSR',colnames(CLR))]<-0
save(CLR, file=file.path(base.dir,'newcells',sprintf('PLSR-lymphocytes-%s',FILES[[1]])))


