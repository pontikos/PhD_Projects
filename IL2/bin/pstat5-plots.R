### Other scatter clusters 

### 

for (k in 1:4) {
setwd(file.path(base.dir,'SPADE',paste(individual,date,sep='_'),'scatter',k))
print(nclust <- nrow(clusters.table <- read.table('clusters.table',header=TRUE)))
graph <- read.graph('mst.gml',format='gml')
layout.table <- SPADE.layout.arch(graph)
pdf(sprintf('~nikolas/Thesis/figures/spade-%d-pstat5mfi.pdf',k),width=10,height=10)
plot.mst.fold(files, layout.table=layout.table, boundaries=boundaries)
dev.off()
}


#color <- rgb(0,0,0,alpha=clust.cell.type[,'Memory Eff'])
#plot(graph, layout=layout.table, vertex.shape="circle", vertex.color=color, vertex.frame.color=color, edge.color='grey', vertex.size=4, vertex.label=NA, edge.arrow.size=0.25, edge.arrow.width=1)
#color <- rgb(1,0,0,alpha=clust.cell.type[,'Memory Treg']) 
#par(new=T)
#plot(graph, layout=layout.table, vertex.shape="circle", vertex.color=color, vertex.frame.color=color, edge.color='grey', vertex.size=4, vertex.label=NA, edge.arrow.size=0.25, edge.arrow.width=1)
#color <- rgb(0,0.5,0,alpha=clust.cell.type[,'Naive Eff'])
#par(new=T)
#plot(graph, layout=layout.table, vertex.shape="circle", vertex.color=color, vertex.frame.color=color, edge.color='grey', vertex.size=4, vertex.label=NA, edge.arrow.size=0.25, edge.arrow.width=1)
#color <- rgb(0,0,1,alpha=clust.cell.type[,'Naive Treg'])
#par(new=T)
#plot(graph, layout=layout.table, vertex.shape="circle", vertex.color=color, vertex.frame.color=color, edge.color='grey', vertex.size=4, vertex.label=NA, edge.arrow.size=0.25, edge.arrow.width=1)


#
setwd(file.path(base.dir,'SPADE',paste(individual,date,sep='_'),'Not.Lymphocytes'))
#setwd(file.path(base.dir,'SPADE',paste(individual,date,sep='_'),'scatter',))
#setwd(file.path(base.dir,'SPADE',paste(individual,date,sep='_')))
#
print(nclust <- nrow(clusters.table <- read.table('clusters.table',header=TRUE)))
graph <- read.graph('mst.gml',format='gml')
layout.table <- SPADE.layout.arch(graph)
#
par(mfrow=c(2,2))
figure.labels <- iter(paste(letters,')',sep=''))
for ( f in sprintf("%s.RData",paste(individual,c('0U','01U','10U','1000U'),date,sep='_')) ) {
clust.cell.type <- data.frame(matrix(0,nrow=nclust,ncol=length(c('unknown',CELL.TYPES))))
colnames(clust.cell.type) <- c('unknown',CELL.TYPES) 
load(f)
load(file.path(base.dir,'CLR',f))
CLR <- CLR[as.logical(CLR[,'Single cells']),]
cell.types <- as.character(apply(CLR[,CELL.TYPES],1,function(x) names(which(x==1))))
#cell.types[!cell.types %in% CELL.TYPES] <- 'unknown'
dim(temp <- table(fcs.data[,'cluster'],cell.types))
colnames(temp) <- c('unknown',CELL.TYPES)
clust.cell.type[rownames(temp),] <- temp
clust.cell.type <- t(apply(clust.cell.type,1,function(x)x/sum(x)))
clust.cell.type[is.nan(clust.cell.type)] <- 0
colnames(clust.cell.type) <- c('unknown',CELL.TYPES)
#unknown <- rownames(clust.cell.type[which(clust.cell.type[,'unknown']==1),])
#length(pstat5.4 <- tapply(fcs.data$PSTAT5, fcs.data$cluster, median))
if (grepl('_0U_',f)) {
#color <- rgb(0,0,0,alpha=clust.cell.type[,'Memory Eff'])
#plot(graph, layout=layout.table, vertex.shape="circle", vertex.color=color, vertex.frame.color=color, edge.color='grey', vertex.size=4, vertex.label=NA, edge.arrow.size=0.25, edge.arrow.width=1)
#color <- rgb(1,0,0,alpha=clust.cell.type[,'Memory Treg']) 
#par(new=T)
#plot(graph, layout=layout.table, vertex.shape="circle", vertex.color=color, vertex.frame.color=color, edge.color='grey', vertex.size=4, vertex.label=NA, edge.arrow.size=0.25, edge.arrow.width=1)
#color <- rgb(0,0.5,0,alpha=clust.cell.type[,'Naive Eff'])
#par(new=T)
#plot(graph, layout=layout.table, vertex.shape="circle", vertex.color=color, vertex.frame.color=color, edge.color='grey', vertex.size=4, vertex.label=NA, edge.arrow.size=0.25, edge.arrow.width=1)
#color <- rgb(0,0,1,alpha=clust.cell.type[,'Naive Treg'])
#par(new=T)
#plot(graph, layout=layout.table, vertex.shape="circle", vertex.color=color, vertex.frame.color=color, edge.color='grey', vertex.size=4, vertex.label=NA, edge.arrow.size=0.25, edge.arrow.width=1)
plot(layout.table, col='grey', pch=21, yaxt='n', xaxt='n', ann=FALSE,frame.plot=FALSE)
title(nextElem(figure.labels), adj=0)
points(layout.table, pch=20, col=rgb(0,0,0,alpha=clust.cell.type[,'Memory Eff']))
points(layout.table, pch=20, col=rgb(1,0,0,alpha=clust.cell.type[,'Memory Treg']))
points(layout.table, pch=20, col=rgb(0,0.5,0,alpha=clust.cell.type[,'Naive Eff']))
points(layout.table, pch=20, col=rgb(0,0,1,alpha=clust.cell.type[,'Naive Treg']))
} else {
plot.mst(fcs.data,param='PSTAT5',layout.table=layout.table)
title(nextElem(figure.labels), adj=0)
}
} 

pc <- princomp(clusters.table)
plot(pc$scores[,1:2])

x <- layout.table[,1]
x <- layout.table[,2]
x <- pc$scores[,1]
x <- pc$scores[,2]
plot(NULL,xlim=range(x),ylim=range(clusters.table))
for (i in 1:ncol(clusters.table)) {
    y <- clusters.table[,i]
    #points(x,y,col=i,pch=20)
    lines(lowess(x,y), col=i, lwd=2)
}
legend('topright',colnames(clusters.table),col=1:ncol(clusters.table),lwd=2)

purple <- numeric(nrow(CLR))
table(i <- ((fcs.data$cluster %in% k) & (rowSums(CLR[,CELL.TYPES])==0)))
purple[i] <- 1
CLR <- cbind(CLR, purple=purple)

plotClusters(fcs.data[as.logical(CLR[,'Single cells']),colnames(clusters.table)], plot.points=clusters.table[k,], outliers=TRUE, posteriors=CLR[,c(CELL.TYPES,'purple')], chulls=FALSE, plot.points.col=5, plot.points.pch='.')

plotClusters(fcs.data[as.logical(!CLR[,'Lymphocytes']),colnames(clusters.table)], plot.points=clusters.table[k,], outliers=TRUE, plot.points.col='purple', plot.points.pch=20)

plot.mst(fcs.data, param='CD25')

#n <- intersect(names(pstat5.1),names(pstat5.4))
n <- intersect(unknown.1, unknown.2)

y <- pstat5.4[n]-pstat5.1[n] 

par(mfrow=c(1,1))
plot(density(y))
res <- kmeans(y,2)
abline(v=max(res$centers))

print(resp.clusters <- names(which(res$cluster==which.max(res$centers))))
dim(resp.fcs.data <- fcs.data[fcs.data$cluster %in% resp.clusters,])
responsive <- clusters.table[resp.clusters,]
X <- fcs.data[,c('CD3','CD4','CD45RA','CD25','CD56','CD8')]
cl <- fcs.data[,'cluster']
table(cl[!cl %in% resp.clusters] <- 0)
plotClusters(X, outliers=TRUE, plot.points=resp.fcs.data, plot.points.col=as.factor(cl), classification=cl)

plotClusters(X, classification=as.numeric(clusters@exprs[,'cluster'] %in% k), outliers=TRUE, plot.points=X, plot.points.col='purple', chulls=FALSE)

plot( log(clusters.table[n,'CD25']), y )

summary( lm(y ~ ., data=log(clusters.table[n,])) )


### 2D binning
load('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/RData/pstat5-join/CB00406Q_2012-06-12.RData')
fcs.data <- applyTransforms(fcs.data, transforms)
fcs.data <- baseline.relative.pstat5(fcs.data)

pstat5 <- fcs.data[,'diff.PSTAT5.4']
#pstat5 <- fcs.data[,'PSTAT5']
b <- binning(x=as.matrix(fcs.data[,c('FSCA','SSCA')]),y=pstat5,nbins=c(512,512))

palette <- colorRampPalette(c("blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red"))
colorscale <- palette(100)
boundaries <- quantile(pstat5,probs=c(.02,.98))
grad <- seq(boundaries[[1]],boundaries[[2]],length.out=length(colorscale)) 

col <- colorscale[findInterval(b$means, grad, all.inside=TRUE)]
plot(b$x,col=col,pch=20,cex=.25)

ggplot(data=data.frame(b$x),aes(x=FSCA,y=SSCA))+geom_point(col=col,alpha=.25)

print(f <- sprintf('CB00406Q_%s_2012-06-12.RData', '0U'))
print(load(f))
b0 <- binning(x=as.matrix(fcs.data[,c('FSCA','SSCA')]),y=fcs.data[,'PSTAT5'],nbins=c(512,512))
par(mfrow=c(2,2))
for (dose in DOSES) {
    print(f <- sprintf('CB00406Q_%s_2012-06-12.RData', dose))
    print(load(f))
    #fcs.data <- applyTransforms(fcs.data, transforms)
    b <- binning(x=as.matrix(fcs.data[,c('FSCA','SSCA')]),y=fcs.data[,'PSTAT5'],nbins=c(512,512))
    col <- colorscale[findInterval(b$means-b0$means, grad, all.inside=TRUE)]
    plot(b$x,col=col,pch=20,cex=.25)
}







