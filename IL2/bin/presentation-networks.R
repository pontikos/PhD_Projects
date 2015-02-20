source('~nikolas/IL2/bin/functions-graph.R')
library(geometry)


lymph.graph.all <- create.graph(lymph.clusters.stats, all.channels)
lymph.graph.core <- create.graph(lymph.clusters.stats, core.channels)
lymph.graph.func <- create.graph(lymph.clusters.stats, func.channels)


lymph.sumdist.all <- lapply( lymph.clusters.stats, function(x) sum(dist(getChannels(x, all.channels))) )

lymph.ch.all <- lapply( lymph.clusters.stats, function(x) convhulln(getChannels(x, all.channels), 'FS') )
lymph.ch.vol.all <- lapply( lymph.ch.all, function(x) x$vol )
lymph.ch.area.all <- lapply( lymph.ch.all, function(x) x$area )

day.vol <- unflatten.day(lymph.ch.vol.all, doses=doses, f=function(x)apply(x,2,unlist))
day.area <- unflatten.day(lymph.ch.area.all, doses=doses, f=function(x)apply(x,2,unlist))


plot.doses.per.day(list(day1=day.vol$day1/day.area$day1, day2=day.vol$day2/day.area$day2), file.name='vol-area-dose-day', dir='~nikolas/IL2/Plots/lymphocytes/graph/convex-hull/')
plot.doses.per.day(day.vol, file.name='vol-dose-day', dir='~nikolas/IL2/Plots/lymphocytes/graph/convex-hull/')
plot.doses.per.day(day.area, file.name='area-dose-day', dir='~nikolas/IL2/Plots/lymphocytes/graph/convex-hull/')
plot.doses.per.day(unflatten.day(lymphocytes.sumdist.all, doses=doses, f=function(x)apply(x,2,unlist)), file.name='sumdist-dose-day', dir='~nikolas/IL2/Plots/lymphocytes/graph/sumdist/')

plot.doses.per.day.best.fit(day.vol, file.name='vol-dose-day', dir='~nikolas/IL2/Plots/lymphocytes/graph/convex-hull/')
plot.doses.per.day.best.fit(day.area, file.name='area-dose-day', dir='~nikolas/IL2/Plots/lymphocytes/graph/convex-hull/')
plot.doses.per.day.best.fit(unflatten.day(lymphocytes.sumdist.all, doses=doses, f=function(x)apply(x,2,unlist)), file.name='sumdist-dose-day-bestfit', dir='~nikolas/IL2/Plots/lymphocytes/graph/sumdist/')


plot.doses.per.day.best.fit(list(day1=day.vol$day1/day.area$day1, day2=day.vol$day2/day.area$day2), file.name='vol-area-dose-day-bestfit', dir='~nikolas/IL2/Plots/lymphocytes/graph/convex-hull/')






cor(day.vol[[1]]/day.area[[1]], day.vol[[2]]/day.area[[2]])**2

cor(dose.vol[[1]]/dose.area[[1]], dose.vol[[2]]/dose.area[[2]])**2

#convex hull
plot.agreement(dose0=dose.vol[[1]], dose1=dose.vol[[2]], day0=day.vol[[1]], day1=day.vol[[2]], 'volume', dir='~nikolas/IL2/Plots/lymph/graph/convex-hull/')
plot.agreement(dose0=dose.area[[1]], dose1=dose.area[[2]], day0=day.area[[1]], day1=day.area[[2]], 'area', dir='~nikolas/IL2/Plots/lymph/graph/convex-hull/')
plot.agreement(dose0=dose.vol[[1]]/dose.area[[1]], dose1=dose.vol[[2]]/dose.area[[1]], day0=day.vol[[1]]/day.area[[1]], day1=day.vol[[2]]/day.area[[2]], 'vol-area', dir='~nikolas/IL2/Plots/lymph/graph/convex-hull/')


#
doses <- unflatten.dose(lymph.clusters.stats, do.unlist=FALSE)
print(mapply(function(dose0,dose1) emdw(A=dose0[,2:6],B=dose1[,2:6],wA=dose0[,'pct'],wB=dose1[,'pct'],dist='manhattan'), doses$resting, doses$stimulated))
#
days <- unflatten.day(lymph.clusters.stats, do.unlist=FALSE)
print(mapply(function(day0,day1) emdw(A=day0[,2:6],B=day1[,2:6],wA=day0[,'pct'],wB=day1[,'pct'],dist='manhattan'), days$day1, days$day2))


desired_samples <- 1000
k <- 200
#k <- 10
#k <- 3

#
lymph.density <- lapply(lymph.fcs, multi.density, apprx_mult=2, channels=core.channels)
# downsample
lymph.down <- mapply(downsample.FCS, lymph.fcs, lymph.density, desired_samples=desired_samples)
#plot
marginal.2D.density(lymph.fcs[[s]], lymph.down[[s]], dir=sprintf('~nikolas/IL2/Plots/downsample-%d/%s/',desired_samples,s))
#
lymph.all.dist <- lapply(lymph.down, function(x) dist(getChannels(x, channels=all.channels)))
lymph.core.dist <- lapply(lymph.down, function(x) dist(getChannels(x, channels=core.channels)))
lymph.func.dist <- lapply(lymph.down, function(x) dist(getChannels(x, channels=func.channels)))

lymph.hclust <- lapply(lymph.dist, hclust)
# cluster
lymph.clust <- lapply(lymph.hclust,cutree,k=k)
lymph.clust <- lapply(lymph.down, function(x) pam(getChannels(x, core.channels), k=k)$clustering)
lymph.clust <- lapply(lymph.down, function(x) SPADE.cluster(x, k=k)$assign)
#lymph.clust.centers <- lapply(lymph.down, function(x) SPADE.cluster(x, k=k)$centers)
lymph.assign <- mapply(function(x,y,z) SPADE.assignToCluster(getChannels(x, core.channels),getChannels(y, core.channels),z), lymph.fcs, lymph.down, lymph.clust)
#lymph.fcs.clust <- mapply(function(x,y) build.flowFrame(cbind(getChannels(x,channels=channels),col=as.factor(y))), lymph.fcs, lymph.assign)
#compute stats for all markers core + functional
lymph.clusters.stats <- mapply(compute.cluster.stats, lymph.fcs, lymph.assign, MoreArgs=list(channels=all.channels), SIMPLIFY=FALSE, USE.NAMES=TRUE) 
lymph.clusters.down.stats <- mapply(compute.cluster.stats, lymph.down, lymph.clust, MoreArgs=list(channels=all.channels), SIMPLIFY=FALSE, USE.NAMES=TRUE) 

#
for (individual in rep.individuals) {
    for (dose in doses) {
        individual.dose <- paste(individual, dose, sep='.')
        print(individual.dose)
        print(s <- grep(individual.dose, names(lymph.fcs), value=TRUE))
        s1 <- s[[1]]
        s2 <- s[[2]] 
        #plot
        marginal.2D.density(lymph.fcs[[s1]], lymph.clusters.stats[[s1]], dir=sprintf('~nikolas/IL2/Plots/downsample-%d/%s/clusters-%d',desired_samples,individual.dose,k), channels=all.channels)
        marginal.2D.density(lymph.fcs[[s1]], lymph.clusters.down.stats[[s1]], dir=sprintf('~nikolas/IL2/Plots/downsample-%d/%s/down-clusters-%d',desired_samples,individual.dose,k), channels=all.channels)
        #marginal.2D.density(lymph.fcs.clust[[s]], lymph.fcs.clust[[s]], dir=sprintf('~nikolas/IL2/Plots/downsample-%d/%s/clusters-%d',down,s,k), channels=channels)
        marginal.2D.density(lymph.clusters.stats[[s1]], lymph.clusters.stats[[s2]], dir=sprintf('~nikolas/IL2/Plots/day-to-day/downsample-%d/%s/clusters-%d',desired_samples,individual.dose,k), channels=all.channels)
        marginal.2D.density(lymph.fcs[[s1]], lymph.fcs[[s2]], dir=sprintf('~nikolas/IL2/Plots/day-to-day/%s',individual.dose), channels=all.channels)
    }
}

                    

f <- lymph.fcs.clust[[1]]
pdf('~nikolas/IL2/Plots/plot3.pdf')
plot(f@exprs[,c('CD45RA','CD25')], pch='.', col=as.factor(f@exprs[,'col']))
plot(f@exprs[,c('CD4','CD25')], pch='.', col=as.factor(f@exprs[,'col']))
plot(f@exprs[,c('CD45RA','CD4')], pch='.', col=as.factor(f@exprs[,'col']))
dev.off()



X <- lymph.clusters.stats[grep('CB00165D.0U', names(lymph.clusters.stats))]
x <- X[order(X$name), c('FSC-A', 'SSC-A', 'CD25', 'CD4', 'CD45RA','pSTAT5','count', 'pct') ]
#d <- as.matrix(dist(x))

pdf('~nikolas/IL2/Plots/plot1.pdf')
s <- lymph.clusters.stats3[which(lymph.clusters.stats3$name=='2013-03-07.CB00165D.0U'),]
plot(getChannels(lymph.fcs[["2013-03-07.CB00165D.0U"]], c('CD25', 'CD45RA')), pch='', col=lymph.assign3[['2013-03-07.CB00165D.0U']])
points(s[,c('CD25','CD45RA')], pch=20, cex=s$pct*5, col='purple')
dev.off()


plot(x[,c('CD25','pSTAT5')], pch=20, col=as.factor(sort(X$name)), cex=X$pct*2)
plot(ecdf(dist(x[1:200,c('CD25','pSTAT5')])), col='black')
lines(ecdf(dist(x[200:400,c('CD25','pSTAT5')])), col='red')
plot(density(dist(x[1:200,c('CD25','pSTAT5')])), col='black')
lines(density(dist(x[200:400,c('CD25','pSTAT5')])), col='red')
dev.off()

X <- lymph.clusters.stats[grep('CB00165D', lymph.clusters.stats$name),]
x <- X[order(X$name), c('name','FSC-A', 'SSC-A', 'CD25', 'CD4', 'CD45RA','pSTAT5','count', 'pct') ]
#d <- as.matrix(dist(x))

pdf('~nikolas/IL2/Plots/plot1.pdf')
plot(x[,c('CD25','CD45RA')], pch=20, col=as.factor(sort(x$name)), cex=x$pct*2)
plot(x[,c('CD25','pSTAT5')], pch=20, col=as.factor(sort(x$name)), cex=x$pct*2)
plot(ecdf(dist(x[1:200,c('CD25','pSTAT5')])), col='black')
lines(ecdf(dist(x[200:400,c('CD25','pSTAT5')])), col='red')
plot(density(dist(x[1:200,c('CD25','pSTAT5')])), col='black')
lines(density(dist(x[200:400,c('CD25','pSTAT5')])), col='red')
dev.off()

marginal.2D.density(lymph.fcs[['2013-03-07.CB00165D.0U']], lymph.down[['2013-03-07.CB00165D.0U']], dir='~nikolas/IL2/Plots/marginal-density/day2/')
marginal.2D.density(lymph.down[['2012-11-29.CB00165D.0U']], lymph.down[['2013-03-07.CB00165D.0U']], dir='~nikolas/IL2/Plots/marginal-density/day-diff-down/')
marginal.2D.density(lymph.fcs[['2012-11-29.CB00165D.0U']], lymph.fcs[['2013-03-07.CB00165D.0U']], dir='~nikolas/IL2/Plots/marginal-density/day-diff/')
#
marginal.2D.density(lymph.down1000[['2012-11-29.CB00165D.0U']], lymph.down1000[['2013-03-07.CB00165D.0U']], dir='~nikolas/IL2/Plots/marginal-density/day-diff-down-1000/')
#
marginal.2D.density(flat.rep.fcs[[1]], flat.rep.fcs[[1]], dir='~nikolas/IL2/Plots/lymph-gate/', channels=c('FSC-A', 'CD4$'))


pdf('~nikolas/IL2/Plots/lymph-down-distances.pdf')
plot(density(lymph.dist[[1]]),ylim=c(0,.6), col='white')
#day1 solid lines
for (n in grep('2012-11-29.CB00165D\\.0U$', names(lymph.dist), value=T)) lines(density(lymph.dist[[n]]), col='black')
for (n in grep('2012-11-29.CB00165D\\.1000U$', names(lymph.dist), value=T)) lines(density(lymph.dist[[n]]), col='red')
#day2 dashed lines
for (n in grep('2013-03-07.CB00165D\\.0U$', names(lymph.dist), value=T)) lines(density(lymph.dist[[n]]), col='black', lty=2)
for (n in grep('2013-03-07.CB00165D\\.1000U$', names(lymph.dist), value=T)) lines(density(lymph.dist[[n]]), col='red', lty=2)
dev.off()

#pdf('~nikolas/IL2/Plots/graphs.pdf')
graph.adjacency(lymph.dist[['2012-11-29.CB00165D.0U']])->a1
graph.adjacency(lymph.dist[['2012-11-29.CB00165D.1000U']])->a2
graph.adjacency(lymph.dist[['2013-03-07.CB00165D.0U']])->b1
graph.adjacency(lymph.dist[['2013-03-07.CB00165D.1000U']])->b2
print(centralization.degree(a1)$centralization/centralization.degree(a2)$centralization)
print(centralization.degree(b1)$centralization/centralization.degree(b2)$centralization)

graph.adjacency(lymph.dist[['2012-10-16.CB00406Q.0U']])->a1
graph.adjacency(lymph.dist[['2012-10-16.CB00406Q.1000U']])->a2
graph.adjacency(lymph.dist[['2013-01-22.CB00406Q.0U']])->b1
graph.adjacency(lymph.dist[['2013-01-22.CB00406Q.1000U']])->b2
print(centralization.degree(a1)$centralization/centralization.degree(a2)$centralization)
print(centralization.degree(b1)$centralization/centralization.degree(b2)$centralization)



#summary stats
f <- diameter
f <- closeness
f <- function(x) determinant(graph.laplacian(x))$modulus
f <- function(x) mean(E(x)$weight)
f <- function(x) median(E(x)$weight)
f <- function(x) min(E(x)$weight)
f <- function(x) diameter(x)
f <- function(x) transitivity(x, type='undirected') #returns 0
f <- function(x) betweenness(x, directed=FALSE) #returns 0
f <- function(x) max(shortest.paths(x))

cor( sapply(graph1, f), sapply(graph2, f) ) ** 2
#cor( sapply(graph1, function(x) sum(E(x))), sapply(graph2, function(x) sum(E(x))) )
#cor( sapply(graph1, function(x) centralization.degree(x)$centralization), sapply(graph2, function(x) centralization.degree(x)$centralization) )
#cor( sapply(graph1, function(x) centralization.closeness(x)), sapply(graph2, function(x) centralization.closeness(x)) )
#across doses
qplot( x=sapply(graph[['0U']][[1]], f), y=sapply(graph[['1000U']][[1]], f) ) + geom_abline(slope=1,intercept=0)
#across days
qplot( x=sapply(graph[['0U']][[1]], f), y=sapply(graph[['0U']][[2]], f) ) + geom_abline(slope=1,intercept=0)

graph<-lymph.graph.all

#vertex pct 
graph.summary(graph, function(x) median(V(x)$pct), 'median of vertex pct', dir='~nikolas/IL2/Plots/lymph/graph/summary/pct', file.name='median')
graph.summary(graph, function(x) sd(V(x)$pct), 'sd of vertex pct', dir='~nikolas/IL2/Plots/lymph/graph/summary/pct', file.name='sd')
#edge weight properties
graph.summary(lymph.graph.all, function(x) sum(E(x)$weight), 'Sum of edge weights', dir='~nikolas/IL2/Plots/lymph/graph/summary/weight', file.name='sum')
graph.summary(lymph.graph.all, function(x) sd(E(x)$weight), 'sd of edge weights', dir='~nikolas/IL2/Plots/lymph/graph/summary/weight', file.name='sd')
graph.summary(lymph.graph.all, function(x) Mode(E(x)$weight), 'Mode of edge weights', dir='~nikolas/IL2/Plots/lymph/graph/summary/weight', file.name='mode')
graph.summary(lymph.graph.all, function(x) diameter(x), 'Diameter', dir='~nikolas/IL2/Plots/lymph/graph/summary/weight/', file.name='diameter')
#always the same value
#graph.summary(lymph.graph.all, function(x) centralization.degree(x)$centralization, 'Centralisation', dir='~nikolas/IL2/Plots/lymph/graph/summary/centralisation/', file.name='centralisation')
#same as diameter
#graph.summary(graph, function(x) max(shortest.paths(x)), 'Max shortest path')
graph.summary(graph, function(x) median(shortest.paths(x)), 'Median shortest path', dir='~nikolas/IL2/Plots/lymph/graph/summary/weight/', file.name='shortest-path')
graph.summary(graph, function(x) determinant(graph.laplacian(x, sparse=FALSE))$modulus, 'Determinant of graph Laplacian', dir='~nikolas/IL2/Plots/lymph/graph/summary/weight/', file.name='graph-laplacian')


#multiple stats
f <- function(x) closeness(x)
f <- function(x) closeness.centralisation(x)
f <- function(x) betweenness.norm(x)
f <- function(x) degree.norm(x)

f <- function(x) E(x)$weight/sum(E(x)$weight)
graph.f.density(graph=lymph.graph.all, dir='~nikolas/IL2/Plots/lymph/graph-all-distances/', f=f)
graph.f.density(graph=lymph.graph.core, dir='~nikolas/IL2/Plots/lymph/graph-core-distances/', f=f)
graph.f.density(graph=lymph.graph.func, dir='~nikolas/IL2/Plots/lymph/graph-func-distances/', f=f)

f <- function(x) V(x)$pct
graph.f.density(graph=lymph.graph.all, dir='~nikolas/IL2/Plots/lymph/graph-all-density/', f=f)
graph.f.density(graph=lymph.graph.core, dir='~nikolas/IL2/Plots/lymph/graph-core-density/', f=f)

f <- function(x) shortest.paths(x)
graph.f.density(graph=lymph.graph.all, dir='~nikolas/IL2/Plots/lymph/graph-all-shortest-paths/', f=f)
graph.f.density(graph=lymph.graph.core, dir='~nikolas/IL2/Plots/lymph/graph-core-shortest-paths/', f=f)

#down
#rep.lymph.down.fcs <- load.FCS(ext='.lymphocytes.fcs.downsample.fcs')
#flat.lymph.down.fcs <- flatten(rep.lymph.down.fcs)
#
#update.FCS('.fcs')
#rep.down.fcs <- load.FCS(ext='.fcs.downsample.fcs')


### EFFECT OF DOWNSAMPLING

# channel downsample
channel.downsample(rep.lymph.fcs, rep.lymph.down.fcs,  dir='~nikolas/IL2/Plots/downsampling-mean/lymphocytes/', fun=mean)
channel.downsample(rep.lymph.fcs, rep.lymph.down.fcs,  dir='~nikolas/IL2/Plots/downsampling-density/lymphocytes/', fun=density)





