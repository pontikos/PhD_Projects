


fcs.files <- c(
      '/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KM00783A/23_July_2102_pSTAT5_D-GAP-KM00783A_0_U.fcs',
      '/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KM00782Z/23_July_2102_pSTAT5_D-GAP-KM00782Z_0_U.fcs',
      '/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KM00784B/23_July_2102_pSTAT5_D-GAP-KM00784B_0_U.fcs',
      '/dunwich/scratch/nikolas/FCS.Tony/DGAP_Treg_KA866533H/DGAP_Treg_KA866533H_0U.fcs')
K <- 3
par(mfrow=c(2,2))
title(main=K)
fun <- function(f) {
dim(d <- flowCore::read.FCS(f))
X <- d@exprs[,c(1,4)]
k <- kmeans(X,8)$cluster
smoothPlot(X,classification=k,chulls=FALSE)
}
lapply(fcs.files, fun)

ch <- function(x) (x$betweenss/(length(x$centers)-1))/(x$tot.withinss/(length(x$cluster)-length(x$centers)))
d <- flowCore::read.FCS('/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KM00782Z/23_July_2102_pSTAT5_D-GAP-KM00782Z_0_U.fcs')
X <- d@exprs[,c(1,4)]
res <- kmeans(X,K)
ch(res)

n <- 1
K.max <- 10
kmeans.res <- lapply(1:K.max, function(k) kmeans(X,k,nstart=n))
ch.values <- sapply(kmeans.res, ch)
ch.max <- which.max(ch.values)
plot(1:K.max, ch.values, col=ifelse(1:K.max==ch.max,'red','black'), xlab='Number of clusters', ylab='Variance Ratio Criterion')
lines(1:K.max, ch.values)

par(mfrow=c(3,3))
for (k in 2:K.max) smoothPlot(X,classification=kmeans.res[[k]]$cluster,chulls=FALSE,main=k)

par(mfrow=c(1,1))
k <- 6
smoothPlot(X,classification=kmeans.res[[k]]$cluster,chulls=FALSE,main=k,ellipse.lwd=5)

