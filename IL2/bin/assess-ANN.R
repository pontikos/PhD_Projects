source('~nikolas/Projects/IL2/bin/common.R')

#individual <- 'CB00165D'
#individual <- 'CB00366X'
#individual <- 'CB01483L'
individual <- 'KM00744H'
date <- '2012-07-02'

individual <- 'CB00010K'
date <- '2012-11-13'
#date <- '2012-11-29'
#date <- '2012-11-07'
#date <- '2012-11-13'

l <- lapply(gsub('\\.','',DOSES), function(dose) read.FCS(sprintf('/chiswick/data/store/facs/Tony-FCS/PSTAT5/CD25/CD45RA/CD4/FOXP3/%s_%s_%s.fcs',individual,dose,date), channels=c('FSCW','SSCW','FSCH','SSCH',"FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3", 'PSTAT5'), TRANS=logicleTransform(w=1)))
names(l) <- DOSES
CLR <- read.csv(sprintf('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CLR/%s_0U_%s.clr',individual,date))
join <- c("FSCW", "SSCW", "FSCH", "SSCH", "FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3")
euclidean.dist <- function(x1,x2) sqrt(sum((x1-x2)**2))

DOSES <- gsub('\\.','',DOSES)
l <- list()
for (dose in DOSES) {
    load(sprintf('/chiswick/data/store/facs/Tony-RData/PSTAT5/CD25/CD45RA/CD4/FOXP3/%s_%s_%s.RData',individual,dose,date))
    #fcs.data <- applyTransforms(fcs.data, transforms)
    l[[dose]] <- applyTransforms(fcs.data, transforms)
}
print(load(sprintf('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR/%s_%s.RData',individual,date)))

load('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/KM00744H_2012-07-02.RData')
fcs.data <- applyTransforms(fcs.data, transforms)


par(mfrow=c(3,4), cex.lab=3)
figure.labels <- iter(paste(letters,')',sep=''))
for (chan in c("FSCA", "SSCA", "CD4", "CD45RA", "CD25", "FOXP3", 'PSTAT5')) {
    xquant <- quantile(transforms[[chan]](l[['0U']][,chan]),probs=seq(0,1,.01))
    print(xlim <- c(xquant[['1%']],xquant[['99%']]))
    plot(density(transforms[[chan]](l[['0U']][,chan])), xlab=chan, ylab='', main='', col='white', xlim=xlim, cex.lab=3)
    title(nextElem(figure.labels), adj=0)
    i <- 1
    for (dose in DOSES) {
        lines(density(transforms[[chan]](l[[as.character(dose)]][,chan])), lwd=i, col=blues4[i])
        i <- i+1
    }
    legend('topright',DOSES, lwd=1:4, col=blues4[1:4])
}


dim( X1 <- l[[1]][which(as.logical(CLR[,'CD4'])),join] )
dim( X1 <- l[[1]][,join] )
dim( X2 <- l[[2]][,join] )

X1 <- applyTransforms(X1, transforms)
X2 <- applyTransforms(X2, transforms)

load(sprintf('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR/%s_%s.RData',individual,date))

# nearest neigbour in X2 to X1
nn <- nn2(X2,query=X1,k=1)
nn.standard <- nn

nn <- nn2(apply(X2,2,scale),query=apply(X1,2,scale),k=1)
nn.scale <- nn

nn <- nn2(X2,query=X1,k=1, searchtype='radius', radius=.1)
nn.radius <- nn

# random test
i <- sample(1:nrow(X1),1)
nn$nn.dist[i] - euclidean.dist(X1[i,],X2[nn$nn.idx[i],])

#
f <- function(nn, n=1000) {
    # order by decreasing distance
    o <- order(nn$nn.dists,decreasing=TRUE)
    nn.dists <- nn$nn.dists[o]
    x1 <- X1[o,]
    x2 <- X2[nn$nn.idx[o],]
    par(mfrow=c(2,5))
    for (chan in join) {
        xlim <- ylim <- range(c(x1[,chan],x2[,chan]))
        #plot(x1[1:n,chan],x2[1:n,chan], pch=20, main=chan, xlab='x1', ylab='x2')
        #smoothScatter(x1[,chan],x2[,chan], main=chan, xlab='x1', ylab='x2')
        smoothPlot(cbind(x1[,chan],x2[,chan]), main=chan, xlab='', ylab='',outliers=TRUE, xlim=xlim, ylim=ylim)
        abline(b=1,a=0)
    }
}


zoom <- 5
pdf('~nikolas/GoogleDrive/PhD/Thesis/IL2/figures/ann-join-0U-01U.pdf', height=2*zoom, width=4*zoom)
f(nn.standard)
dev.off()

f(nn.scale)
f(nn.radius)



plot(density(nn.dists))


#nn$nn.dists
X2 <- as.matrix(f2[nn$nn.idx, fp])
X <- cbind(f1, X2)
colnames(X) <- c(join,paste(fp,1:2,sep='.'))

return( mean( X[,'PSTAT5.2'] - X[,'PSTAT5.1'] ) )


### join on all markers but one and see the agreement on the marker which has not been joined?
### compare the distributions


x <- rnorm(100, 1)
y <- rnorm(100, 10)

plot(density(c(x,y)))
points(cbind(y,0),col='red',pch=20,cex=.5)
points(cbind(x,0),col='blue',pch=20,cex=.5)

nn <- nn2(t(t(y)),query=t(t(x)),k=1)
points(cbind(y[nn$nn.idx],0),col='green',pch=20,cex=.5)

load('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/CB00165D_2012-11-29.RData')
f1 <- baseline.relative.pstat5(fcs.data)
load('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/CB00165D_2013-03-07.RData')
f2 <- baseline.relative.pstat5(fcs.data)
nn2(f2[,-grep('PSTAT5',colnames(f1))],query=f1[,-grep('PSTAT5',colnames(f1))],k=1)->nn

smoothPlot(cbind(f1[,'diff.PSTAT5.4'],f2[nn$nn.idx,'diff.PSTAT5.4']))
abline(b=1,a=0)


### We want to see if there is a difference in pSTAT wherever we transform/join or join/transform.
### It would appear the difference in negligeable.

join <- c("FSCW", "SSCW", "FSCH", "SSCH", "FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3")
print(load(sprintf('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR/%s_%s.RData',individual,date)))

### pSTAT5 density: transform then join
l <- lapply(gsub('\\.','',DOSES), function(dose) read.FCS(sprintf('/chiswick/data/store/facs/Tony-FCS/CD25-CD3-CD4-CD45RA-CD56-FOXP3-PSTAT5/%s_%s_%s.fcs',individual,dose,date), channels=c('FSCW','SSCW','FSCH','SSCH',"FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3", 'PSTAT5'), TRANS=logicleTransform(w=1)))
names(l) <- DOSES
#
X <- l[[1]][,join]
fp <- 'PSTAT5'
for (i in 1:length(DOSES)) {
    nn <- nn2(l[[i]][,join],query=l[[1]][,join],k=1)
    X2 <- as.matrix(l[[i]][nn$nn.idx, 'PSTAT5'])
    colnames(X2) <- paste(fp, i, sep='.')
    X <- cbind(X, X2)
}
#
X <- baseline.relative.pstat5(X)
X.trans.join <- X


### pSTAT5 density: join then transform
l <- lapply(gsub('\\.','',DOSES), function(dose) read.FCS(sprintf('/chiswick/data/store/facs/Tony-FCS/CD25-CD3-CD4-CD45RA-CD56-FOXP3-PSTAT5/%s_%s_%s.fcs',individual,dose,date), channels=c('FSCW','SSCW','FSCH','SSCH',"FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3", 'PSTAT5'), TRANS=NULL))
#
X <- l[[1]][,join]
fp <- 'PSTAT5'
for (i in 1:length(DOSES)) {
    nn <- nn2(l[[i]][,join],query=l[[1]][,join],k=1)
    X2 <- as.matrix(l[[i]][nn$nn.idx, 'PSTAT5'])
    colnames(X2) <- paste(fp, i, sep='.')
    X <- cbind(X, X2)
} 
#do transform
X <- apply(X, 2, logicleTransform(w=1))
X <- baseline.relative.pstat5(X)
X.join.trans <- X

print(load(sprintf('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/%s_%s.RData',individual,date)))
fcs.data <- apply(fcs.data, 2, logicleTransform(w=1))
fcs.data <- baseline.relative.pstat5(fcs.data)
X.join.trans2 <- fcs.data

#
gate <- 'Single cells'
plot(density(X.join.trans[which(as.logical(CLR[,gate])),'diff.PSTAT5.4']))
lines(density(X.trans.join[which(as.logical(CLR[,gate])),'diff.PSTAT5.4']), col='red')
lines(density(X.join.trans2[which(as.logical(CLR[,gate])),'diff.PSTAT5.4']), col='green')

#
plot(density(X.join.trans[,'diff.PSTAT5.4']))
lines(density(X.trans.join[,'diff.PSTAT5.4']), col='red')
lines(density(X.join.trans2[,'diff.PSTAT5.4']), col='green')


plot(density(X.join.trans[,'PSTAT5.4']))
lines(density(X.join.trans[,'PSTAT5.1']),col='red')


plot(normalised.density(logicleTransform(w=.6)(l[[4]][,'PSTAT5'])))
lines(normalised.density(logicleTransform(w=.6)(l[[1]][,'PSTAT5'])),col='red')


plot(normalised.density(logicleTransform(w=1)(l[[1]][,'PSTAT5'])), col='green')
sapply(seq(.1,2,.1), function(w) lines(normalised.density(logicleTransform(w=w)(l[[1]][,'PSTAT5']))))

r <- range((l[[1]][,'PSTAT5']))
curve(logicleTransform(w=0)(x), from=r[1], to=5000, lwd=.5)
curve(logicleTransform(w=.1)(x), from=r[1], to=5000, add=TRUE, lwd=1)
curve(logicleTransform(w=.5)(x), from=r[1], to=5000, add=TRUE, lwd=1.5)
curve(logicleTransform(w=1)(x), from=r[1], to=5000, add=TRUE, lwd=2)





