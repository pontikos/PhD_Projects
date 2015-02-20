source('~nikolas/Projects/IL2/bin/common.R')
source('~nikolas/Projects/IL2/bin/functions.R')
source('~nikolas/bin/FCS/normalise-functions.R')
library(scales)
library(reshape2)
library(flowBeads)
library(spade)
library(fastcluster)
library(ggplot2)
library(lubridate)


### used for presentation


# density day effect
channel.repeatability(rep.fcs, '~nikolas/IL2/Plots/ungated/day-effect-density/normalised/gaussNorm', fun=density, normalise=gaussNorm.normalize)
channel.repeatability(rep.fcs, '~nikolas/IL2/Plots/ungated/day-effect-density/normalised/quantiles', fun=density, normalise=quantile.normalize, quantiles=c(.25,.5,.75))
channel.repeatability(rep.fcs, '~nikolas/IL2/Plots/ungated/day-effect-density/', fun=density, normalise=quantile.normalize, quantiles=seq(0,1,.1))
channel.repeatability(rep.fcs, '~nikolas/IL2/Plots/ungated/day-effect-density/', fun=density, normalise=my.normalize)
channel.repeatability(rep.fcs, '~nikolas/IL2/Plots/ungated/day-effect-density/', fun=density, normalise=peak.normalize)
channel.repeatability(rep.fcs, '~nikolas/IL2/Plots/ungated/day-effect-density/', fun=density)
channel.repeatability(rep.fcs, '~nikolas/IL2/Plots/ungated/day-effect-scale-density/', fun=function(x) density(scale(x)))
channel.repeatability(rep.nonlymph.fcs, '~nikolas/IL2/Plots/nonlymphocytes/day-effect-density/', fun=density, normalise=TRUE)
channel.repeatability(rep.lymph.fcs, '~nikolas/IL2/Plots/lymphocytes/day-effect-density/', fun=density, normalise=quantile.normalize, quantiles=c(.25,.75))
channel.repeatability(rep.lymph.fcs, '~nikolas/IL2/Plots/lymphocytes/day-effect-density-dose-1000/', fun=density, dose=rep('1000U',2))

#does gaussNorm improve pSTAT5 repeatability?
channel.repeatability(rep.fcs, '~nikolas/IL2/Plots/ungated/day-effect-density/', fun=density, channel='pSTAT5')
channel.repeatability(rep.fcs, '~nikolas/IL2/Plots/ungated/day-effect-density/normalised/gaussNorm', fun=density, normalise=gaussNorm.normalize, channel='pSTAT5')
# yes slightly r^2 of median before norm 0.0009578821, r^2 of median after norm 0.002628568 

# mean day effect
channel.repeatability(rep.fcs, '~nikolas/IL2/Plots/day-effect-mean/', fun=mean)
channel.repeatability(rep.lymph.fcs, '~nikolas/IL2/Plots/day-effect-mean/lymphocytes', fun=mean)

# median day effect
channel.repeatability(rep.fcs, '~nikolas/IL2/Plots/day-effect-median/', fun=median)
channel.repeatability(rep.lymph.fcs, '~nikolas/IL2/Plots/day-effect-median/lymphocytes', fun=median)
channel.repeatability(rep.nonlymph.fcs, '~nikolas/IL2/Plots/day-effect-median/nonlymphocytes', fun=median)

# var day effect
channel.repeatability(rep.fcs, '~nikolas/IL2/Plots/day-effect-var/', fun=var)
channel.repeatability(rep.lymph.fcs, '~nikolas/IL2/Plots/day-effect-var/lymphocytes', fun=var)

# min day effect
channel.repeatability(rep.fcs, '~nikolas/IL2/Plots/day-effect-min/', fun=min)
channel.repeatability(rep.lymph.fcs, '~nikolas/IL2/Plots/day-effect-min/lymphocytes', fun=min)

# min day effect
channel.repeatability(rep.fcs, '~nikolas/IL2/Plots/day-effect-max/', fun=max)
channel.repeatability(rep.lymph.fcs, '~nikolas/IL2/Plots/day-effect-max/lymphocytes', fun=max)

# density dose effect
channel.repeatability(rep.fcs, '~nikolas/IL2/Plots/dose-effect-density/', day=c('day1','day1'), dose=c('0U','1000U'), fun=density)
channel.repeatability(rep.lymph.fcs, '~nikolas/IL2/Plots/dose-effect-density/lymphocytes/', day=c('day1','day1'), dose=c('0U','1000U'), fun=density)
channel.repeatability(rep.nonlymph.fcs, '~nikolas/IL2/Plots/dose-effect-density/nonlymphocytes/', day=c('day1','day1'), dose=c('0U','1000U'), fun=density)

# ecdf dose effect
channel.repeatability(rep.fcs, '~nikolas/IL2/Plots/dose-effect-ecdf/', day=rep('day1',2), dose=c('0U','1000U'), fun=ecdf)
channel.repeatability(rep.lymph.fcs, '~nikolas/IL2/Plots/dose-effect-ecdf/lymphocytes/', day=rep('day1',2), dose=c('0U','1000U'), fun=ecdf)

# ecdf dose effect
#channel.repeatability(rep.fcs, '~nikolas/IL2/Plots/dose-effect-ecdf/', day=rep('day1',2), dose=c('0U','1000U'), fun=ecdf)
channel.repeatability(rep.lymph.down.fcs, '~nikolas/IL2/Plots/dose-effect-ecdf/lymphocytes-down/', day=rep('day1',2), dose=c('0U','1000U'), fun=ecdf)

# median dose effect
channel.repeatability(rep.fcs, '~nikolas/IL2/Plots/dose-effect-median/', day=rep('day1',2), dose=c('U0','U1000'), fun=median)
channel.repeatability(rep.lymph.fcs, '~nikolas/IL2/Plots/dose-effect-median/lymphocytes/', day=rep('day1',2), dose=c('0U','1000U'), fun=median)
channel.repeatability(rep.nonlymph.fcs, '~nikolas/IL2/Plots/dose-effect-median/nonlymphocytes/', day=rep('day1',2), dose=c('0U','1000U'), fun=median)


# channel boxplot
channel.boxplot(flat.fcs, '~nikolas/IL2/Plots/boxplot/0U/', dose='^0U$')
channel.boxplot(flat.fcs, '~nikolas/IL2/Plots/boxplot/1000U/', dose='^1000U$')
channel.boxplot(flat.lymph.fcs, '~nikolas/IL2/Plots/boxplot/lymphocytes/0U/', dose='^0U$')
channel.boxplot(flat.lymph.fcs, '~nikolas/IL2/Plots/boxplot/lymphocytes/1000U/', dose='^1000U$')

# beads
load.Beads <- function(dir) {
    beads.fcs <- list()
    for (f in list.files(file.path(dir), recursive=T, pattern='beads.fcs', full.names=TRUE)) {
        print(f)
        bead.data <- gateBeads(BeadFlowFrame(f),K=6)
        beads.fcs[[f]] <- bead.data
    }
    beads.fcs
} 
beads <- load.Beads(dir='~/dunwich/FCS.Tony/repeats')
#beads cannot be used for normalisation
channel.timeseries(flat.rep.fcs, beads, dir='~nikolas/IL2/Plots/time-series/')
channel.timeseries(nonlymph.fcs, beads, dir='~nikolas/IL2/Plots/time-series-nonlymph/') 
channel.timeseries(lymph.fcs, beads, dir='~nikolas/IL2/Plots/time-series-lymph/') 


#

#\includegraphics[scale=.75]{figures/alexa488-beads.pdf}
beads.mfi <- data.frame()
for (f in list.files(file.path('~/dunwich/FCS.Tony/'), recursive=T, pattern='beads.fcs', full.names=TRUE)) {
    print(f)
    print(getDate(f))
    bead.data <- flowCore::read.FCS(f)
    print(bead.data)
    trans <- logicleTransform(w=.75)
    invtrans <- inverseLogicleTransform(trans=trans)
    x <- trans(bead.data@exprs[,'Alexa Fluor 488-A'])
    mfi <- invtrans(sort(cluster::clara(x, 6)$medoids))
    beads.mfi <- rbind(beads.mfi, data.frame(date=getDate(f), t(mfi)))
}

pstat5.mfi <- data.frame()
for (f in list.files(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/RData/'), pattern='.*_1000U_.*.RData', full.names=TRUE)) {
    print(f)
    load(f)
    f <- basename(f)
    #print(load(file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR',f)))
    print(load(file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CLR',f)))
    mfi <- sapply(CELL.TYPES, function(gate) median(fcs.data[which(as.logical(CLR[,gate])),'PSTAT5']))
    f <- gsub('.RData','',basename(f))
    print(individual <- unlist(strsplit(f,'_'))[[1]])
    print(date <- unlist(strsplit(f,'_'))[[3]])
    pstat5.mfi <- rbind(pstat5.mfi, data.frame(individual=individual,date=date,t(mfi)))
}

cd25.mfi <- data.frame()
for (f in list.files(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/RData/'), pattern='.*_1000U_.*.RData', full.names=TRUE)) {
    print(f)
    load(f)
    f <- basename(f)
    #print(load(file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR',f)))
    print(load(file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CLR',f)))
    mfi <- sapply(CELL.TYPES, function(gate) median(fcs.data[which(as.logical(CLR[,gate])),'CD25']))
    f <- gsub('.RData','',basename(f))
    print(individual <- unlist(strsplit(f,'_'))[[1]])
    print(date <- unlist(strsplit(f,'_'))[[3]])
    cd25.mfi <- rbind(cd25.mfi, data.frame(individual=individual,date=date,t(mfi)))
}

x <- melt(cd25.mfi[,3:6])
y <- melt(pstat5.mfi[,3:6])
#plot(logicleTransform(w=1)(x[,2]),logicleTransform(w=1)(y[,2]),pch=20,col=x$variable,xlab='CD25 MFI',ylab='pSTAT5 MFI')
plot((x[,2]),(y[,2]),pch=20,col=x$variable,xlab='CD25 MFI',ylab='pSTAT5 MFI')
#abline(line(x[,2], y[,2]))



#plot(as.Date(pstat5.mfi$date), pstat5.mfi[,2], col='white', ylim=c(0,max(na.omit(pstat5.mfi[,2:5]))), ylab='', xlab='')
#sapply(2:5, function(gate) points(as.Date(pstat5.mfi$date),pstat5.mfi[,gate],col=gate-1, pch=20) )
#sapply(2:5, function(gate) { lo<-loess( pstat5.mfi[,gate] ~ as.numeric(pstat5.mfi$date) ); lines(as.Date(pstat5.mfi$date),lo$fitted,col=gate-1)} )
#legend('topright', CELL.TYPES, text.col=1:4, bty='n') 
#sapply(2:4, function(i) points(as.Date(beads.mfi[,'date']),beads.mfi[,i]))

pstat5.mfi <- melt(pstat5.mfi)
pstat5.mfi$variable <- gsub('\\.',' ',pstat5.mfi$variable)

#
pdf('~nikolas/Thesis/figures/pstat5-beads.pdf',width=10,height=5)
ggplot(data=pstat5.mfi, aes(x=as.Date(date),y=value, group=variable, col=variable))+theme_bw()+geom_point()+
#geom_smooth(alpha=.1)+
stat_summary(fun.y=mean, geom='line') + ylab('pSTAT5 MFI at 1000U')+xlab('')+
geom_point(data=melt(beads.mfi[,c('date','X1','X2','X3')]), aes(x=as.Date(date),y=value), colour='orange') + 
#geom_smooth(data=melt(beads.mfi[,c('date','X1','X2','X3')]), aes(x=as.Date(date),y=value, group=variable),colour='orange',alpha=.1)+
stat_summary(fun.y=mean, geom='line', data=melt(beads.mfi[,c('date','X1','X2','X3')]), aes(x=as.Date(date),y=value, group=variable),colour='orange')+
scale_x_date(limits=range(as.Date(beads.mfi$date))) + scale_colour_manual(breaks=c(CELL.TYPES), values=c('black','red','green','blue'), guide=guide_legend(title='Cell type'))
dev.off()


# dose response normalisation
q.transforms <- quantile.transforms(flat.rep.fcs)

abc <- plot.dose.response(flat.rep.fcs, q.transforms, dir='~nikolas/IL2/Plots/ungated/dose-response')
colnames(abc) <- c('individual', 'normalised', 'abc.pct.day1', 'abc.pct.day2')
abc.nonlymph <- plot.dose.response(nonlymph.fcs, q.transforms, dir='~nikolas/IL2/Plots/nonlymph/dose-response')
colnames(abc.nonlymph) <- c('individual', 'normalised', 'abc.pct.day1', 'abc.pct.day2')
abc.lymph <- plot.dose.response(lymph.fcs, q.transforms, dir='~nikolas/IL2/Plots/lymph/dose-response')
colnames(abc.lymph) <- c('individual', 'normalised', 'abc.pct.day1', 'abc.pct.day2')

#
abc <- plot.dose.response(lymph.fcs, NULL, dose='1U', dir='~nikolas/IL2/Plots/lymphocytes/dose-response-1U')
d <- merge(abc, pch, by='individual')
plot.dose.response.agreement(d, dir='~nikolas/IL2/Plots/lymphocytes/dose-response-1U', main='dose 0.1U')
abc <- plot.dose.response(lymph.fcs, NULL, dose='10U', dir='~nikolas/IL2/Plots/lymphocytes/dose-response-10U')
d <- merge(abc, pch, by='individual')
plot.dose.response.agreement(d, dir='~nikolas/IL2/Plots/lymphocytes/dose-response-10U', main='dose 10U')
abc <- plot.dose.response(lymph.fcs, NULL, dose='1000U', dir='~nikolas/IL2/Plots/lymphocytes/dose-response-1000U')
d <- merge(abc, pch, by='individual')
plot.dose.response.agreement(d, dir='~nikolas/IL2/Plots/lymphocytes/dose-response-1000U', main='dose 1000U')



d <- merge(abc, pch, by='individual')
plot.dose.response.agreement(d, dir='~nikolas/IL2/Plots/ungated/dose-response')
d.nonlymph <- merge(abc.nonlymph, pch, by='individual')
plot.dose.response.agreement(d.nonlymph, dir='~nikolas/IL2/Plots/nonlymph/dose-response/')
d.lymph <- merge(abc.lymph, pch, by='individual')
plot.dose.response.agreement(d.lymph, dir='~nikolas/IL2/Plots/lymph/dose-response/')


x <- flowSet( lymph.fcs[grep(".CB00165D.0U", names(lymph.fcs))] )
wp <- warpSet(x, stains=c('Pacific Blue-A'))

print(densityplot(~., x, xlim=c(-1,4)))

f <- lymph.fcs[grep('.CB01494Y.0U', names(lymph.fcs))]
pdf('~nikolas/IL2/Plots/rplots.pdf')
plot(density(getChannels(f[[1]], 'cd45ra')))
lines(density(getChannels(f[[2]], 'cd45ra')),col='red')
lines(density(quantile.normalize(f[[1]], f[[2]], channel='cd45ra', quantiles=c(.25,.5,.75))),col='green')
dev.off()


#### 
setwd('~nikolas/dunwich/Projects/IL2/clean.Tony/CB00010K')

a <- read.FCS('121112_IL2_sens_T1D_TC14_pSTAT5_NKCD8-I022267C_CB00010K_0U.fcs',channels=c('pSTAT5'))
b <- read.FCS('121112_IL2_sens_T1D_TC14_pSTAT5_NKCD8-I022267C_CB00010K_1000U.fcs',channels=c('pSTAT5'))

x <- as.numeric(getChannels(a, channels='pSTAT5'))
x2 <- as.numeric(getChannels(b, channels='pSTAT5'))

plot(density(x), ylim=c(0,1), col='darkgreen')
p <- sort(pam(rep(pam(sample(x,round(length(x)/1000)),k=2)$medoids,5000),k=2)$medoids)
abline(v=p, col='darkgreen',lty=2)
lines(density(x2), ylim=c(0,1), col='red')
p2 <- sort(pam(rep(pam(sample(x2,round(length(x2)/1000)),k=2)$medoids,5000),k=2)$medoids)
abline(v=p2, col='red',lty=2)
f <- function(x) cbind(1,x)%*%coefficients(lm(p2~p))
lines(density(f(x)),col='darkgreen',lty=2)
p3 <- sort(pam(rep(pam(sample(f(x),round(length(f(x))/1000)),k=2)$medoids,5000),k=2)$medoids)
abline(v=p3, col='darkgreen',lty=2)

p<-sort(p)
p2<-sort(p2)



###
plot(NULL, xlim=c(-.75,2), ylim=c(0,0.05))
for (f in list.files(path='~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/', pattern='.*.RData', full.names=TRUE) ) {
    load(f)
    fcs.data <- baseline.relative.pstat5(fcs.data)
    lines(normalised.density(fcs.data[,'diff.PSTAT5.4']), lwd=.25)
}


sliding.window.peaks <- function(d, span=40) {
  y <- d$y
  x <- d$x
  ## sliding a window of size span and returns locations where the middle point in the window is maximum.  
  ## returns the indexes of the peaks
  ind <- c()
  for( i in 1:(length(y)-span)) {
    mid <- i+span%/%2
    if ( y[mid]==max(y[i:(i+span)]) & y[mid]!=y[i] & y[mid]!=y[i+span] ) ind <- c(ind, mid)
  }
  return(x[ind])
}


# returns the top K sliding window peaks
top.sliding.window.peaks <- function(d, K, span=40) {
  y <- d$y
  x <- d$x
  ## sliding a window of size span and returns locations where the middle point in the window is maximum.  
  ## returns the indexes of the peaks
  ind <- c()
  for( i in 1:(length(y)-span)) {
    mid <- i+span%/%2
    if ( y[mid]==max(y[i:(i+span)]) & y[mid]!=y[i] & y[mid]!=y[i+span] ) ind <- c(ind, mid)
  }
  peaks <- cbind(x=x[ind],y=y[ind])
  top.peaks <- peaks[order(peaks[,'y'],decreasing=TRUE)[1:K],]
  top.peaks <- top.peaks[order(top.peaks[,'x']),]
  return(top.peaks)
}



pstat5.peaks <- function(gate=NULL,chan='PSTAT5.4',trans=logicleTransform(w=.1)) {
    fcs.files <- list.files(path='~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/RData/pstat5-join/', pattern='.*.RData', full.names=TRUE)[1:3]
    #fcs.files <- "~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All//CB01513T_2012-11-29.RData"
    length(fcs.files)
    pstat5.dens <- pstat5.norm.dens <- matrix(0, ncol=512, nrow=length(fcs.files))
    pstat5.peaks <- data.frame()
    #plot(NULL, xlim=c(0.5,3), ylim=c(0,0.01))
    for (i in  1:length(fcs.files)) {
        print(f <- fcs.files[[i]])
        print(load(f))
        if (!is.null(gate)) {
            print(load(file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR',basename(f))))
            fcs.data <- fcs.data[which(as.logical(CLR[,gate])),]
        }
        fcs.data <- baseline.relative.pstat5(fcs.data)
        x <- fcs.data[,chan]
        x <- x[percentile.filter(x)]
        x <- sort(x)
        x <- trans(x)
        d <- normalised.density(x,from=.5,to=3,n=512)
        #unnormalised density
        pstat5.dens[i,] <- d$y
        dens <- splinefun(x=d$x,y=d$y)
        sw.peaks.x <- sort(sliding.window.peaks(d))
        #this can potentially return more than 2 peaks so always pick the one closest to one
        #and the highest of the remaining
        i1 <- which.min(abs(sw.peaks.x-1))
        sw.peaks1 <- sw.peaks.x[i1]
        sw.peaks2 <- sw.peaks.x[(i1+1):length(sw.peaks.x)]
        sw.peaks2 <- sw.peaks2[which.max(dens(sw.peaks2))]
        sw.peaks.x <- c(sw.peaks1, sw.peaks2)
        #clara is like pam but works on large datasets
        kmeans.res <- kmeans(x,centers=c(1,2))
        clara.res <- clara(x,k=2,samples=10,trace=0)
        mclust.res <- Mclust(x,G=2,modelNames='V')
        mclust.res$cluster <- mclust.res$classification
        #try all three approaches and pick the two highest peaks
        #kmeans.peaks.x <- sort( as.numeric( tapply( x, kmeans.res$cluster,function(x) sliding.window.peaks(density(x)) ) ) ) 
        kmeans.peaks.x <- tapply( x, kmeans.res$cluster,function(x) {
                                 swp <- sliding.window.peaks(density(x))
                                 return(swp[which.max(dens(swp))]) })
        #clara.peaks.x <- sort( as.numeric( tapply( x, clara.res$cluster,function(x) x[which.max(dens(x))] ) ) ) 
        clara.peaks.x <- tapply( x, clara.res$cluster,function(x) {
                                 swp <- sliding.window.peaks(density(x))
                                 return(swp[which.max(dens(swp))]) })
        #mclust.peaks.x <- sort( as.numeric( tapply( x, mclust.res$cluster,function(x) x[which.max(dens(x))] ) ) ) 
        mclust.peaks.x <- tapply( x, mclust.res$cluster,function(x) {
                                 swp <- sliding.window.peaks(density(x))
                                 return(swp[which.max(dens(swp))]) })
        #we will chose the two highest peaks (making sure they are not too close to each other)
        x.peaks <- rbind(kmeans.peaks.x, clara.peaks.x, mclust.peaks.x, sw.peaks.x)
        x.peak1 <- x.peaks[which.min(abs(x.peaks[,1]-1)),1]
        x.peak2 <- x.peaks[which.max(dens(x.peaks[,2])),2]
        x.peaks <- c(x.peak1, x.peak2)
        #x.peaks <- x.peaks[,which.min(colSums(sqrt((1:2-x.peaks)**2)))]
        #which x is closest to the local peak?
        y.peaks <- dens(x.peaks)
        pstat5.peaks <- rbind(pstat5.peaks, data.frame(peak1.x=x.peaks[[1]],peak1.y=y.peaks[[1]],peak2.x=x.peaks[[2]],peak2.y=y.peaks[[2]]))
        #map peaks to 1 and 2
        x.norm <- (x-x.peaks[[1]])/diff(x.peaks)+1
        print((x.peaks-x.peaks[[1]])/diff(x.peaks)+1)
        #normalised density
        pstat5.norm.dens[i,] <- normalised.density(x.norm,from=.5,to=3,n=512)$y
        #points(x.peaks, dens(x.peaks), pch=20, cex=2, col=1:2)
    }
    individual <- do.call('rbind', (strsplit(gsub('.RData','',basename(fcs.files)),'_')))[,1]
    date <- do.call('rbind', (strsplit(gsub('.RData','',basename(fcs.files)),'_')))[,2]
    return(list(individual=individual, date=date, x=d$x, dens=pstat5.dens, norm.dens=pstat5.norm.dens, peaks=pstat5.peaks))
} 
plot.pstat5.peaks <- function(d, main='',xlab='pSTAT5') {
    plot(NULL, xlim=c(0.5,3), ylim=range(d$dens), main=main, xlab=xlab, ylab='')
    #densities
    apply(d$dens, 1, function(y) lines(d$x, y, lwd=.5, col=alpha('black',.5)))
    #peaks
    apply(d$peaks, 1, function(d) {
          points(d['peak1.x'],d['peak1.y'],col='black',pch=20, cex=2)
          points(d['peak2.x'],d['peak2.y'],col='red',pch=20, cex=2)
        })
}
plot.norm.pstat5.peaks <- function(d, main='',xlab='pSTAT5') {
    plot(NULL, xlim=c(0.5,3), ylim=range(d$norm.dens), main=main, xlab=xlab, ylab='')
    #transformed densities
    apply(d$norm.dens, 1, function(y) {
          lines(d$x, y, lwd=.5, col=alpha('black',.5))
          f <- splinefun(d$x, y)
          points(1,f(1),col='black',pch=20, cex=2)
          points(2,f(2),col='red',pch=20, cex=2)
        })
}


ungated <- pstat5.peaks()
lymphocytes <- pstat5.peaks('Lymphocytes')
single.cells <- pstat5.peaks('Single cells')
cd4 <- pstat5.peaks('CD4',chan='PSTAT5.3')

pdf('~nikolas/GoogleDrive/PhD/Thesis/IL2/figures/pstat5-peak-normalisation.pdf',width=10,height=10)
figure.labels <- iter(paste(letters,')',sep=''))
par(mfrow=c(2,2))
#ungated unormalised
plot.pstat5.peaks(ungated,xlab='pSTAT5 at 1000U')
title(nextElem(figure.labels), adj=0)
#gated unormalised
plot.norm.pstat5.peaks(ungated,xlab='pSTAT5 at 1000U')
title(nextElem(figure.labels), adj=0) 
#cd4 lymphocytes unormalised
plot.pstat5.peaks(cd4,xlab='pSTAT5 at 10U')
title(nextElem(figure.labels), adj=0)
#cd4 lymphocytes normalised
plot.norm.pstat5.peaks(cd4,xlab='pSTAT5 at 10U')
title(nextElem(figure.labels), adj=0)
dev.off()

fcs.files <- list.files(path='~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/', pattern='.*.RData', full.names=TRUE)
norm.dens.y <- norm.dens.x <- dens.y <- dens.x <- matrix(0, ncol=512, nrow=length(fcs.files))
for (i in  1:length(fcs.files)) {
    print(f <- fcs.files[[i]])
    individual <- unlist(strsplit(gsub('.RData','',basename(f)), '_'))[[1]]
    date <- unlist(strsplit(gsub('.RData','',basename(f)), '_'))[[2]]
    print(load(file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All-pstat5-normalised/',basename(f))))
    pstat5.norm <- pstat5.normalisation[['Lymphocytes']]
    x <- pstat5.norm$x
    dens <- normalised.density(lgcl(x),n=512)
    dens.x[i,] <- dens$x
    dens.y[i,] <- dens$y
    norm.dens <- normalised.density(lgcl(pstat5.norm$pstat5.norm.transform(x)),n=512)
    norm.dens.x[i,] <- norm.dens$x
    norm.dens.y[i,] <- norm.dens$y
}

par(mfrow=c(2,1))
plot(NULL, xlim=range(dens.x), ylim=range(dens.y))
sapply(1:length(fcs.files), function(i) lines(dens.x[i,],dens.y[i,],lwd=.25,col=alpha('black',.5)))
plot(NULL, xlim=range(norm.dens.x), ylim=range(norm.dens.y))
sapply(1:length(fcs.files), function(i) lines(norm.dens.x[i,],norm.dens.y[i,],lwd=.25,col=alpha('black',.5)))


plot(normalised.density(lgcl(fcs.data[,'PSTAT5.1'])), col='white')
sapply( 1:4, function(i) lines(normalised.density(lgcl(fcs.data[,paste('PSTAT5',i,sep='.')])),col=blues4[[i]],lwd=i) ) 
plot(normalised.density(lgcl(pstat5.norm[,'PSTAT5.1'])), col='white')
sapply( 1:4, function(i) lines(normalised.density(lgcl(pstat5.norm[,paste('PSTAT5',i,sep='.')])),col=blues4[[i]],lwd=i) ) 


### example of shortcoming of pam for clustering
pdf('~nikolas/GoogleDrive/PhD/Thesis/IL2/figures/kmeans-fail.pdf')
print(load('~/kmeans-fail.RData'))
plot(normalised.density(x), xlab='pSTAT5', main='')
res.clara <- clara(x,k=2)
abline(v=max(x[res.clara$cluster==1]))
res.mclust <- Mclust(x,G=2)
abline(v=max(x[res.mclust$classification==1]),col='red')
dev.off()


### peak normalised pSTAT5 MFI per cell type
### ungated normalisation
fcs.files <- list.files(path='~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/', pattern='.*.RData', full.names=TRUE)
pstat5.mfi <- data.frame()
for (i in  1:length(fcs.files)) {
    print(f <- fcs.files[[i]])
    individual <- unlist(strsplit(gsub('.RData','',basename(f)), '_'))[[1]]
    date <- unlist(strsplit(gsub('.RData','',basename(f)), '_'))[[2]]
    print(load(file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All-pstat5-normalised/',basename(f))))
    print(load(file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR',basename(f))))
    for (gate in CELL.TYPES) {
        print(gate)
        d <- pstat5[which(as.logical(CLR[,gate])),]
        mfi <- sapply(paste('PSTAT5', 1:4, sep='.'), function(chan) {
            x <- d[,chan]
            return( median(x) )
            })
        names(mfi) <- paste('PSTAT5',1:4,sep='.')
        #d.norm <- pstat5.normalisation$Lymphocytes.10U$pstat5.norm[which(as.logical(CLR[,gate])),]
        d.norm <- pstat5.normalisation$Lymphocytes$pstat5.norm[which(as.logical(CLR[,gate])),]
        norm.mfi <- sapply(paste('PSTAT5', 1:4, sep='.'), function(chan) {
             x <- d.norm[,chan]
             return( median(x) )
            })
        names(norm.mfi) <- paste('norm','PSTAT5',1:4,sep='.')
        #d.norm <- pstat5.normalisation$Lymphocytes.1000U$pstat5.norm[which(as.logical(CLR[,gate])),]
        d.norm <- pstat5.normalisation$CD4$pstat5.norm[which(as.logical(CLR[,gate])),]
        norm2.mfi <- sapply(paste('PSTAT5', 1:4, sep='.'), function(chan) {
             x <- d.norm[,chan]
             return( median(x) )
            })
        names(norm2.mfi) <- paste('norm2','PSTAT5',1:4,sep='.')
        pstat5.mfi <- rbind( pstat5.mfi,
                            data.frame(individual=as.character(individual),
                                       date=date,
                                       cell.type=gate,
                                       t(mfi),
                                       t(norm.mfi),
                                       t(norm2.mfi),
                                       auc.pstat5=sum((mfi)),
                                       auc.norm.pstat5=sum((norm.mfi)),
                                       auc.norm2.pstat5=sum((norm2.mfi)))
                            )
    }
}

nn.peak.pstat5mfi <- pstat5.mfi


### definition of pSTAT5 AUC
zoom <- 5
pdf('~nikolas/GoogleDrive/PhD/Thesis/IL2/figures/pstat5-auc-celltypes.pdf', height=1*zoom, width=2*zoom)
#
par(mfrow=c(1,2))
figure.labels <- iter(paste(letters,')',sep=''))
plot(NULL, xlim=c(0,3), ylim=range(pstat5.mfi[,paste('PSTAT5',1:4,sep='.')]), xaxt='n', xlab='dose', ylab='pSTAT5 MFI')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
for (cell.type in CELL.TYPES) {
    apply(pstat5.mfi[grep(cell.type,pstat5.mfi$cell.type),paste('PSTAT5',1:4,sep='.')], 1, function(x) lines(0:3, x, col=i, lwd=.25))
    i <- i+1
}
legend('topleft', CELL.TYPES, text.col=1:4, bty='n')
plot(NULL, xlim=c(0,3), ylim=range(pstat5.mfi[,paste('PSTAT5',1:4,sep='.')]), xaxt='n', xlab='dose', ylab='pSTAT5 MFI')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
colMedians <- function(x) apply(x, 2, median)
for (cell.type in CELL.TYPES) {
    lines(0:3, colMedians( pstat5.mfi[grep(cell.type,pstat5.mfi$cell.type),paste('PSTAT5',1:4,sep='.')] ), col=i, lwd=2)
    i <- i+1
}
legend('topleft', CELL.TYPES, text.col=1:4, bty='n')
dev.off()


### boxplot of pSTAT5 AUC per cell type: unormalised vs normalised
pdf('~nikolas/GoogleDrive/PhD/Thesis/IL2/figures/pstat5-auc-boxplots-celltypes.pdf',width=10,height=10)
par(mfrow=c(2,2))
for (cell.type in CELL.TYPES) {
    d <- pstat5.mfi[which(pstat5.mfi$cell.type==cell.type),]
    # across all doses
    boxplot( pSTAT5.MFI ~ type, data=rbind( data.frame('pSTAT5 MFI'=d[,'auc.pstat5'], type='unnormalised')
                                            ,data.frame('pSTAT5 MFI'=d[,'auc.norm.pstat5'], type='normalised\nin lymphocytes')
                                            ,data.frame('pSTAT5 MFI'=d[,'auc.norm2.pstat5'], type='normalised\nin CD4+')
                                            ),
            col=c('white','red','blue')
            ,outline=FALSE
            ,ylab='pSTAT5 AUC'
            ,main=cell.type )
}
dev.off()

REPEATS <- rbind(cbind(REPEATS,date=REPEATS$day1),cbind(REPEATS,date=REPEATS$day2))

### repeatability of pSTAT5 AUC per cell type with normalisation in ungated (red), normalisation in CD4+ (blue) and without (black) normalisation
pdf('~nikolas/GoogleDrive/PhD/Thesis/IL2/figures/pstat5-auc-repeatability-celltypes.pdf',width=10,height=10)
par(mfrow=c(2,2))
for (cell.type in CELL.TYPES) {
    dim(repeats <- merge(pstat5.mfi[which(pstat5.mfi$cell.type==cell.type),], REPEATS))
    #
    pstat5 <- paste('PSTAT5',1:4,sep='.')
    norm.pstat5 <- paste('norm','PSTAT5',1:4,sep='.')
    norm2.pstat5 <- paste('norm2','PSTAT5',1:4,sep='.')
    #
    pstat5 <- cbind(rowSums(abs(repeats[c(TRUE,FALSE),pstat5])),rowSums(abs(repeats[c(FALSE,TRUE),pstat5])))
    norm.pstat5 <- cbind(rowSums(abs(repeats[c(TRUE,FALSE),norm.pstat5])),rowSums(abs(repeats[c(FALSE,TRUE),norm.pstat5])))
    norm2.pstat5 <- cbind(rowSums(abs(repeats[c(TRUE,FALSE),norm2.pstat5])),rowSums(abs(repeats[c(FALSE,TRUE),norm2.pstat5])))
    #
    xlim <- range(c(pstat5, norm2.pstat5, norm.pstat5))
    xlim <- range(c(pstat5, norm.pstat5))
    plot(pstat5, xlim=xlim, ylim=xlim, pch=repeats[c(TRUE,FALSE),'pch'], xlab='day 1: pSTAT5 AUC', ylab='day 2: pSTAT5 AUC', col='black', main=cell.type)
    points( norm.pstat5, pch=repeats[c(TRUE,FALSE),'pch'], col='red')
    points( norm2.pstat5, pch=repeats[c(TRUE,FALSE),'pch'], col='blue')
    rp <- vector('expression',2)
    rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(pstat5)[1,2]**2,3))**2,3))[2]
    rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(norm.pstat5)[1,2]**2,3))**2,3))[2]
    rp[3] <- substitute(expression(r^2 == r2),list(r2=(round(cor(norm2.pstat5)[1,2]**2,3))**2,3))[2]
    legend('topleft', legend=rp, text.col=c('black','red','blue'), bty='n')
    abline(b=1, a=0)
}
dev.off()


### boxplot of pSTAT5 MFI per cell type: unormalised vs normalised
for (cell.type in c('Memory Eff', 'Memory Treg', 'Naive Eff', 'Naive Treg')) {
    d <- pstat5.mfi[which(pstat5.mfi$cell.type==cell.type),]
    pdf(sprintf('~nikolas/GoogleDrive/PhD/Thesis/IL2/figures/pstat5-mfi-%s.pdf',gsub(' ','-',cell.type)),width=10,height=10)
    par(mfrow=c(2,2))
    for (i in 2:4) {
    #
    pstat5 <- paste('PSTAT5',i,sep='.')
    norm.pstat5 <- paste('norm','PSTAT5',i,sep='.')
    boxplot( pSTAT5.MFI ~ type, data=rbind(data.frame('pSTAT5 MFI'=d[,pstat5],type='unnormalised'), data.frame('pSTAT5 MFI'=d[,norm.pstat5],type='normalised')), main=cell.type )
    }
    # across all doses
    boxplot( pSTAT5.MFI ~ type, data=rbind( data.frame('pSTAT5 MFI'=rowSums(abs(d[,paste('PSTAT5',1:4,sep='.')])),type='unnormalised'),
                                            data.frame('pSTAT5 MFI'=rowSums(abs(d[,paste('norm2','PSTAT5',1:4,sep='.')])),type='normalised in ungated'),
                                            data.frame('pSTAT5 MFI'=rowSums(abs(d[,paste('norm','PSTAT5',1:4,sep='.')])),type='normalised in CD4+')),
            col=c('black','red','blue'),
            main=cell.type )
    dev.off()
}



### repeatability per cell type with (red) and without (black) normalisation
for (cell.type in c('Memory Eff', 'Memory Treg', 'Naive Eff', 'Naive Treg')) {
    dim(repeats <- merge(pstat5.mfi[which(pstat5.mfi$cell.type==cell.type),], REPEATS))
    #pdf(sprintf('~nikolas/GoogleDrive/PhD/Thesis/IL2/figures/repeatability-norm-pstat5-mfi-%s.pdf',gsub(' ','-',cell.type)),width=10,height=10)
    par(mfrow=c(2,2))
    for (i in 1:4) {
    #
    pstat5 <- paste('PSTAT5',i,sep='.')
    norm.pstat5 <- paste('norm','PSTAT5',i,sep='.')
    xlim <- range(c(repeats[,pstat5],repeats[,norm.pstat5]))
    plot(repeats[c(TRUE,FALSE),pstat5],repeats[c(FALSE,TRUE),pstat5], xlim=xlim, ylim=xlim, pch=repeats[c(TRUE,FALSE),'pch'],
         xlab=paste('day 1:', 'pSTAT5 MFI at dose',DOSES[[i]]), ylab=paste('day 2:','pSTAT5 MFI at dose',DOSES[[i]]), col='black')
    points( repeats[c(TRUE,FALSE),norm.pstat5], repeats[c(FALSE,TRUE),norm.pstat5], pch=repeats[c(TRUE,FALSE),'pch'], col='red')
    rp <- vector('expression',2)
    rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(repeats[c(TRUE,FALSE),pstat5],repeats[c(FALSE,TRUE),pstat5])**2,3))**2,3))[2]
    rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(repeats[c(TRUE,FALSE),norm.pstat5],repeats[c(FALSE,TRUE),norm.pstat5])**2,3))**2,3))[2]
    legend('topleft', legend=rp, text.col=c('black','red'), bty='n')
    abline(b=1, a=0)
    }
    #dev.off()
}



### repeatability per cell type with (red) and without (black) normalisation
for (cell.type in c('Memory Eff', 'Memory Treg', 'Naive Eff', 'Naive Treg')) {
    dim(repeats <- merge(pstat5.mfi[which(pstat5.mfi$cell.type==cell.type),], REPEATS))
    pdf(sprintf('~nikolas/GoogleDrive/PhD/Thesis/IL2/figures/repeatability-norm-pstat5-mfi-%s.pdf',gsub(' ','-',cell.type)),width=10,height=10)
    par(mfrow=c(2,2))
    for (i in 2:4) {
    #
    pstat5 <- paste('PSTAT5',i,sep='.')
    norm.pstat5 <- paste('norm','PSTAT5',i,sep='.')
    xlim <- range(c(repeats[,pstat5],repeats[,norm.pstat5]))
    plot(repeats[c(TRUE,FALSE),pstat5],repeats[c(FALSE,TRUE),pstat5], xlim=xlim, ylim=xlim, pch=repeats[c(TRUE,FALSE),'pch'],
         xlab=paste('day 1:', 'pSTAT5 MFI at dose',DOSES[[i]]), ylab=paste('day 2:','pSTAT5 MFI at dose',DOSES[[i]]), col='black')
    points( repeats[c(TRUE,FALSE),norm.pstat5], repeats[c(FALSE,TRUE),norm.pstat5], pch=repeats[c(TRUE,FALSE),'pch'], col='red')
    rp <- vector('expression',2)
    rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(repeats[c(TRUE,FALSE),pstat5],repeats[c(FALSE,TRUE),pstat5])**2,3))**2,3))[2]
    rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(repeats[c(TRUE,FALSE),norm.pstat5],repeats[c(FALSE,TRUE),norm.pstat5])**2,3))**2,3))[2]
    legend('topleft', legend=rp, text.col=c('black','red'), bty='n')
    abline(b=1, a=0)
    }
    #
    pstat5 <- paste('PSTAT5',1:4,sep='.')
    norm.pstat5 <- paste('norm','PSTAT5',1:4,sep='.')
    pstat5.day1 <- rowSums(repeats[c(TRUE,FALSE),pstat5])
    pstat5.day2 <- rowSums(repeats[c(FALSE,TRUE),pstat5])
    norm.pstat5.day1 <- rowSums(repeats[c(TRUE,FALSE),norm.pstat5])
    norm.pstat5.day2 <- rowSums(repeats[c(FALSE,TRUE),norm.pstat5])
    xlim <- range(c(pstat5.day1, pstat5.day2, norm.pstat5.day1, norm.pstat5.day2))
    plot(pstat5.day1, pstat5.day2, xlim=xlim, ylim=xlim, pch=repeats[c(TRUE,FALSE),'pch'], xlab='day 1: pSTAT5 AUC', ylab='day 2: pSTAT5 AUC', col='black')
    points( norm.pstat5.day1, norm.pstat5.day2, pch=repeats[c(TRUE,FALSE),'pch'], col='red')
    rp <- vector('expression',2)
    rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(pstat5.day1,pstat5.day2)**2,3))**2,3))[2]
    rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(norm.pstat5.day1,norm.pstat5.day2)**2,3))**2,3))[2]
    legend('topleft', legend=rp, text.col=c('black','red'), bty='n')
    abline(b=1, a=0)
    dev.off()
}


### pSTAT5 MFI response per cell-subset unormalised (top row) vs normalised (bottom)
zoom <- 5
pdf('~nikolas/GoogleDrive/PhD/Thesis/IL2/figures/norm-pstat5-response-cellsubsets.pdf', height=2*zoom, width=2*zoom)
#
par(mfrow=c(2,2))
figure.labels <- iter(paste(letters,')',sep=''))
plot(NULL, xlim=c(0,3), ylim=range(pstat5.mfi[,paste('PSTAT5',1:4,sep='.')]), xaxt='n', xlab='dose', ylab='pSTAT5 MFI')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
for (cell.type in CELL.TYPES) {
    apply(pstat5.mfi[grep(cell.type,pstat5.mfi$cell.type),paste('PSTAT5',1:4,sep='.')], 1, function(x) lines(0:3, x, col=i, lwd=.25))
    i <- i+1
}
legend('topleft', CELL.TYPES, text.col=1:4, bty='n')
plot(NULL, xlim=c(0,3), ylim=range(pstat5.mfi[,paste('PSTAT5',1:4,sep='.')]), xaxt='n', xlab='dose', ylab='pSTAT5 MFI')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
colMedians <- function(x) apply(x, 2, median)
for (cell.type in CELL.TYPES) {
    lines(0:3, colMedians( pstat5.mfi[grep(cell.type,pstat5.mfi$cell.type),paste('PSTAT5',1:4,sep='.')] ), col=i, lwd=2)
    i <- i+1
}
legend('topleft', CELL.TYPES, text.col=1:4, bty='n')
#
ylim <- range(pstat5.mfi[grep(cell.type,pstat5.mfi$cell.type),paste('norm','PSTAT5',1:4,sep='.')])
plot(NULL, xlim=c(0,3), ylim=ylim, xaxt='n', xlab='dose', ylab='normalised pSTAT5 MFI')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
for (cell.type in CELL.TYPES) {
    apply(pstat5.mfi[grep(cell.type,pstat5.mfi$cell.type),paste('norm','PSTAT5',1:4,sep='.')], 1, function(x) lines(0:3, x, col=i, lwd=.25))
    i <- i+1
}
legend('topleft', CELL.TYPES, text.col=1:4, bty='n')
plot(NULL, xlim=c(0,3), ylim=ylim, xaxt='n', xlab='dose', ylab='normalised pSTAT5 MFI')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
colMedians <- function(x) apply(x, 2, median)
for (cell.type in CELL.TYPES) {
    lines(0:3, colMedians( pstat5.mfi[grep(cell.type,pstat5.mfi$cell.type),paste('norm','PSTAT5',1:4,sep='.')] ), col=i, lwd=2)
    i <- i+1
}
legend('topleft', CELL.TYPES, text.col=1:4, bty='n')
dev.off()



#
#plot(pstat5.mfi.repeats[c(TRUE,FALSE),'PSTAT5.4'],pstat5.mfi.repeats[c(FALSE,TRUE),'PSTAT5.4'],pch=pstat5.mfi.repeats[c(TRUE,FALSE),'pch'])
cor(pstat5.mfi.repeats[c(TRUE,FALSE),'PSTAT5.4'],pstat5.mfi.repeats[c(FALSE,TRUE),'PSTAT5.4'])**2
#
#plot(pstat5.mfi.repeats[c(TRUE,FALSE),'norm.PSTAT5.4'],pstat5.mfi.repeats[c(FALSE,TRUE),'norm.PSTAT5.4'],pch=pstat5.mfi.repeats[c(TRUE,FALSE),'pch'])
cor(pstat5.mfi.repeats[c(TRUE,FALSE),'norm.PSTAT5.4'],pstat5.mfi.repeats[c(FALSE,TRUE),'norm.PSTAT5.4'])**2
#
#plot(pstat5.mfi.repeats[c(TRUE,FALSE),'PSTAT5.3'],pstat5.mfi.repeats[c(FALSE,TRUE),'PSTAT5.3'],pch=20)
cor(pstat5.mfi.repeats[c(TRUE,FALSE),'PSTAT5.3'],pstat5.mfi.repeats[c(FALSE,TRUE),'PSTAT5.3'])**2
#
#plot(pstat5.mfi.repeats[c(TRUE,FALSE),'norm.PSTAT5.3'],pstat5.mfi.repeats[c(FALSE,TRUE),'norm.PSTAT5.3'],pch=20)
cor(pstat5.mfi.repeats[c(TRUE,FALSE),'norm.PSTAT5.3'],pstat5.mfi.repeats[c(FALSE,TRUE),'norm.PSTAT5.3'])**2
#
#plot(pstat5.mfi.repeats[c(TRUE,FALSE),'PSTAT5.2'],pstat5.mfi.repeats[c(FALSE,TRUE),'PSTAT5.2'],pch=pstat5.mfi.repeats[c(TRUE,FALSE),'pch'])
cor(pstat5.mfi.repeats[c(TRUE,FALSE),'PSTAT5.2'],pstat5.mfi.repeats[c(FALSE,TRUE),'PSTAT5.2'])**2
#
#plot(pstat5.mfi.repeats[c(TRUE,FALSE),'norm.PSTAT5.2'],pstat5.mfi.repeats[c(FALSE,TRUE),'norm.PSTAT5.2'],pch=pstat5.mfi.repeats[c(TRUE,FALSE),'pch'])
cor(pstat5.mfi.repeats[c(TRUE,FALSE),'norm.PSTAT5.2'],pstat5.mfi.repeats[c(FALSE,TRUE),'norm.PSTAT5.2'])**2
}





w=(m-log10(t/abs(r))) / 2

dy <- diff(d$y)
dy2 <- diff(dy)
dy3 <- diff(dy2)
dy4 <- diff(dy3)
n <- length(dy4)

X <- cbind(dy[1:n],dy2[1:n],dy3[1:n],dy4[1:n])

order(abs(X[,1]))
(order(X[,2]))
order(abs(X[,3]))
(order(X[,4],decreasing=T))

#kmeans on whole data or pam on subset

#sliding window approach on density function

#





