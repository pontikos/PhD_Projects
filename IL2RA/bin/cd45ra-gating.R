# Here we are interested in gating the NON T REGS on CD45RA to extract MEMORY and NAIVE subsets
# The memory cell phenotypes we are interested in are:

PHENO.ratio <- '% memory'
PHENO.cd25.mfi <- PHENO.mfi <- PHENO.mef <- 'memory CD25 MEF'

library(flowBeads)
library(mixtools)
library(cluster)
library(iterators)
source('~nikolas/bin/FCS/fcs.R')

par(cex.lab=1.5, cex.main=2)

nontregs.name <- 'lymphocytes-cd4-cd127hi.fcs'
memory.name <- 'lymphocytes-cd4-cd127hi-127hi-cd45raneg.fcs'


fcsFiles <- unique(unlist(lapply(strsplit(list.files(pattern='cad.*.fcs', path='~nikolas/dunwich/Projects/IL2RA/FCS.Gated/Calli/all.gated'), '-'), function(x) x[[1]])))

#sanity check that components are always ordered
for (fcsFile in fcsFiles) { 
    print( load(file.path('~/dunwich/Projects/IL2RA/FCS.Gated/Calli/NON_T_REGS/Fitted.Mixtures/', paste('mm',fcsFile,'obj',sep='.'))) )
    if (mm$mu[1] > mm$mu[2]) print(fcsFile)
}

#sanity check that components are always ordered
#if not order them
for (fcsFile in fcsFiles) { 
    print( load(file.path('~/dunwich/Projects/IL2RA/FCS.Gated/Calli/NON_T_REGS/Fitted.Mixtures/', paste('np','mm',fcsFile,'obj',sep='.'))) )
    if (np.mm$muhat[1] > np.mm$muhat[2]) {
        print(fcsFile)
        np.mm$posteriors <- np.mm$posteriors[,2:1]
        np.mm$mu <- np.mm$mu[,2:1]
        np.mm$muhat <- np.mm$muhat[2:1]
        np.mm$lambda <- np.mm$lambda[,2:1]
        np.mm$lambdahat <- np.mm$lambdahat[2:1]
        print( save(np.mm,file=file.path('~/dunwich/Projects/IL2RA/FCS.Gated/Calli/NON_T_REGS/Fitted.Mixtures/', paste('np','mm',fcsFile,'obj',sep='.'))) )
    }
}


#beads
beads <- b <- read.csv('~nikolas/Projects/flowBeads/Beads.Stats/beads.fcs2.kmedoids.stats') 
colMedoids <- function(X) apply(X, 2, function(x) median(x))
b$date <- as.Date(b$date)
#beads.mef <- colMeans(b[,paste('mfi',2:6,sep='.')])
data(dakomef)
beads.mef <- dakomef[,'APC'][2:6]
beads <- do.call('rbind',by(b,b$date,function(x) colMeans(x[,paste('mfi',2:6,sep='.')])))


pct <- seq(.5,.99,.01)*100
d <- data.frame()
for (fcsFile in fcsFiles) { 
    nontregs.num <- nrow(nontregs<-read.FCS(nontregs.fcsFile<-paste(file.path('~nikolas/dunwich/Projects/IL2RA/FCS.Gated/Calli/all.gated',fcsFile), nontregs.name, sep='-'),channel=c('CD25','CD45RA'),TRANS=log10))
    memory.num <- nrow(memory<-read.FCS(memory.fcsFile<-paste(file.path('~nikolas/dunwich/Projects/IL2RA/FCS.Gated/Calli/all.gated',fcsFile),memory.name, sep='-'),channel=c('CD25','CD45RA'),TRANS=log10))
    cd45ra.max <- max(memory[,'CD45RA'])
    nontregs.cd45ra <- nontregs[,'CD45RA']
    memory.cd45ra <- memory[,'CD45RA']
    manual.memory.ratio <- 100*memory.num/nontregs.num
    MEF <- function(x) {
        beads.mfi <- as.numeric(beads[getDate(memory.fcsFile),paste('mfi',2:6,sep='.')])
        return(cbind(1,x)%*%coefficients(lm(log10(beads.mef) ~ log10(beads.mfi)))) 
    }
    manual.memory.cd25.mfi <- MEF(mean(memory[,'CD25']))
    #X <- nontregs.cd45ra[which(nontregs.cd45ra>0)] 
    #res.pam <- pam(X,2)
    #res <- list( mu=res.pam$medoids[,1], lambda=as.numeric(prop.table(table(res.pam$clustering))), sigsqrd=as.numeric(by(X,res.pam$clustering,var)) )
    #m <- mixtools::normalmixEM2comp( X, mu=res$mu, sigsqrd=res$sigsqrd, lambda=res$lambda ) 
    print( load(file.path('~/dunwich/Projects/IL2RA/FCS.Gated/Calli/NON_T_REGS/Fitted.Mixtures/', paste('mm',fcsFile,'obj',sep='.'))) )
    print( load(file.path('~/dunwich/Projects/IL2RA/FCS.Gated/Calli/NON_T_REGS/Fitted.Mixtures/', paste('np','mm',fcsFile,'obj',sep='.'))) )
    mm.ratio <- round(100*mm$lambda[1])
    # TODO the issue is that mm$x and nontregs are not the same length AAAARG
    #mm.cd25.mfi <- weighted.mean(nontregs[,'CD25'],mm$posterior[,1])
    spmm.ratio <- round(100*np.mm$lambdahat[1])
    m <- mm
    i <- order(m$mu)[1]
    # thresholding on posterior does poorly on individual d day 1 when we have a trimodal distribution
    cd45ra.post.gates <- sapply( pct/100, function(thresh) max(m$x[which(m$posterior[,i]>=thresh)]) )
    names(cd45ra.post.gates) <- paste('g','post',pct,sep='.')
    #
    cd45ra.pct.gates <- sapply( pct/100, function(thresh) qnorm(thresh,mean=m$mu[i],sd=m$sigma[i]) )
    names(cd45ra.pct.gates) <- paste('g','pct',pct,sep='.')
    #
    pct.memory.ratios <- sapply( cd45ra.pct.gates, function(cd45ra.gate) 100*(length(which(memory.cd45ra<cd45ra.gate)))/nontregs.num )
    names(pct.memory.ratios) <- paste('pct','ratio',pct,sep='.')
    post.memory.ratios <- sapply( cd45ra.post.gates, function(cd45ra.gate) 100*(length(which(memory.cd45ra<cd45ra.gate)))/nontregs.num )
    names(post.memory.ratios) <- paste('post','ratio',pct,sep='.')
    #
    pct.memory.cd25.mfi <- sapply( cd45ra.pct.gates, function(cd45ra.gate) MEF(mean(memory[which(memory.cd45ra<cd45ra.gate),'CD25'])) )
    names(pct.memory.cd25.mfi) <- paste('pct','cd25.mfi',pct,sep='.')
    post.memory.cd25.mfi <- sapply( cd45ra.post.gates, function(cd45ra.gate) MEF(mean(memory[which(memory.cd45ra<cd45ra.gate),'CD25'])) )
    names(post.memory.cd25.mfi) <- paste('post','cd25.mfi',pct,sep='.')
    #
    d <- rbind(d,data.frame(fcsFile=fcsFile,
                            date=getDate(memory.fcsFile),
                            manual.memory.ratio=manual.memory.ratio,
                            mm.ratio=mm.ratio,
                            spmm.ratio=spmm.ratio,
                            manual.memory.cd25.mfi=manual.memory.cd25.mfi,
                            cd45ra.max=cd45ra.max,
                            #fixed.gate=100*(length(which(memory.cd45ra<0.6964311)))/nontregs.num,
                            t(cd45ra.post.gates),
                            t(cd45ra.pct.gates),
                            t(pct.memory.ratios),
                            t(pct.memory.cd25.mfi),
                            t(post.memory.ratios),
                            t(post.memory.cd25.mfi)
                            ))
}

d.backup <- d



# best gate position agreement with manual
# select threshold for pct.thresh and post.thresh
zoom <- 5
pdf('~nikolas/GoogleDrive/PhD/Thesis/IL2RA/figures/cd45ra-gate-agreement.pdf',width=2*zoom,height=1*zoom)
par(mfrow=c(1,2))
figure.labels <- iter(paste(letters,')',sep=''))
ylab <- 'MSD with position of manual gate'
# threshold pctile
ms.diff <- sapply(paste('g','pct',pct,sep='.'), function(i) mean((d[,'cd45ra.max']-d[,i])**2))
plot(pct, ms.diff,
     ylab=ylab,
     #xlab='percentile of first component',
     xlab='pct.thresh',
     main='')
title(nextElem(figure.labels), adj=0)
i <- as.numeric(which.min(ms.diff))
points(pct[[i]],ms.diff[[i]], pch=20)
text(pct[[i]], min(ms.diff), as.character(pct[[i]]),pos=3) 
print(min.diff.pct.gate <- paste('g','pct',pct,sep='.')[[i]])
# threshold posterior
ms.diff <- sapply(paste('g','post',pct,sep='.'), function(i) mean((d[,'cd45ra.max']-d[,i])**2))
plot(pct, ms.diff,
     ylab=ylab,
     #xlab='posterior threshold of first component',
     xlab='post.thresh',
     main='')
title(nextElem(figure.labels), adj=0)
i <- as.numeric(which.min(ms.diff))
points(pct[[i]],ms.diff[[i]], pch=20)
text(pct[[i]], min(ms.diff), as.character(pct[[i]]),pos=3) 
print(min.diff.post.gate <- paste('g','post',pct,sep='.')[[i]])
dev.off()


# two methods of picking a threshold once a 2 comp mm has been fitted
# a) using the percentile of the first component
# b) use the posterior of the first component
#individual d
# noisy one: trimodal
fcsFile.noisy <- 'cad64_2008mar26_treg_i007576j_017.fcs'
fcsFile.good <- 'cad116_2008oct09_treg_i009546a_013.fcs'
pdf('~nikolas/GoogleDrive/PhD/Thesis/IL2RA/figures/cd45ra-threshold-example.pdf')
par(mfrow=c(2,2))
figure.labels <- iter(paste(letters,')',sep=''))
for (fcsFile in c(fcsFile.good, fcsFile.noisy)) {
#manual cd45ra threshold
manual.cd45ra.max <- d[grep(gsub('.fcs','',fcsFile),gsub('.fcs','',d$fcsFile)),'cd45ra.max']
print( load(file.path('~/dunwich/Projects/IL2RA/FCS.Gated/Calli/NON_T_REGS/Fitted.Mixtures/', paste('mm',fcsFile,'obj',sep='.'))) )
#threshold pctile: percentile of first component
plot(density(mm$x),xlim=c(0,2),ylim=c(0,1), main='', xlab=expression('Log '[10]*' CD45RA intensity') )
title(nextElem(figure.labels), adj=0)
f <- function(x,i,mm) mm$lambda[i]*dnorm(x,mean=mm$mu[i],sd=mm$sigma[i])
curve(f(x,1,mm), from=-1, to=2, col='red', add=TRUE, lwd=2, lty=2)
curve(f(x,2,mm), from=-1, to=2, col='green', add=TRUE, lwd=2, lty=2)
pct.thresh <- as.numeric(gsub('g.pct','',min.diff.pct.gate))
auto.cd45ra.max <- qnorm(pct.thresh, mean=mm$mu[1], sd=mm$sigma[1])
segments(x0=auto.cd45ra.max, x1=auto.cd45ra.max, y0=0, y1=mm$lambda[1]*dnorm(auto.cd45ra.max,mean=mm$mu[1],sd=mm$sigma[1]),col='red',lty=2,lwd=2)
#auto gate
segments(x0=-1, y0=0.1, x1=auto.cd45ra.max, y1=0.1, lwd=4, col='red')
#manual gate
segments(x0=-1, y0=0, x1=manual.cd45ra.max, y1=0, lwd=4, col='black')
# threshold posterior: 
plot(density(mm$x),xlim=c(0,2),ylim=c(0,1), main='', xlab=expression('Log '[10]*' CD45RA intensity') )
title(nextElem(figure.labels), adj=0)
posterior <- function(x,i,mm) mm$lambda[i]*dnorm(x,mean=mm$mu[i],sd=mm$sigma[i]) / (mm$lambda[1]*dnorm(x,mean=mm$mu[1],sd=mm$sigma[1])+mm$lambda[2]*dnorm(x,mean=mm$mu[2],sd=mm$sigma[2]))
curve(posterior(x,1,mm), from=-1, to=2, col='red', add=TRUE, lty=2, lwd=2)
curve(posterior(x,2,mm), from=-1, to=2, col='green', add=TRUE, lty=2, lwd=2)
post.thresh <- as.numeric(gsub('g.post','',min.diff.post.gate))
auto.cd45ra.max <- seq(0,2,.05)[which.min(abs(posterior(seq(0,2,.05),1,mm)-post.thresh))]
segments(x0=auto.cd45ra.max,y0=0,x1=auto.cd45ra.max,y1=post.thresh,col='red', lty=2, lwd=2)
#auto gate
segments(x0=-1, y0=0.1, x1=auto.cd45ra.max, y1=0.1, lwd=4, col='red')
#manual gate
segments(x0=-1, y0=0, x1=manual.cd45ra.max, y1=0, lwd=4, col='black')
}
dev.off()


# these are the six cases when the posterior fails to reach 99%
zoom <- 3
pdf('~nikolas/GoogleDrive/PhD/Thesis/IL2RA/figures/cd45ra-posterior-threshold-fail.pdf',width=2*zoom,height=3*zoom)
par(mfrow=c(3,2))
figure.labels <- iter(paste(letters,')',sep=''))
for (fcsFile in d[which(!is.finite(d[,'g.post.99'])),'fcsFile']) {
    manual.cd45ra.max <- d[which(d$fcsFile==fcsFile),'cd45ra.max']
    post.max <- pct[which(!is.finite(as.numeric(d[which(d$fcsFile==fcsFile),paste('g','post',pct,sep='.')])))][1]
    print( load(file.path('~/dunwich/Projects/IL2RA/FCS.Gated/Calli/NON_T_REGS/Fitted.Mixtures/', paste('mm',fcsFile,'obj',sep='.'))) )
    plot(density(mm$x),xlim=c(0,2),ylim=c(0,1), main=paste('max posterior=',post.max-1), xlab=expression('Log '[10]*' CD45RA intensity') )
    title(nextElem(figure.labels), adj=0)
    posterior <- function(x,i,mm) mm$lambda[i]*dnorm(x,mean=mm$mu[i],sd=mm$sigma[i]) / (mm$lambda[1]*dnorm(x,mean=mm$mu[1],sd=mm$sigma[1])+mm$lambda[2]*dnorm(x,mean=mm$mu[2],sd=mm$sigma[2]))
    curve(posterior(x,1,mm), from=-1, to=2, col='red', add=TRUE, lty=2, lwd=2)
    curve(posterior(x,2,mm), from=-1, to=2, col='green', add=TRUE, lty=2, lwd=2)
    abline(h=post.max/100,lty=2,lwd=.5,col='red')
    thresh <- as.numeric(gsub('g.post','',min.diff.post.gate))
    max.post <- which.min(abs(posterior(seq(0,2,.05),1,mm)-thresh))
    auto.cd45ra.max <- seq(0,2,.05)[max.post]
    segments(x0=auto.cd45ra.max,y0=0,x1=auto.cd45ra.max,y1=thresh,col='red', lty=2, lwd=2)
    #auto gate
    segments(x0=-1, y0=0.1, x1=auto.cd45ra.max, y1=0.1, lwd=4, col='red')
    #manual gate
    segments(x0=-1, y0=0, x1=manual.cd45ra.max, y1=0, lwd=4, col='black')
    print(max.post)
}
dev.off()


for (fcsFile in d[which(!is.finite(d[,'g.post.99'])),'fcsFile']) {
    manual.cd45ra.max <- d[which(d$fcsFile==fcsFile),'cd45ra.max']
    print( load(file.path('~/dunwich/Projects/IL2RA/FCS.Gated/Calli/NON_T_REGS/Fitted.Mixtures/', paste('np','mm',fcsFile,'fcs','obj',sep='.'))) )
    plot(density(np.mm$data),xlim=c(0,2),ylim=c(0,1), main='', xlab=expression('Log '[10]*' CD45RA intensity') )
    title(nextElem(figure.labels), adj=0)
    points(np.mm$data, np.mm$posterior[,1], col='red')
    points(np.mm$data, np.mm$posterior[,2], col='green')
}



#we need this merge to get the individual names
x <- read.csv('~/Projects/IL2RA/Calli_CD25bright_CBR200.csv')
x$fcsFile <- tolower(x$fcsFile)
print(dim(x <- merge(x, d, by='fcsFile')))
x$fcsFile <- gsub('.fcs','',x$fcsFile)
x$date <- x$date.y
#recalled individuals
r <- read.csv('~/Projects/IL2RA/CellPhenotypes/recalled.individuals.pch') 
#date mismatch always fucks things up on merge
print(dim(r <- merge(r[,c('individual','fcsFile','pch')], x, all.x=TRUE, all.y=FALSE)))
r <- r[order(r$individual,as.Date(r$date)),]
r$pch <- as.character(r$pch)
#get pch for individuals
print(dim(x<-merge(x, r[c(TRUE,FALSE),c('individual','pch')], all.x=TRUE)))
print(table(x$pch <- ifelse(is.na(x$pch),'x',x$pch)))
d <- x


# sensitivity of cell phenotype to threshold gate position
pdf( '~nikolas/GoogleDrive/PhD/Thesis/IL2RA/figures/cd45raneg-memory-gate-phenotype-sensitivity.pdf')
par(mfrow=c(2,2),cex.lab=1.5)
figure.labels <- iter(paste(letters,')',sep=''))
#distribution of memory CD25 MFI does not change much with gate position
boxplot(d[,grep('^pct.cd25.mfi',colnames(d))],names=pct, xlab='pct.thresh', ylab=PHENO.mef)
title(nextElem(figure.labels), adj=0)
# on the other hand as the threshold increases CD45RAneg (memory) cell pct increases
boxplot(d[,grep('^pct.ratio',colnames(d))],names=pct, xlab='pct.thresh', ylab=PHENO.ratio) 
title(nextElem(figure.labels), adj=0)
#distribution of memory CD25 MFI does not change much with gate position
boxplot(d[,grep('^post.cd25.mfi',colnames(d))],names=pct, xlab='post.thresh', ylab=PHENO.mef)
title(nextElem(figure.labels), adj=0)
# on the other hand as the threshold increases CD45RAneg (memory) cell pct increases
boxplot(d[,grep('^post.ratio',colnames(d))],names=pct, xlab='post.thresh', ylab=PHENO.ratio) 
title(nextElem(figure.labels), adj=0)
dev.off()


# memory cd25 mfi phenotype agreement
# % memory phenotype agreement
zoom <- 5
pdf( '~nikolas/GoogleDrive/PhD/Thesis/IL2RA/figures/memory-auto-manual-agreement-thresholds.pdf',height=2*zoom,width=2*zoom)
par(mfrow=c(2,2))
figure.labels <- iter(paste(letters,')',sep=''))
# memory cd25 mfi phenotype agreement
x$auto.memory.cd25.mfi <- x[,paste('pct','cd25.mfi',pct.thresh*100,sep='.')]
x$pct.cd25.mfi <- x$auto.memory.cd25.mfi
xlim <- range(x[,c('manual.memory.cd25.mfi','auto.memory.cd25.mfi')])
main <- round(cor(x$manual.memory.cd25.mfi,x$auto.memory.cd25.mfi)**2,digits=3)
plot( x$auto.memory.cd25.mfi, x$manual.memory.cd25.mfi, xlab=paste('pct.threshold', PHENO.mef), ylab=paste('manual', PHENO.mef),
     xlim=xlim, ylim=xlim, main=bquote(r^2 == .(main)), pch=x$pch )
title(nextElem(figure.labels), adj=0)
abline(b=1,a=0)
x$auto.memory.cd25.mfi <- x[,paste('post','cd25.mfi',post.thresh*100,sep='.')]
x$post.cd25.mfi <- x$auto.memory.cd25.mfi
xlim <- range(x[,c('manual.memory.cd25.mfi','auto.memory.cd25.mfi')])
main <- round(cor(x$manual.memory.cd25.mfi,x$auto.memory.cd25.mfi)**2,digits=3)
plot( x$auto.memory.cd25.mfi, x$manual.memory.cd25.mfi, xlab=paste('post.threshold', PHENO.mef), ylab=paste('manual', PHENO.mef),
     xlim=xlim, ylim=xlim, main=bquote(r^2 == .(main)), pch=x$pch )
title(nextElem(figure.labels), adj=0)
abline(b=1,a=0) 
# % memory phenotype agreement
x$auto.memory.ratio <- x[,paste('pct','ratio',pct.thresh*100,sep='.')]
x$pct.ratio <- x[,paste('pct','ratio',pct.thresh*100,sep='.')]
xlim <- range(x[,c('manual.memory.ratio','auto.memory.ratio')])
main <- round(cor(x$manual.memory.ratio,x$auto.memory.ratio)**2,digits=3)
plot( x$auto.memory.ratio, x$manual.memory.ratio, xlab=paste('pct.threshold',PHENO.ratio), ylab=paste('manual',PHENO.ratio),
     xlim=xlim, ylim=xlim, main=bquote(r^2 == .(main)), pch=x$pch )
title(nextElem(figure.labels), adj=0)
abline(b=1,a=0)
x$auto.memory.ratio <- x[,paste('post','ratio',post.thresh*100,sep='.')]
x$post.ratio <- x[,paste('post','ratio',post.thresh*100,sep='.')]
xlim <- range(x[,c('manual.memory.ratio','auto.memory.ratio')])
main <- round(cor(x$manual.memory.ratio,x$auto.memory.ratio)**2,digits=3)
plot( x$auto.memory.ratio, x$manual.memory.ratio, xlab=paste('post.threshold',PHENO.ratio), ylab=paste('manual',PHENO.ratio),
     xlim=xlim, ylim=xlim, main=bquote(r^2 == .(main)), pch=x$pch )
title(nextElem(figure.labels), adj=0)
abline(b=1,a=0)
dev.off()

# % memory phenotype agreement
zoom <- 5
pdf( '~nikolas/GoogleDrive/PhD/Thesis/IL2RA/figures/memory-auto-manual-agreement-weights.pdf',height=1*zoom,width=2*zoom)
par(mfrow=c(1,2))
figure.labels <- iter(paste(letters,')',sep=''))
# % memory phenotype agreement
xlim <- range(x[,c('manual.memory.ratio','mm.ratio')])
main <- round(cor(x$manual.memory.ratio,x$mm.ratio)**2,digits=3)
plot( x$mm.ratio, x$manual.memory.ratio, xlab=paste('mm',PHENO.ratio), ylab=paste('manual',PHENO.ratio),
     xlim=xlim, ylim=xlim, main=bquote(r^2 == .(main)), pch=x$pch )
title(nextElem(figure.labels), adj=0)
abline(b=1,a=0)
xlim <- range(x[,c('manual.memory.ratio','spmm.ratio')])
main <- round(cor(x$manual.memory.ratio,x$spmm.ratio)**2,digits=3)
plot( x$spmm.ratio, x$manual.memory.ratio, xlab=paste('spmm.',PHENO.ratio), ylab=paste('manual',PHENO.ratio),
     xlim=xlim, ylim=xlim, main=bquote(r^2 == .(main)), pch=x$pch )
title(nextElem(figure.labels), adj=0)
abline(b=1,a=0)
dev.off()


#CD45RA- gate over time
par(mfrow=c(1,1))
d <- d[order(as.Date(d$date)),]
plot(as.Date(d$date), d$cd45ra.max, pch='-', cex=2, xlab='', ylab='CD45RA- threshold')
manual.cd45ra.gate.range <- sapply(unique(as.Date(d$date)), function(x) range(d[which(as.Date(d$date)==x),'cd45ra.max']))
segments(x0=unique(as.Date(d$date)), y0=manual.cd45ra.gate.range[1,], y1=manual.cd45ra.gate.range[2,])
abline(h=median(d$cd45ra.max),lty=2, col='black')
d$pct.gate <- d[,paste('g','pct',pct.thresh*100,sep='.')]
points(as.Date(d$date), d$pct.gate, col='red')
abline(h=median(d$pct.gate),lty=2, col='red')
d$post.gate <- d[,paste('g','post',post.thresh*100,sep='.')]
points(as.Date(d$date), d$post.gate, col='blue')
abline(h=median(d$post.gate),lty=2, col='blue')

par(mfrow=c(1,1))
d <- d[order(as.Date(d$date)),]
plot(as.Date(d$date), d$cd45ra.max-d$post.gate, cex=2, xlab='', ylab='CD45RA- threshold', col='blue', pch=d$pch)
points(as.Date(d$date), d$cd45ra.max-d$pct.gate, cex=2, col='red', pch=d$pch)


#recalled individuals
r <- read.csv('~/Projects/IL2RA/CellPhenotypes/recalled.individuals.pch') 
#date mismatch always fucks things up on merge
print(dim(r <- merge(r[,c('individual','fcsFile','pch')], x, all.x=TRUE, all.y=FALSE)))
r <- r[order(r$individual,as.Date(r$date)),]
r$pch <- as.character(r$pch)
phenos <- c('manual.memory.ratio','post.ratio','pct.ratio', 'manual.memory.cd25.mfi', 'post.cd25.mfi', 'pct.cd25.mfi', 'spmm.ratio', 'mm.ratio')
R <- na.omit(cbind(r[c(TRUE,FALSE),c('individual','pch',phenos)],r[c(FALSE,TRUE),phenos]))
colnames(R) <- c("individual", "pch", paste(phenos, 'day1', sep='.'), paste(phenos, 'day2', sep='.'))

#repeatability
zoom <- 5
pdf('~/GoogleDrive/PhD/Thesis/IL2RA/figures/repeatability-memory-thresholds.pdf',height=1*zoom,width=2*zoom)
par(mfrow=c(1,2))
figure.labels <- iter(paste(letters,')',sep=''))
# cd25.mfi
xlim <- range(R[,grep('cd25.mfi',colnames(R))])
plot(R[,'manual.memory.cd25.mfi.day1'],R[,'manual.memory.cd25.mfi.day2'], pch=as.character(R$pch), xlim=xlim, ylim=xlim, xlab=paste('day 1:',PHENO.mef), ylab=paste('day 2:',PHENO.mef), cex=1.5, cex.lab=1.5)
title(nextElem(figure.labels), adj=0)
abline(b=1,a=0)
points(R[,'pct.cd25.mfi.day1'],R[,'pct.cd25.mfi.day2'], pch=as.character(R$pch), col='red', cex=1.5)
points(R[,'post.cd25.mfi.day1'],R[,'post.cd25.mfi.day2'], pch=as.character(R$pch), col='blue', cex=1.5)
rp <- vector('expression',3)
rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(R[,'manual.memory.cd25.mfi.day1'],R[,'manual.memory.cd25.mfi.day2'])**2,3))))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(R[,'pct.cd25.mfi.day1'],R[,'pct.cd25.mfi.day2'])**2,3))))[2]
rp[3] <- substitute(expression(r^2 == r2),list(r2=(round(cor(R[,'post.cd25.mfi.day1'],R[,'post.cd25.mfi.day2'])**2,3))))[2]
legend('bottomright', legend=rp, text.col=c('black','red','blue'))
# ratio
xlim <- range(R[,grep('ratio',colnames(R))])
plot(R[,'manual.memory.ratio.day1'],R[,'manual.memory.ratio.day2'], pch=as.character(R$pch), xlim=xlim, ylim=xlim, xlab=paste('day 1:',PHENO.ratio), ylab=paste('day 2:',PHENO.ratio), cex=1.5, cex.lab=1.5)
title(nextElem(figure.labels), adj=0)
abline(b=1,a=0)
points(R[,'pct.ratio.day1'],R[,'pct.ratio.day2'], pch=as.character(R$pch), col='red', cex=1.5)
points(R[,'post.ratio.day1'],R[,'post.ratio.day2'], pch=as.character(R$pch), col='blue', cex=1.5)
rp <- vector('expression',3)
rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(R[,'manual.memory.ratio.day1'],R[,'manual.memory.ratio.day2'])**2,3))))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(R[,'pct.ratio.day1'],R[,'pct.ratio.day2'])**2,3))))[2]
rp[3] <- substitute(expression(r^2 == r2),list(r2=(round(cor(R[,'post.ratio.day1'],R[,'post.ratio.day2'])**2,3))))[2]
legend('bottomright', legend=rp, text.col=c('black','red','blue'))
dev.off()
#
pdf('~/GoogleDrive/PhD/Thesis/IL2RA/figures/repeatability-memory-weights.pdf')
par(mfrow=c(1,1))
# ratio
xlim <- range(R[,grep('ratio',colnames(R))])
plot(R[,'manual.memory.ratio.day1'],R[,'manual.memory.ratio.day2'], pch=as.character(R$pch), xlim=xlim, ylim=xlim, xlab=paste('day 1:',PHENO.ratio), ylab=paste('day 2:',PHENO.ratio), cex=1.5, cex.lab=1.5)
abline(b=1,a=0)
points(R[,'mm.ratio.day1'],R[,'mm.ratio.day2'], pch=as.character(R$pch), col='red', cex=1.5)
points(R[,'spmm.ratio.day1'],R[,'spmm.ratio.day2'], pch=as.character(R$pch), col='blue', cex=1.5)
print(manual.rsquared <- round(cor(R[,'manual.memory.ratio.day1'],R[,'manual.memory.ratio.day2'])**2,3))
print(mm.rsquared <- round(cor(R[,'mm.ratio.day1'],R[,'mm.ratio.day2'])**2,3))
print(spmm.rsquared <- round(cor(R[,'spmm.ratio.day1'],R[,'spmm.ratio.day2'])**2,3))
rp <- vector('expression',3)
rp[1] <- substitute(expression(r^2 == r2),list(r2=manual.rsquared))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=mm.rsquared))[2]
rp[3] <- substitute(expression(r^2 == r2),list(r2=spmm.rsquared))[2]
legend('bottomright', legend=rp, text.col=c('black','red','blue'))
dev.off()



v <- read.csv('~nikolas/Projects/IL2RA/CellPhenotypes/vincent.cell.phenotypes')
merge(x,v,by='fcsFile')->v
plot(v$manual.memory.ratio, v$memory.freqpar, xlab='vincent', ylab='regated')
abline(b=1,a=0)

nontregs.cd45ra <- nontregs[,'CD45RA']
memory.cd45ra <- memory[,'CD45RA']
#
plot(density(nontregs.cd45ra))
lines(density(memory.cd45ra),col='red')
#
plot(normalised.density(nontregs.cd45ra[which(nontregs.cd45ra>0)]))
lines(normalised.density(memory.cd45ra[which(memory.cd45ra>0)]),col='red')

X <- nontregs.cd45ra[which(nontregs.cd45ra>0)] 
res.pam <- pam(X,2)
res <- list( mu=res.pam$medoids[,1], lambda=as.numeric(prop.table(table(res.pam$clustering))), sigsqrd=as.numeric(by(X,res.pam$clustering,var)) )
m <- mixtools::normalmixEM2comp( X, mu=res$mu, sigsqrd=res$sigsqrd, lambda=res$lambda )


nontregs.num <- nrow(nontregs<-read.FCS(nontregs.fcsFile<-paste(file.path('~nikolas/dunwich/Projects/IL2RA/FCS.Gated/Calli/all.gated',fcsFile), nontregs.name, sep='-'),channel='CD45RA',TRANS=log10))
#plot(density(nontregs[,'CD45RA'])) 
print( load(file.path('~/dunwich/Projects/IL2RA/FCS.Gated/Calli/NON_T_REGS/Fitted.Mixtures/', paste('np','mm',fcsFile,'obj',sep='.'))) )
plot(density(np.mm, component=1, block=1, scale=TRUE))
plot(np.mm, hist=FALSE)



#association tests
library(nlme)
library(lme4)
library(influence.ME)
library(beeswarm)
rm(box)
library(whisker)

x$rs12722495 <- as.numeric(x$DIL9620.CD25.rs12722495)-1
x$rs2104286 <- as.numeric(x$DIL8103.CD25.rs2104286)-1
x$rs11594656 <- as.numeric(x$DIL10847.CD25.rs11594656)-1


format.pval.col <- function(pval) {
    if (as.numeric(pval)<.05)
        sprintf('\\textcolor{red}{%s}',format.pval(pval))
    else
        format.pval(pval)
}

f <- function(phenos,covariate,x,...) {
    i <- 1
    X <- data.frame()
    fm <- formula(paste('pheno ~',covariate,sep=''))
    fm2 <- formula(paste('pheno ~',covariate,'+ (1|individual)',sep=''))
    m <- list()
    for (pheno in phenos) {
        print(dim(x1 <- na.omit(x[,c(pheno,covariate,'individual','pch')])))
        x1$col <- i
        colnames(x1) <- c('pheno',covariate,'individual','pch','col')
        m1<-lme(fm, random=~ 1|as.factor(individual), data=x1)
        m2<-lmer(fm2, data=x1)
        m[[pheno]] <- list(
                m=m1,
                m2=m2,
                cd=cooks.distance(influence(m2,obs=TRUE)),
                lower=round(as.numeric(intervals(m1)$fixed[,1][2]),digits=3),
                intercept=round(as.numeric(intervals(m1)$fixed[,2][1]),digits=3),
                fixed=round(as.numeric(intervals(m1)$fixed[,2][2]),digits=3),
                upper=round(as.numeric(intervals(m1)$fixed[,3][2]),digits=3),
                pvalue=format.pval.col(as.numeric(summary(m1)$tTable[,'p-value'][[2]])),
                col=i
            )
        X <- rbind(X, x1)
        i <- i + 1
    }
    #par(oma = c(4, .5, 1, 1), mfrow=c(1,1))
    par(mfrow=c(1,1))
    if (covariate=='age') plot(fm, data=X, col=X$col, pch=X$pch, ...)
    else beeswarm(fm, data=X, pwcol=X$col, pwpch=X$pch, cex=1.5, corral='wrap', ...)
    i <- 1
    for (m1 in m) {
        abline(b=m1$fixed, a=m1$intercept, col=unique(i), lwd=1.5)
        i <- i+1
    }
    #par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new=TRUE)
    #legend('bottom', phenos, text.col=1:length(phenos), horiz=TRUE, xpd=TRUE, inset=c(0,0), bty='n')
    legend('topleft', gsub('.ratio|.memory','',phenos), text.col=1:length(phenos))
    m
} 

covariates <- c('rs12722495', 'rs2104286', 'rs11594656', 'sex', 'age')

phenos <- c("manual.memory.cd25.mfi", "pct.cd25.mfi", "post.cd25.mfi")
m <- list()
for (covariate in covariates) {
pdf(sprintf('~nikolas/GoogleDrive/Phd/Thesis/IL2RA/figures/%s-cd25mfi.pdf',covariate))
m[[covariate]] <- f(phenos,covariate,x,ylab=PHENO.cd25.mfi,cex=.75)
dev.off()
}
params <- unlist(m)
names(params) <- gsub('\\.','_',names(params))
names(params) <- gsub('cd25_mfi_','',names(params))
names(params) <- gsub('memory_','',names(params))
names(params) <- unlist(lapply(strsplit(names(params),'_'), function(x) if (length(x)==3) paste(x[[2]],x[[1]],x[[3]],sep='_') else paste(x,collapse='_')))
cat(whisker.render(template <- readLines('~/GoogleDrive/PhD/Thesis/IL2RA/cd45ra-assoc-whisker-template.tex'),params))

phenos <- c("manual.memory.ratio", "pct.ratio", "post.ratio", "mm.ratio", "spmm.ratio")
m <- list()
for (covariate in covariates) {
pdf(sprintf('~nikolas/GoogleDrive/Phd/Thesis/IL2RA/figures/%s-ratio.pdf',covariate))
m[[covariate]] <- f(phenos,covariate,x,ylab=PHENO.ratio,cex=.75)
dev.off()
}
params <- unlist(m)
names(params) <- gsub('\\.','_',names(params))
names(params) <- gsub('ratio_','',names(params))
names(params) <- gsub('memory_','',names(params))
names(params) <- unlist(lapply(strsplit(names(params),'_'), function(x) if (length(x)==3) paste(x[[2]],x[[1]],x[[3]],sep='_') else paste(x,collapse='_')))
cat(whisker.render(template <- readLines('~/GoogleDrive/PhD/Thesis/IL2RA/cd45ra-assoc-whisker-template.tex'),params))

covariate <- 'rs2104286'
pdf(sprintf('~nikolas/GoogleDrive/Phd/Thesis/IL2RA/figures/%s-ratio-cooks-distance.pdf',covariate))
X <- na.omit(x[,c(covariate,'pch')])
xlim <- range(c(m[[covariate]][['post.ratio']]$cd, m[[covariate]][['manual.memory.ratio']]$cd))
plot( m[[covariate]][['post.ratio']]$cd, m[[covariate]][['manual.memory.ratio']]$cd , pch=X$pch, xlim=xlim, ylim=xlim, ylab="manual % memory", xlab="post.thresh % memory",
     main="Cook's distance")
abline(b=1,a=0)
dev.off()

# if d is set to the same value as with manual then association disappears
x2 <- x
x2[which(x2$pch=='d' & x2$date=='2008-03-27'),'post.ratio'] <- x2[which(x2$pch=='d' & x2$date=='2008-03-27'),'manual.memory.ratio']
phenos <- c("manual.memory.ratio", "pct.ratio", "post.ratio")
m <- list()
for (covariate in covariates) {
m[[covariate]] <- f(phenos,covariate,x2,ylab=PHENO.ratio,cex=.75)
}
params <- unlist(m)
names(params) <- gsub('\\.','_',names(params))
names(params) <- gsub('ratio_','',names(params))
names(params) <- gsub('memory_','',names(params))
names(params) <- unlist(lapply(strsplit(names(params),'_'), function(x) if (length(x)==3) paste(x[[2]],x[[1]],x[[3]],sep='_') else paste(x,collapse='_')))
cat(whisker.render(template <- readLines('~/GoogleDrive/PhD/Thesis/IL2RA/cd45ra-assoc-whisker-template.tex'),params))


M <- m[[covariate]]$mm.ratio$m 
M <- m[[covariate]]$spmm.ratio$m 
M2 <- m[[covariate]]$mm.ratio$m2 
M2 <- m[[covariate]]$spmm.ratio$m2 

covariate <- 'rs2104286'
X <- na.omit(x[,c(covariate,'pch')])
beeswarm( m[[covariate]][['post.ratio']]$cd ~ X[,covariate], pwpch=X$pch )

covariate <- 'rs2104286'
X <- na.omit(x[,c(covariate,'pch')])
plot( m[[covariate]][['post.ratio']]$cd, m[[covariate]][['manual.memory.ratio']]$cd , pch=X$pch )
abline(b=1,a=0)

covariate <- 'age'


plot(M2)
class(M)
plot( M, pch=x$pch )

plot(hist(residuals(M)))

phenos <- c("manual.memory.ratio", "pct.ratio", "post.ratio", "mm.ratio", "spmm.ratio")
covariate <- 'rs2104286'
plot(NULL, xlim=c(0,200), ylim=c(0,.1))
for (pheno in phenos) {
    infl <- influence( m[[covariate]][[pheno]]$m2 ,obs=TRUE)
    cd <- cooks.distance(infl)
    points(cd, col=m[[covariate]][[pheno]]$col, pch=x$pch)
}

infl <- influence(M2,obs=TRUE)
cd <- cooks.distance(infl)
plot(cooks.distance(infl), pch=x$pch)

plot(infl, which='cook')


