library(flowCore)
library(tools)
library(mixtools)
library(flowStats)
library(ncdfFlow)
library(ks)
library(cluster)
library(scales)
library(reshape2)
library(flowBeads)
library(spade)
library(fastcluster)
library(ggplot2)
library(lubridate)
library(numDeriv)

setwd('~nikolas/Projects/IL2')


# my read.FCS function
source('~nikolas/bin/FCS/fcs.R') 
source('~nikolas/Projects/IL2/bin/functions.R') 
source('~nikolas/Projects/IL2/bin/downsample.R')

doses <- c('0U', '1U', '10U', '1000U')
doses <- c('0U', '1000U')
core.channels <- c('FSC-A', 'SSC-A', 'CD25', 'CD4$', 'CD45RA')
core.channels <- gsub('\\$', '', core.channels)
func.channels <- 'pSTAT5'
all.channels <- c(core.channels, func.channels)
all.channels <- gsub('\\$', '', all.channels)
dir <- '~nikolas/dunwich/Projects/IL2/FCS.repeats/'
rep.individuals <- sort(c('CB01494Y', 'CB01495Z', 'CB01498C', 'CB00406Q', 'CB00165D', 'CB01503H', 'CB00396E', 'CB01484M', 'CB00366X', 'CB01504J'))


pch <- data.frame(cbind(individual=rep.individuals, pch=letters[1:length(rep.individuals)],
                        #day1=sapply(rep.individuals, function(x) min( getDate(rep.fcs[[x]][['day1']][[1]]), getDate(rep.fcs[[x]][['day2']][[1]]) )),
                        #day2=sapply(rep.individuals, function(x) max( getDate(rep.fcs[[x]][['day1']][[1]]), getDate(rep.fcs[[x]][['day2']][[1]]) )),
                        col=sample(rainbow(length(rep.individuals)))
                        ),stringsAsFactors=F)
#pch$day.diff <- as.Date(pch[,'day2'])-as.Date(pch[,'day1'])
print(pch)

pch <- structure(list(individual = c("CB00165D", "CB00366X", "CB00396E", "CB00406Q", 
"CB01484M", "CB01494Y", "CB01495Z", "CB01498C", "CB01503H", "CB01504J"
), pch = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j"), 
    day1 = c("2012-11-29", "2012-11-07", "2012-09-25", "2012-10-16", 
    "2012-09-25", "2012-10-09", "2012-10-09", "2012-10-16", "2012-11-07", 
    "2012-11-07"), day2 = c("2013-03-07", "2013-03-27", "2013-03-11", 
    "2013-01-22", "2013-03-11", "2013-01-29", "2013-01-29", "2013-01-22", 
    "2013-03-07", "2013-03-27"), col = c("#0066FFFF", "#FF9900FF", 
    "#00FFFFFF", "#FF0099FF", "#33FF00FF", "#CCFF00FF", "#CC00FFFF", 
    "#3300FFFF", "#00FF66FF", "#FF0000FF"),
    day.diff = structure(c(98, 140, 167, 98, 167, 112, 112, 98, 120, 140),
    units = "days", class = "difftime")),
    .Names = c("individual", "pch", "day1", "day2", "col", "day.diff"),
    row.names = c("CB00165D", "CB00366X", "CB00396E", "CB00406Q", "CB01484M", "CB01494Y", "CB01495Z", "CB01498C", "CB01503H", "CB01504J"),
    class = "data.frame")

rep.fcs <- load.FCS('.fcs', doses=doses)
flat.rep.fcs <- flatten(rep.fcs)
#rep.density <- lapply(flat.rep.fcs, multi.density)
#rep.down <- mapply(downsample.FCS, flat.rep.fcs, rep.density)


rep.nonlymph.fcs <- load.FCS('.nonlymphocytes2.fcs')
nonlymph.fcs <- flatten(rep.nonlymph.fcs)


#channels <- c('FSC-A', 'SSC-A', 'CD25', 'CD4$', 'CD45RA', 'pstat5')
#channels <- c('CD25', 'CD4$', 'CD45RA')
#update.FCS('.lymphocytes.fcs')
rep.lymph.fcs <- load.FCS('.lymphocytes2.fcs', doses=doses)
lymph.fcs <- flatten(rep.lymph.fcs)
#lymph <- mapply(downsample.FCS, flat.rep.fcs, rep.density)

cd4.norm2 <- mapply(function(cd4,pstat5,pstat5.norm) cd4*cov(cd4,pstat5)/cov(cd4,pstat5.norm), cd4, pstat5, pstat5.norm )  

analyse.channel <- function(fcs, centers, channel='pstat5', X=seq(0,4,length.out=1000), normalised=FALSE) {
    d <- list()
    d$x <- lapply(fcs, function(x) lgcl(getChannels(x, channel)) )
    d$peaks <- lapply( d$x, function(x) sort(kmeans(x,centers=centers)$centers)) 
    #d$dens.deriv <- lapply(d$x, function(x) { h <- hns(x,deriv.order=1); d <- kdde(x, h=h, deriv.order=1); return(splinefun(d$eval.points, d$estimate)) }) 
    #d$peaks <- lapply( d$x, function(x) sort(pam(sample(x,1000),k=length(centers))$medoids)) 
    #if (normalised) d$x <- mapply( function(x, p) cbind(1,x) %*% coefficients(lm(centers ~ p)), d$x, d$peaks )  
    d$dens <- lapply( d$x, dens.fun )
    d$cdf <- lapply( d$x, ecdf )
    #d$dens.median <- apply(do.call('rbind',lapply(d$dens, function(dens) dens(X))),2,median)
    #d$cdf.median <- apply(do.call('rbind',lapply(d$cdf, function(cdf) cdf(X))),2,median)
    return(d)
}


plot.analysed.channel <- function(d,X,file.name, what='dens', ...) {
    png(file.name)
    plot(X, d[[what]][[1]](X), col='white', xlim=range(X), ...)
    sapply(d[[what]], function(f) lines(X, f(X), lwd=.25))
    #lines(X, d$dens.median, col='red', lwd=2)
    #mapply(function(x,f) points(x, f(x), pch=20, cex=2, col=1:length(d$peaks)),d$peaks,d$dens)
    dev.off()
}

area.between.curves <- function(d, X, doses=c('0U', '1000U')) {
    abc <- function(f1, f2, x=X) sum(abs(f1(x)-f2(x)))
    ABC <- data.frame()
    for (i in 1:nrow(pch)) {
        individual <- pch[i,]
        day1.cdf <- d$cdf[with(individual, paste(day1, individual, doses, sep='.'))]
        day2.cdf <- d$cdf[with(individual, paste(day2, individual, doses, sep='.'))]
        day1.dens <- d$dens[with(individual, paste(day1, individual, doses, sep='.'))]
        day2.dens <- d$dens[with(individual, paste(day2, individual, doses, sep='.'))]
        ABC <- rbind(ABC, data.frame(individual=individual$individual,
                                     cdf.day1=abc(day1.cdf[[1]],day1.cdf[[2]]),
                                     cdf.day2=abc(day2.cdf[[1]],day2.cdf[[2]]),
                                     dens.day1=abc(day1.dens[[1]],day1.dens[[2]]),
                                     dens.day2=abc(day2.dens[[1]],day2.dens[[2]])) )
    }
    return(ABC)
} 

plot.abc.agreement <- function(day1, day2, file.name, ...) {
    png(file.name)
    r <- range(c(day1,day2))
    plot(day1, day2, pch=pch$pch, xlim=r, ylim=r, xlab='day 1', ylab='day 2', ...)
    abline(b=1,a=0)
    r2 <- round(cor(day1, day2)**2, digits=2)
    legend('topleft', as.expression(bquote(r^2 == .(r2))))
    dev.off()
}


pstat5$spEMsymloc <- lapply(pstat5$x, spEMsymloc, mu=c(0,1.5))


X <- seq(-.5,3,length.out=1000)
centers <- c(.5,2)
pstat5 <- analyse.channel(flat.rep.fcs, centers=centers, channel='pstat5', X=X, normalised=FALSE)
pstat5.norm <- analyse.channel(flat.rep.fcs, centers=centers, channel='pstat5', X=X, normalised=TRUE)
plot.analysed.channel(pstat5, X=X, what='cdf', file.name='pstat5-cdf.png', ylab='', xlab='')
plot.analysed.channel(pstat5, X=X, what='dens', file.name='pstat5-dens.png', ylim=c(0,1.5), ylab='', xlab='')
plot.analysed.channel(pstat5, what='dens.deriv', X=X, file.name='pstat5-densderiv.png', ylim=c(-2,2), ylab='', xlab='')
plot.analysed.channel(pstat5.norm, X=X, file.name='pstat5-norm.png', ylim=c(0,1.5), ylab='', xlab='')
pstat5.abc <- area.between.curves(pstat5, X)
pstat5.norm.abc <- area.between.curves(pstat5.norm, X)
plot.abc.agreement(pstat5.abc$cdf.day1, pstat5.abc$cdf.day2, 'pstat5-abc-cdf.png')
plot.abc.agreement(pstat5.norm.abc$cdf.day1, pstat5.norm.abc$cdf.day2, 'pstat5-norm-abc-cdf.png')
plot.abc.agreement(pstat5.abc$dens.day1, pstat5.abc$dens.day2, 'pstat5-abc-dens.png')
plot.abc.agreement(pstat5.norm.abc$dens.day1, pstat5.norm.abc$dens.day2, 'pstat5-norm-abc-dens.png')

X <- seq(-.5,4.5,length.out=1000)
centers <- c(1,3,5)
cd4 <- analyse.channel(flat.rep.fcs, centers=centers, channel='cd4', X=X, normalised=FALSE)
cd4.norm <- analyse.channel(flat.rep.fcs, centers=centers, channel='cd4', X=X, normalised=TRUE)
plot.analysed.channel(cd4, X=X, file.name='cd4.png', ylim=c(0,1.5), ylab='', xlab='')
plot.analysed.channel(cd4.norm, X=X, file.name='cd4-norm.png', ylim=c(0,1.5), ylab='', xlab='')
cd4.abc <- area.between.curves(cd4, X)
cd4.norm.abc <- area.between.curves(cd4.norm, X)
plot.abc.agreement(cd4.abc$cdf.day1, cd4.abc$cdf.day2, 'cd4-abc-cdf.png')
plot.abc.agreement(cd4.norm.abc$cdf.day1, cd4.norm.abc$cdf.day2, 'cd4-norm-abc-cdf.png')
plot.abc.agreement(cd4.abc$dens.day1, cd4.abc$dens.day2, 'cd4-abc-dens.png')
plot.abc.agreement(cd4.norm.abc$dens.day1, cd4.norm.abc$dens.day2, 'cd4-norm-abc-dens.png')







png('~/cd4.cdf.png')
plot(X, cd4.cdf[[1]](X), col='white', xlim=c(0,4))
sapply(cd4.cdf, function(f) lines(X, f(X), lwd=.5))
lines(X, cd4.cdf.median, col='red', lwd=2)
dev.off()


png('~/cd4.diff.png')
plot(X, cd4.cdf[[1]](X)-cd4.cdf.median, type='l', col='white', xlim=range(X), ylim=c(-.5,.5))
sapply(cd4.cdf, function(f) lines(X, f(X)-cd4.cdf.median))
dev.off()




write.csv(cd4.dens.median, file='cd4.dens.median.csv')
write.csv(cd4.dens[[1]](X), file='cd4.dens.1.csv')

cd4.dens.median <- read.csv('~nikolas/thor/cd4.dens.median.csv')[,2]
cd4.dens.1 <- read.csv('~nikolas/thor/cd4.dens.1.csv')[,2]


png('~/cd4.norm.png')
plot(density(cd4.norm[[1]]), col='white', xlim=c(0,4))
sapply(cd4.norm, function(x) lines(density(x),lwd=.5))
dev.off()


png('~/cd4.norm2.png')
plot(density(cd4.norm2[[1]]), col='white', xlim=c(0,4))
sapply(cd4.norm2, function(x) lines(density(x),lwd=.5))
dev.off()






e <- lapply(pstat5, ecdf)
d <- lapply(pstat5, density, bw=.1)
d.f <- lapply(d, function(d) splinefun(d$x,d$y))

e.norm <- lapply(pstat5.norm, ecdf)
d.norm <- lapply(pstat5.norm, density, bw=.1)
d.f.norm <- lapply(d.norm, function(d) splinefun(d$x,d$y))

abc(d.f[[1]], d.f[[4]])
abc(d.f.norm[[1]], d.f.norm[[4]])


e <- lapply(lymph[names(lymph)[1:4]], function(x) ecdf(lgcl(getChannels(x, 'pstat5'))))
d <- lapply(lymph[names(lymph)[1:4]], function(x) density(lgcl(getChannels(x, 'pstat5')),bw=.1))
d.f <- lapply(d, function(d) splinefun(d$x,d$y))
y.max <- max(sapply(d, function(d) d$y))

lymph.pstat5 <- lapply( lymph[names(lymph)[1:4]], function(x) lgcl(getChannels(x, 'pstat5')) )

cor(ABC.norm[-2,'cdf.day1'], ABC.norm[-2,'cdf.day2'])**2
cor(ABC.unorm[-2,'cdf.day1'], ABC.unorm[-2,'cdf.day2'])**2

#plot(ABC[,'dens.day1'], ABC[,'dens.day2'],pch=pch$pch)
#abline(b=1,a=0)
#cor(ABC[,'dens.day1'], ABC[,'dens.day2'])**2



abc(aday1.dens[[1]],aday1.dens[[4]],x)
abc(aday1.cdf[[1]],aday1.cdf[[4]],x)

abc(aday2.dens[[1]],aday2.dens[[4]],x)
abc(aday2.cdf[[1]],aday2.cdf[[4]],x)

with(individual, paste(day2, individual, doses, sep='.'))

#pdf('~nikolas/lymph-dose-effect.pdf',width=10,height=5)
pdf('~nikolas/ungated-dose-effect.pdf',width=10,height=5)
par(mfrow=c(1,2))
plot(x, e[[1]](x), col='white', xlab='pSTAT5', xlim=c(-.5,3), ylab='') 
mapply(function(e,lwd,lty) lines(x,e(x),lwd=lwd,lty=lty),e,seq(1,2.5,.5),c(1,2,2,1))
legend('topleft',doses, lwd=seq(1,2.5,.5), lty=c(1,2,2,1))
plot(x, d.f[[1]](x), col='white', xlab='pSTAT5', xlim=c(-.5,3),ylim=c(0,1),ylab='') 
mapply(function(d,lwd,lty) lines(x,d(x),lwd=lwd,lty=lty),d.f,seq(1,2.5,.5),c(1,2,2,1))
legend('topleft',doses, lwd=seq(1,2.5,.5), lty=c(1,2,2,1))
dev.off()


pdf('~nikolas/ungated-lymph-dose-effect.pdf')
par(mfrow=c(1,1))
x <- seq(-.5, 3, length.out=20000)
plot(x, d.f[[1]](x), col='white', xlab='pSTAT5', xlim=c(-.5,3),ylim=c(0,1),ylab='') 
d <- lapply(flat.rep.fcs[names(flat.rep.fcs)[c(1,4)]], function(x) density(lgcl(getChannels(x, 'pstat5')),bw=.1))
d.f <- lapply(d, function(d) splinefun(d$x,d$y))
y.max <- max(sapply(d, function(d) d$y))
mapply(function(d,lwd,lty) lines(x,d(x),lwd=lwd,lty=lty),d.f,c(1,2.5),1)
d <- lapply(lymph[names(lymph)[c(1,4)]], function(x) density(lgcl(getChannels(x, 'pstat5')),bw=.1))
d.f <- lapply(d, function(d) splinefun(d$x,d$y))
y.max <- max(sapply(d, function(d) d$y))
#r <- mapply(function(a,b) length(a)/length(b), lymph[1:4], flat.rep.fcs[1:4])
mapply(function(d,lwd,lty) lines(x,d(x),lwd=lwd,lty=lty,col='red'),d.f,c(1,2.5),1) 
legend('topleft',doses[c(1,4)], lwd=c(1,2.5), lty=1)
dev.off()




lymph <- sapply(flat.rep.fcs, function(x) { cd4 <- lgcl(getChannels(x, 'cd4')); return(x[2 <  cd4 & cd4 < 2.75,]) })


d <- lapply(flat.rep.fcs[names(flat.rep.fcs)[1:4]], function(x) density(lgcl(getChannels(x, 'cd4')),bw=.1))

pdf('~nikolas/cd4.pdf')
par(mfrow=c(1,1))
plot(d[[1]], col='white')
sapply(d, function(d) lines(d))
abline(v=c(2,2.75))
dev.off()





pdf(file.path(dir, sprintf('%s.pdf',individual)))
par(mfrow=c(2,2), mar=c(2,3,2,2), mgp=c(2,1,0), las=1)
f <- function(normalised) {
abc <- data.frame(individual, normalised)
for (day in sapply(flow.data[grep(paste(individual, '0U', sep='.'), names(flow.data))] , getDate)) {
    cat(individual, day, normalised, '\n')
    x <- seq(-.5, 4, length.out=20000)
    e0 <- ecdf(getChannels(flow.data[[grep(paste(day, individual, '0U', sep='.'), names(flow.data))]], 'pstat5'))(x)
    if (normalised) 
        e1 <- ecdf(q.transforms[[paste(day,individual,sep='.')]](getChannels(flow.data[[grep(paste(day, individual, dose, sep='.'), names(flow.data))]], 'pstat5')))(x)
    else
        e1 <- ecdf(getChannels(flow.data[[grep(paste(day, individual, dose, sep='.'), names(flow.data))]], 'pstat5'))(x)
    #plot(x, e0, xlab='pStat5 signal', col='white', main=sprintf("Day %s",day),ylab="", xlim=c(-.5,3))
    plot(x, e0, col='white', xlim=c(-.5,3), xlab='', ylab='')
    #axis(side=1)
    #axis(side=2)
    lines(x, e0, col=ifelse(normalised, 'pink', 'grey'), lwd=2)
    lines(x, e1, col=ifelse(normalised, 'red', 'black'), lwd=2)
    ecdf.diff <- e0-e1
    lines(x, ecdf.diff, lty=2, col='black',lwd=2)
    addpoints(x, e0,q0,q1,0)
    addpoints(x, e1,q0,q1,1)
    #text <- c("0U", "1000U")
    #legend("right",lwd=rep(2,length(text)),col=ifelse(normalised, c('green','darkgreen'), c('pink','red')),legend=text)
    ecdf.diff.area <- abs(sum(ecdf.diff)*diff(x)[1])
    abc.pct <- round(ecdf.diff.area,digits=2)
    abc <- cbind(abc, abc.pct)
    #title(
    text(0,1,sprintf('%.2f', abc.pct)) 
    }
return(abc)
}
abc <- rbind(abc, f(FALSE))
#abc <- rbind(abc, f(TRUE))
dev.off()
}
colnames(abc) <- c('individual', 'normalised', 'abc.pct.day1', 'abc.pct.day2')
return(abc)
}


path <- '~nikolas/dunwich/Projects/IL2/FCS.repeats/'

print(files <- list.files(path=, pattern=".*.fcs", full.names=TRUE)) 
meta <- data.frame(do.call('rbind', strsplit(file_path_sans_ext(basename(files)), '-')))
colnames(meta) <- c('individual', 'day', 'dose') 
#rownames(meta) <- file_path_sans_ext(basename(files))
rownames(meta) <- basename(files)
a <- new("AnnotatedDataFrame", data=meta, varMetadata=data.frame(row.names=c('individual', 'day', 'dose'), labelDescription=c('CB number', 'date', 'IL2 dose in U')))

f  <- read.ncdfFlowSet(files=files,ncdfFile="ncfsTest.nc",isWriteSlice= TRUE,phenoData=a) 
#f  <- read.flowSet(files=files,phenoData=a, path=path)

wf <- workFlow(f)

f@colnames <- c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W", "Alexa Fluor 488-A", "PerCP-Cy5-5-A", "APC-A", "Alexa Fluor 700-A", "APC-Cy7-A", "CD45RA", "Qdot 605-A", "PE YG-A", "PE-Cy7 YG-A", "Time")

tl <- transformList('APC-A', logicleTransform(), transformationId="logicleTransform")
add(wf, tl)

pars <- c('APC-A')
print(densityplot(day~., f, channels=pars, groups=individual, scales=list(y=list(draw=F))))

norm <- normalization(normFun=function(x, parameters, ...)
                      warpSet(x, parameters, ...),
                      parameters=pars,
                      arguments=list(grouping="GroupID", monwrd=TRUE),
                      normalizationId="Warping")


add(wf, norm)

head(cb <- read.csv('individuals-date.txt',header=FALSE,col.names=c('cb','date')))
#head /ipswich/data/Immunochip/support/cbr/cbr-immunochip-2012-06-07.tab
head(cb.ichip <-read.csv('cbr-immunochip-2012-06-07.csv',header=TRUE))
cb.ichip$cb <- cb.ichip$ipid

m <- (merge(cb, cb.ichip, by='cb'))

#/ipswich/data/Immunochip/support/cbr/cbr-immunochip-2012-06-07.tab 
#/ipswich/data/Immunochip/support/cbr/ic-clean-lookup-2012-10-26.tab 
#/chiswick/data/store/immunochip/vasculitis 
#/ipswich/data/Immunochip/


