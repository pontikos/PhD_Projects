source('~nikolas/bin/FCS/fcs.R')
source('~nikolas/Projects/IL2/bin/common.R')

dens.fun <- function(x,...) {
    d <- density(x,...)
    return(splinefun(d$x,d$y))
}

analyse.channel <- function(fcs, centers, channel='PSTAT5', X=seq(0,4,length.out=1000), normalised=FALSE) {
    d <- list()
    d$x <- lapply(fcs, function(x) x[,'PSTAT5'])
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

area.between.curves <- function(d, X, doses=c('0U', '1000U')) {
    abc <- function(f1, f2, x=X) mean((f1(x)-f2(x)))
    ABC <- data.frame()
    for (i in 1:nrow(REPEATS)) {
        individual <- REPEATS[i,]
        day1.cdf <- d$cdf[with(individual, paste('day1', individual, doses, sep='.'))]
        day2.cdf <- d$cdf[with(individual, paste('day2', individual, doses, sep='.'))]
        day1.dens <- d$dens[with(individual, paste('day1', individual, doses, sep='.'))]
        day2.dens <- d$dens[with(individual, paste('day2', individual, doses, sep='.'))]
        ABC <- rbind(ABC, data.frame(individual=individual$individual,
                                     cdf.day1=abc(day1.cdf[[1]],day1.cdf[[2]]),
                                     cdf.day2=abc(day2.cdf[[1]],day2.cdf[[2]]),
                                     dens.day1=abc(day1.dens[[1]],day1.dens[[2]]),
                                     dens.day2=abc(day2.dens[[1]],day2.dens[[2]])) )
    }
    return(ABC)
} 


mean.nn <- function(f1, f2, join, fp) {
    X1 <- f1[,join]
    X2 <- f2[,join]
    # nearest neigbour in X2 to X1
    nn <- nn2(X2,query=X1,k=1)
    #nn$nn.dists
    X2 <- as.matrix(f2[nn$nn.idx, fp])
    X <- cbind(f1, X2)
    colnames(X) <- c(join,paste(fp,1:2,sep='.'))
    return( mean( X[,'PSTAT5.2'] - X[,'PSTAT5.1'] ) )
}


nn.mean.diff <- function(l, join=c("FSCW", "SSCW", "FSCH", "SSCH", "FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3"), fp='PSTAT5', doses=c('0U', '1000U')) {
    ABC <- data.frame()
    for (i in 1:nrow(REPEATS)) {
        individual <- REPEATS[i,]
        f1 <- l[with(individual, paste('day1', individual, doses, sep='.'))]
        f2 <- l[with(individual, paste('day2', individual, doses, sep='.'))]
        ABC <- rbind(ABC, data.frame(individual=individual$individual,
                                     mean.nn.day1=mean.nn(f1[[1]], f1[[2]], join, fp),
                                     mean.nn.day2=mean.nn(f2[[1]], f2[[2]], join, fp)))
    }
    return(ABC)
}


plot.abc.agreement <- function(day1, day2, file.name, ...) {
    #png(file.name)
    r <- range(c(day1,day2))
    plot(day1, day2, pch=REPEATS$pch, xlim=r, ylim=r, xlab='day 1', ylab='day 2', ...)
    abline(b=1,a=0)
    r2 <- round(cor(day1, day2)**2, digits=2)
    legend('topleft', as.expression(bquote(r^2 == .(r2))))
    #dev.off()
}



l <- list()
for (i in 1:nrow(REPEATS)) {
    individual <- REPEATS[i,'individual']
    day1 <- REPEATS[i,'day1']
    day2 <- REPEATS[i,'day2']
    for (dose in DOSES) {
        #
        individual_dose_day1 <- paste(individual, dose, day1, sep='_')
        load( f1 <- file.path(BASE.DIR, sprintf('%s.RData',individual_dose_day1)) )
        l[[paste('day1',individual,dose,sep='.')]] <- fcs.data
        #
        individual_dose_day2 <- paste(individual, dose, day2, sep='_')
        load( f2 <- file.path(BASE.DIR, sprintf('%s.RData',individual_dose_day2)) )
        l[[paste('day2',individual,dose,sep='.')]] <- fcs.data
    }
}


X <- seq(-.5,3,length.out=1000)
pstat5 <- analyse.channel(l, centers=2, X=X)

par(mfrow=c(2,2))
for (dose in DOSES) {
pstat5.abc <- area.between.curves(pstat5, X=X, doses=c('0U',dose))
plot.abc.agreement(pstat5.abc$cdf.day1, pstat5.abc$cdf.day2, 'pstat5-abc-cdf.png', main=dose)
}

par(mfrow=c(2,2))
for (dose in DOSES) {
pstat5.abc <- area.between.curves(pstat5, X=X, doses=c('0U',dose))
plot.abc.agreement(pstat5.abc$dens.day1, pstat5.abc$dens.day2, 'pstat5-abc-dens.png', main=dose)
}

par(mfrow=c(2,2))
for (dose in DOSES) {
pstat5.abc <- nn.mean.diff(l, doses=c('0U',dose))
plot.abc.agreement(pstat5.abc$mean.nn.day1, pstat5.abc$mean.nn.day2, 'pstat5-abc-dens.png', main=dose)
}





#a
plot(density(l[["day1.CB00165D.0U"]][,'PSTAT5']))
lines(density(l[["day1.CB00165D.1000U"]][,'PSTAT5']))
plot(density(l[["day2.CB00165D.0U"]][,'PSTAT5']))
lines(density(l[["day2.CB00165D.1000U"]][,'PSTAT5']))
#h
plot(density(l[["day1.CB01498C.0U"]][,'PSTAT5']))
lines(density(l[["day1.CB01498C.1000U"]][,'PSTAT5']))
plot(density(l[["day2.CB01498C.0U"]][,'PSTAT5']))
lines(density(l[["day2.CB01498C.1000U"]][,'PSTAT5']))
#g
plot(density(l[["day1.CB01495Z.0U"]][,'PSTAT5']))
lines(density(l[["day1.CB01495Z.1000U"]][,'PSTAT5']))
plot(density(l[["day2.CB01495Z.0U"]][,'PSTAT5']))
lines(density(l[["day2.CB01495Z.1000U"]][,'PSTAT5']))



print(dim(read.csv('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CLR/CB00010K_2012-11-13.clr')->g))
#load('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/CB01495Z_2012-10-09.RData')
load('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/CB00010K_2012-11-13.RData')

X <- sapply(paste('PSTAT5',1:4,sep='.'), function(x) fcs.data[,x]-fcs.data[,'PSTAT5.1'])

channels <- c("FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3")
Y <- t(apply(X[,paste('PSTAT5',1:4,sep='.')], 1, function(x) coefficients(line(0:3, x))))


plot(NULL, xlim=c(0,3), ylim=range(X[,paste('PSTAT5',1:4,sep='.')]))
apply(X[1:20,paste('PSTAT5',1:4,sep='.')], 1, function(x) lines(1:4,log10(x), lwd=.5, alpha=.5))

apply(X[1:20,paste('PSTAT5',1:4,sep='.')], 1, function(x) abline(lm(x ~ 0:3)))


apply( X[which(Y[,2] >= t.99), pstat5], 1, function(x) lines(0:3, x, lwd=.5) )

t.99 <- quantile(Y[,2],probs=seq(0,1,.0001))[['99.99%']]
table(Y[,2] >= t.99)
pstat5 <- paste('PSTAT5',1:4,sep='.')



plotClusters( fcs.data[,channels], posteriors=cbind(as.numeric(g[,grep('Naive', colnames(g))]),as.numeric(Y[,2] >= t.99)) )
plotClusters( fcs.data[,channels], classification=as.numeric(Y[,2] >= t.99) )

t.99 <- quantile(X[,2],probs=seq(0,1,.0001))[['99.99%']]
plotClusters( fcs.data[,channels], classification=as.numeric(X[,2] >= t.99) )


#
load('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/CB00010K_2012-11-13.RData')
X <- sapply(paste('PSTAT5',1:4,sep='.'), function(x) fcs.data[,x]-fcs.data[,'PSTAT5.1'])
t.99 <- quantile(X[,2],probs=seq(0,1,.0001))[['99.99%']]
plotClusters( fcs.data[,channels], classification=as.numeric(X[,2] >= t.99) )

#
load('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/CB01495Z_2012-10-09.RData')
X <- sapply(paste('PSTAT5',1:4,sep='.'), function(x) fcs.data[,x]-fcs.data[,'PSTAT5.1'])
plot(NULL, xlim=c(0,3), ylim=range(X[,paste('PSTAT5',1:4,sep='.')]))
t.99 <- quantile(X[,2],probs=seq(0,1,.0001))[['99.99%']]
table(as.numeric(X[,2] >= t.99))
plotClusters( fcs.data[,channels], classification=as.numeric(X[,2] >= t.99) )

#load('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CLR/pstat5-join/All/CB00165D_2012-11-29.RData')

#fcs.data <- read.FCS('/chiswick/data/store/facs/Tony/Tony/281112_IL2_sens_T1D_TC17_pSTAT5_Treg/CB00165D_I022596K_CB00165D_0U.fcs', channels=c('SSCH','SSCW',"FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3"), TRANS=logicleTransform(w=1))
fcs.data <- read.FCS('/chiswick/data/store/facs/Tony-FCS/PSTAT5/CD25/CD45RA/CD4/FOXP3/CB00165D_0U_2012-11-29.fcs', channels=c('SSCH','SSCW',"FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3", 'PSTAT5'), TRANS=logicleTransform(w=1))
#fcs.data <- read.FCS('/chiswick/data/store/facs/Tony-FCS/PSTAT5/CD25/CD45RA/CD4/FOXP3/CB00165D_0U_2012-11-29.fcs', channels=c('SSCH','SSCW',"FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3", 'PSTAT5'), TRANS=logicleTransform(w=.1))
#fcs.data <- read.FCS('/chiswick/data/store/facs/Tony-FCS/PSTAT5/CD25/CD45RA/CD4/FOXP3/CB00165D_0U_2012-11-29.fcs', channels=c('SSCH','SSCW',"FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3", 'PSTAT5'), TRANS=logicleTransform())
original.fcs.data <- fcs.data




e <- classification.to.ellipse(fcs.data[,channels],clr[,1]) 
fcs.data <- read.FCS('/chiswick/data/store/facs/Tony-FCS/PSTAT5/CD25/CD45RA/CD4/FOXP3/CB01504J_1000U_2013-03-27.fcs',
                     channels=c('SSCH','SSCW',"FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3", 'PSTAT5'), TRANS=logicleTransform(w=1))
x <- fcs.data


plot(density(fcs.data[which(as.logical(clr[,6])),'PSTAT5']))
plot(density(fcs.data[,'PSTAT5']))

individual <- 'CB00165D'
date <- '2012-11-29'
individual <- 'CB00366X'
date <- '2012-11-07'
l <- lapply(gsub('\\.','',DOSES), function(dose) read.FCS(sprintf('/chiswick/data/store/facs/Tony-FCS/PSTAT5/CD25/CD45RA/CD4/FOXP3/%s_%s_%s.fcs',individual,dose,date),
                                                          channels=c('FSCW','SSCW','FSCH','SSCH',"FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3", 'PSTAT5'),
                                                          TRANS=logicleTransform(w=1)))
names(l) <- DOSES
CLR <- read.csv(sprintf('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CLR/%s_0U_%s.clr',individual,date))

#pdf('~nikolas/dunwich/Projects/IL2/Paper/IL2stim/article/figures/dose-effect.pdf',width=2*5,height=3*5)
pdf('~nikolas/GoogleDrive/PhD/Thesis/IL2/figures/dose-effect.pdf',width=2*5,height=3*5)
par(mfrow=c(3,2), cex.lab=3)
figure.labels <- iter(paste(letters,')',sep=''))
for (chan in c("FSCA", "SSCA", "CD4", "CD45RA", "CD25", "FOXP3")) {
    xquant <- quantile(l[['0U']][,chan],probs=seq(0,1,.01))
    print(xlim <- c(xquant[['1%']],xquant[['99%']]))
    plot(density(l[['0U']][,chan]), xlab=chan, ylab='', main='', col='white', xlim=xlim, cex.lab=3)
    title(nextElem(figure.labels), adj=0)
    i <- 1
    for (dose in DOSES) {
        lines(density(l[[as.character(dose)]][,chan]), lwd=i, col=blues4[i])
        i <- i+1
    }
    legend('topright',DOSES, lwd=1:4, col=blues4[1:4])
}
dev.off()


zoom <- 5
pdf('~nikolas/dunwich/Projects/IL2/Paper/IL2stim/article/figures/dose-effect-pstat5.pdf', width=2*zoom, height=1*zoom)
par(mfrow=c(1,2), cex.lab=1.5)
figure.labels <- iter(paste(letters,')',sep=''))
chan <- 'PSTAT5'
xquant <- quantile(l[['0U']][,chan],probs=seq(0,1,.005))
print(xlim <- c(xquant[['0.5%']],xquant[['99.5%']]))
plot(normalised.density(l[['0U']][,chan]), xlab='pSTAT5', ylab='', main='', col='white', xlim=xlim)
title(nextElem(figure.labels), adj=0)
i <- 1
for (dose in DOSES) {
    lines(normalised.density(l[[as.character(dose)]][,chan]), lwd=i, col=blues4[i])
    i <- i+1
}
x.step <- seq(min(xlim), max(xlim), .01)
plot(x.step, ecdf(l[['0U']][,chan])(x.step), xlab='pSTAT5',  ylab='', main='', col='white', xlim=xlim)
title(nextElem(figure.labels), adj=0)
i <- 1
for (dose in DOSES) {
    lines(x.step, ecdf(l[[as.character(dose)]][,chan])(x.step), lwd=i, col=blues4[i])
    i <- i+1
} 
dev.off()

zoom <- 5
pdf('~nikolas/GoogleDrive/PhD/Thesis/IL2//figures/pstat5-baseline-relative.pdf', width=2*zoom, height=1*zoom)
par(mfrow=c(1,2), cex.lab=1.5)
figure.labels <- iter(paste(letters,')',sep=''))
#xquant <- quantile(fcs.data[,'PSTAT5.1'],probs=seq(0,1,.005))
#print(xlim <- c(xquant[['0.5%']],xquant[['99.5%']]))
plot(normalised.density(fcs.data[,'PSTAT5.1']), xlab='pSTAT5', ylab='', main='', col='white', xlim=c(.5,3))
title(nextElem(figure.labels), adj=0)
i <- 1
for (chan in paste('PSTAT5',1:4,sep='.')) {
    lines(normalised.density(fcs.data[,chan]), lwd=i, col=blues4[i])
    i <- i+1
}
legend('topright',DOSES, lwd=1:4, col=blues4[1:4])
plot(normalised.density(fcs.data[,'diff.PSTAT5.2']), xlab='pSTAT5', ylab='', main='', col='white', xlim=c(-.5,2))
title(nextElem(figure.labels), adj=0)
i <- 1
for (chan in paste('diff','PSTAT5',2:4,sep='.')) {
    lines(normalised.density(fcs.data[,chan]), lwd=i, col=blues4[i])
    i <- i+1
}
legend('topright',DOSES, lwd=1:4, col=blues4[1:4])
dev.off()



# make sure right gate is loaded from file
pdf('~nikolas/dunwich/Projects/IL2/Paper/IL2stim/article/figures/dose-effect-pstat5-cellsubsets-density.pdf')
par(mfrow=c(2,2))
figure.labels <- iter(paste(letters,')',sep=''))
for (gate in names(e)) {
    X <- l[['0U']]
    #X <- X[which(filter.mahalanobis.ellipses(X,e[[gate]])),]
    X <- X[which(filter.mahalanobis.ellipses(X,e[[gate]])),]
    plot(normalised.density(l[['0U']][,'PSTAT5']), xlab='PSTAT5', ylab='', main=gate, col='white', cex=3, xlim=c(0.5,3))
    title(nextElem(figure.labels), adj=0)
    i <- 1
    for (dose in DOSES) {
        X <- l[[dose]]
        X <- X[which(filter.mahalanobis.ellipses(X,e[[gate]])),]
        lines(normalised.density(X[,'PSTAT5']), lwd=i, col=blues4[i])
        i <- i+1
    }
    legend('topright',DOSES, lwd=1:4, col=blues4[1:4])
}
dev.off()


#CLR <- read.csv('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CLR/CB00366X_2012-11-07.clr')
#colnames(CLR) <- CLR.CELL.TYPES
#load( '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/CB00366X_2012-11-07.RData' )
# pSTAT5 density across doses
pdf('~nikolas/dunwich/Projects/IL2/Paper/IL2stim/article/figures/dose-effect-pstat5-cellsubsets-density.pdf')
par(mfrow=c(2,2))
figure.labels <- iter(paste(letters,')',sep=''))
for (gate in CELL.TYPES) {
    X <- fcs.data[which(as.logical(CLR[,gate])),]
    print(y.max <- max(sapply(1:4, function(i) normalised.density(X[,paste('PSTAT5',i,sep='.')])$y)))
    plot(normalised.density(X[,'PSTAT5.1']), ylim=c(0,y.max), xlab='pSTAT5', ylab='', main=gate, col='white', cex=3, xlim=c(0.5,3))
    title(nextElem(figure.labels), adj=0)
    sapply( 1:4, function(i) lines(normalised.density(X[,paste('PSTAT5',i,sep='.')]), lwd=i, col=blues4[i]) )
    legend('topright',DOSES, lwd=1:4, col=blues4[1:4])
}
dev.off()

# baseline relative pSTAT5 density across doses
pdf('~nikolas/dunwich/Projects/IL2/Paper/IL2stim/article/figures/dose-effect-pstat5-cellsubsets-density-baseline.pdf')
par(mfrow=c(2,2))
figure.labels <- iter(paste(letters,')',sep=''))
for (gate in CELL.TYPES) {
    X <- fcs.data[which(as.logical(CLR[,gate])),]
    pstat5.baseline <- X[,'PSTAT5.1']
    pstat5.diff <- sapply(1:4, function(i) X[,paste('PSTAT5',i,sep='.')]-pstat5.baseline)
    colnames(pstat5.diff) <- paste('PSTAT5',1:4,sep='.')
    X <- pstat5.diff
    print(y.max <- max(sapply(1:4, function(i) normalised.density(X[,paste('PSTAT5',i,sep='.')])$y)))
    plot(normalised.density(X[,'PSTAT5.1']), ylim=c(0,y.max), xlab='pSTAT5', ylab='', main=gate, col='white', cex=3, xlim=c(-.1,2))
    title(nextElem(figure.labels), adj=0)
    sapply( 1:4, function(i) lines(normalised.density(X[,paste('PSTAT5',i,sep='.')]), lwd=i, col=blues4[i]) )
    legend('topright',DOSES, lwd=1:4, col=blues4[1:4])
}
dev.off()





pdf('~nikolas/dunwich/Projects/IL2/Paper/IL2stim/article/figures/dose-effect-pstat5-cellsubsets-ecdf.pdf')
par(mfrow=c(2,2))
figure.labels <- iter(paste(letters,')',sep=''))
pstat5.diff <- list()
for (gate in names(e)) {
    X <- l[['0U']]
    X <- X[which(filter.mahalanobis.ellipses(X,e[[gate]])),]
    #x.steps <- seq(min(l[['0U']][,'PSTAT5']),max(l[['0U']][,'PSTAT5']),.1)
    x.steps <- seq(.5, 3, .05)
    plot(x.steps, ecdf(X[,'PSTAT5'])(x.steps), xlab='PSTAT5', ylab='', main=gate, col='white', cex=3 , xlim=c(0.5,3))
    title(nextElem(figure.labels), adj=0)
    i <- 1
    X <- l[['0U']]
    X <- X[which(filter.mahalanobis.ellipses(X,e[[gate]])),]
    baseline.pstat5.ecdf <- ecdf(X[,'PSTAT5'])(x.steps)
    for (dose in DOSES) {
        X <- l[[dose]]
        X <- X[which(filter.mahalanobis.ellipses(X,e[[gate]])),]
        pstat5.ecdf <- ecdf(X[,'PSTAT5'])
        lines(x.steps, pstat5.ecdf(x.steps), lwd=i, col=blues4[i])
        pstat5.diff[[gate]] <- c( pstat5.diff[[gate]], sum(baseline.pstat5.ecdf - pstat5.ecdf(x.steps)))
        i <- i+1
    }
} 
ind.pstat5.diff <- do.call('rbind',pstat5.diff)
dev.off()

load('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/CB00165D_2012-11-29.RData')
Y <- fcs.data

pdf('~nikolas/dunwich/Projects/IL2/Paper/IL2stim/article/figures/dose-effect-pstat5-cellsubsets-density-joined.pdf')
par(mfrow=c(2,2))
figure.labels <- iter(paste(letters,')',sep=''))
for (gate in names(e)) {
    print(gate)
    print(dim(X <- Y[which(filter.mahalanobis.ellipses(Y,e[[gate]])),]))
    plot(normalised.density(X[,'PSTAT5.1']), xlab='PSTAT5', ylab='', main=gate, col='white', xlim=c(0.5,3))
    title(nextElem(figure.labels), adj=0)
    for (i in 1:4) {
        lines(normalised.density(X[,paste('PSTAT5',i,sep='.')]), lwd=i, col=blues4[i])
    }
}
dev.off()

pdf('~nikolas/dunwich/Projects/IL2/Paper/IL2stim/article/figures/dose-effect-pstat5-cellsubsets-ecdf-joined.pdf')
par(mfrow=c(2,2))
figure.labels <- iter(paste(letters,')',sep=''))
pstat5.diff <- list()
for (gate in names(e)) {
    print(gate)
    print(dim(X <- Y[which(filter.mahalanobis.ellipses(Y,e[[gate]])),]))
    baseline.pstat5.ecdf <- ecdf(X[,'PSTAT5.1'])(x.steps)
    plot(x.steps,baseline.pstat5.ecdf, xlab='PSTAT5', ylab='', main=gate, col='white', xlim=c(0.5,3))
    title(nextElem(figure.labels), adj=0)
    for (i in 1:4) {
        pstat5.ecdf <- ecdf(X[,paste('PSTAT5',i,sep='.')])
        lines(x.steps,pstat5.ecdf(x.steps), lwd=i, col=blues4[i])
        pstat5.diff[[gate]] <- c( pstat5.diff[[gate]], sum(baseline.pstat5.ecdf - pstat5.ecdf(x.steps)))
    }
}
joined.pstat5.diff <- do.call('rbind',pstat5.diff)
dev.off()


pdf('~nikolas/dunwich/Projects/IL2/Paper/IL2stim/article/figures/pstat5-response-cellsubsets-ecdf-diff.pdf')
plot(NULL, xlim=c(1,4), ylim=range(c(ind.pstat5.diff,joined.pstat5.diff)), xlab='Dose', ylab='ECDF diff with baseline')
i <- 1
for (gate in names(e)) {
    lines(1:4, joined.pstat5.diff[gate,], col=i, lty=2)
    lines(1:4, ind.pstat5.diff[gate,], col=i)
    i <- i+1
}
legend('topleft',names(e), text.col=1:4)
dev.off()



fcs.data <- read.FCS('/chiswick/data/store/facs/Tony-FCS/PSTAT5/CD25/CD45RA/CD4/FOXP3/CB00396E_0U_2013-03-11.fcs', channels=c('SSCH','SSCW',"FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3", 'PSTAT5'), TRANS=logicleTransform(w=1))

X <- fcs.data[,c('FSCA','SSCA')]
d <- kde2D(X)

lymph.dens <- e.lymphocytes$tau*dmvnorm(X, mean=e.lymphocytes$Mu, sigma=e.lymphocytes$Sigma)
d[,'lymph.dens'] <- lymph.dens

new.e.lymphocytes <- e.lymphocytes
#new.e.lymphocytes$Mu['SSCA'] <- sum(e.lymphocytes$Mu['SSCA'] * lymph.dens) / sum(lymph.dens)
#new.e.lymphocytes$Mu['FSCA'] <- sum(e.lymphocytes$Mu['FSCA'] * lymph.dens) / sum(lymph.dens)
new.e.lymphocytes$Mu <- colMeans(X[which(filter.mahalanobis.ellipse(X, e.lymphocytes, p=.999)),])

f(X, c('FSCA','SSCA'), main='', plot.gates=list(e.lymphocytes, new.e.lymphocytes))

centers <- rbind(e.lymphocytes$Mu, c(4.25,4.25), c(3.25,3.25))
res <- kmeans(X, centers=centers)

X <- read.FCS('/chiswick/data/store/facs/Tony-FCS/PSTAT5/CD25/CD45RA/CD4/FOXP3/CB00366X_0U_2012-11-07.fcs', channels=c('SSCH','SSCW',"FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3", 'PSTAT5'), TRANS=logicleTransform(w=1))
X <- X[filter.mahalanobis.ellipses(X,g.cd4),]

f(X, c('CD45RA','SSCA'), main='', plot.gates=list(e.memory,e.naive))

# these don't work :(
res <- mvnormalmixEM(X[,c('CD45RA','SSCA')], mu=list(e.memory$Mu, e.naive$Mu), sigma=list(e.memory$Sigma, e.naive$Sigma), maxit=10) 
res <-  Mclust( X[,c('CD45RA','SSCA')], G=2, prior=priorControl(mean=rbind(e.memory$Mu,e.naive$Mu))) 
# this works but covariance changes
res <- kmeans(X[,c('CD45RA','SSCA')], centers=rbind(e.memory$Mu,e.naive$Mu)) 
ellipses <- list(e.memory, e.naive)
new.ellipses <- kmeans.magnetic.ellipses(X, ellipses)
f(X, c('CD45RA','SSCA'), main='', plot.gates=c(ellipses, new.ellipses))

### compute % of pSTAT5 positive
pct.pstat5.positive <- function(fcs.data) {
    thresh <- quantile(fcs.data[,'PSTAT5.1'], probs=.99)
    pct <- 100*sapply(1:4, function(i) length(which(fcs.data[,paste('PSTAT5',i,sep='.')]>thresh))/length(fcs.data[,'PSTAT5.1']))
    names(pct) <- paste('pct','PSTAT5',1:4,sep='.')
    return(pct)
}

filter.strictly.increasing <- function(fcs.data) apply(t(apply(fcs.data, 1, diff)),1,function(x)all(x>0))

### per cell-type response obtained from gating
gates.pstat5.response <- function(gates.dir) {
    pstat5.response <- data.frame()
    for (f in list.files(path='~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/', pattern='.*.RData', full.names=TRUE) ) {
        print(f)
        print(load(f))
        individual <- unlist(strsplit(gsub('.RData','',basename(f)), '_'))[[1]]
        date <- unlist(strsplit(gsub('.RData','',basename(f)), '_'))[[2]]
        pstat5.baseline <- fcs.data[,'PSTAT5.1']
        pstat5.diff <- sapply(1:4, function(i) fcs.data[,paste('PSTAT5',i,sep='.')]-pstat5.baseline)
        colnames(pstat5.diff) <- paste('PSTAT5',1:4,sep='.')
        print(load(file.path(sprintf('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/%s/CLR',gates.dir),basename(f))))
        pstat5.response <- rbind(pstat5.response, data.frame(individual=individual, date=date, cell.type='ungated', t(colMeans(pstat5.diff)), t(pct.pstat5.positive(fcs.data))))
        pdf(file.path(sprintf('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/%s/Plots/pSTAT5-density',gates.dir),gsub('.RData','.pdf',basename(f))))
        par(mfrow=c(2,2))
        #figure.labels <- iter(paste(letters,')',sep=''))
        #gate memory eff and tregs
        for (cell.type in c('Memory Eff', 'Memory Treg', 'Naive Eff', 'Naive Treg')) {
        cell.type.filter <- which(as.logical(CLR[,cell.type]))
        #plot the pSTAT5 density at each dose of proleukin
        plot(normalised.density(fcs.data[cell.type.filter,'PSTAT5.1']), lwd=1, col=blues4[1], xlim=c(0,3), main=cell.type)
        #title(nextElem(figure.labels), adj=0)
        sapply(1:4, function(i) lines(normalised.density(fcs.data[cell.type.filter,paste('PSTAT5',i,sep='.')]), lwd=i, col=blues4[i]))
        #pstat5.baseline[,cell.type.filter]
        pstat5.response <- rbind( pstat5.response, data.frame( individual=individual, date=date, cell.type=cell.type,
                                                                                      t(colMeans(pstat5.diff[cell.type.filter,])),
                                                                                      t(pct.pstat5.positive(fcs.data[cell.type.filter,])) ) )
        }
        dev.off()
    }
    return(pstat5.response)
} 

### per cell-type response obtained from fixed gates
fixed.manual.gates.response <- gates.pstat5.response('fixed-manual-gates')
write.csv(fixed.manual.gates.response, file='~nikolas/Projects/IL2/pstat5-response-fixed-manual.csv', row.names=FALSE, quote=FALSE)
### per cell-type response obtained from magnetic gates
magnetic.manual.gates.response <- gates.pstat5.response('magnetic-manual-gates')
write.csv(magnetic.manual.gates.response, file='~nikolas/Projects/IL2/pstat5-response-magnetic-manual.csv', row.names=FALSE, quote=FALSE)
### per cell-type response obtained from magnetic gates with fixed terminal gate
magnetic2.manual.gates.response <- gates.pstat5.response('magnetic-manual-gates2')
write.csv(magnetic2.manual.gates.response, file='~nikolas/Projects/IL2/pstat5-response-magnetic-manual2.csv', row.names=FALSE, quote=FALSE)



for (f in list.files(path='~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/', pattern='.*.RData', full.names=TRUE) ) {
    print(f)
    print(load(f))
    individual <- unlist(strsplit(gsub('.RData','',basename(f)), '_'))[[1]]
    date <- unlist(strsplit(gsub('.RData','',basename(f)), '_'))[[2]]
    pstat5.baseline <- fcs.data[,'PSTAT5.1']
    pstat5.diff <- sapply(1:4, function(i) fcs.data[,paste('PSTAT5',i,sep='.')]-pstat5.baseline)
    colnames(pstat5.diff) <- paste('PSTAT5',1:4,sep='.')
    print(load(file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR',basename(f))))

    par(mfrow=c(2,2))
    i <- 1
    for (cell.type in CELL.TYPES) {
        plot(NULL, xlim=c(0,3), ylim=c(0,max(pstat5.diff[,paste('PSTAT5',1:4,sep='.')])), xaxt='n', xlab='dose', ylab='pSTAT5', main=cell.type)
        axis(1, at=0:3, labels=DOSES)
        apply(pstat5.diff[which(as.logical(CLR[,cell.type])),], 1, function(x) lines(0:3, x, col=i, lwd=.01))
        i <- i+1
    }

}


head(read.csv('pstat5-response-fixed-manual.csv')->fixed.manual.gates.response)
head(read.csv('pstat5-response-magnetic-manual.csv')->magnetic.manual.gates.response)
head(read.csv('pstat5-response-magnetic-manual2.csv')->magnetic2.manual.gates.response)

# pSTAT5 MFI response per cell-subset
pdf('~nikolas/dunwich/Projects/IL2/Paper/IL2stim/article/figures/pstat5-response-cellsubsets-manualgates.pdf', height=2*zoom, width=2*zoom)
#
par(mfrow=c(2,2))
figure.labels <- iter(paste(letters,')',sep=''))
plot(NULL, xlim=c(0,3), ylim=range(fixed.manual.gates.response[,paste('PSTAT5',1:4,sep='.')]), xaxt='n', xlab='dose', ylab='pSTAT5 MFI')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
for (cell.type in CELL.TYPES) {
    apply(fixed.manual.gates.response[grep(cell.type,fixed.manual.gates.response$cell.type),paste('PSTAT5',1:4,sep='.')], 1, function(x) lines(0:3, x, col=i, lwd=.25))
    i <- i+1
}
legend('topleft', CELL.TYPES, text.col=1:4, bty='n')
plot(NULL, xlim=c(0,3), ylim=range(fixed.manual.gates.response[,paste('PSTAT5',1:4,sep='.')]), xaxt='n', xlab='dose', ylab='pSTAT5 MFI')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
colMedians <- function(x) apply(x, 2, median)
for (cell.type in CELL.TYPES) {
    lines(0:3, colMedians( fixed.manual.gates.response[grep(cell.type,fixed.manual.gates.response$cell.type),paste('PSTAT5',1:4,sep='.')] ), col=i, lwd=2)
    i <- i+1
}
legend('topleft', CELL.TYPES, text.col=1:4, bty='n')
#
plot(NULL, xlim=c(0,3), ylim=range(magnetic.manual.gates.response[,paste('PSTAT5',1:4,sep='.')]), xaxt='n', xlab='dose', ylab='pSTAT5 MFI')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
for (cell.type in CELL.TYPES) {
    apply(magnetic.manual.gates.response[grep(cell.type,magnetic.manual.gates.response$cell.type),paste('PSTAT5',1:4,sep='.')], 1, function(x) lines(0:3, x, col=i, lwd=.25))
    i <- i+1
}
legend('topleft', CELL.TYPES, text.col=1:4, bty='n')
plot(NULL, xlim=c(0,3), ylim=range(magnetic.manual.gates.response[,paste('PSTAT5',1:4,sep='.')]), xaxt='n', xlab='dose', ylab='pSTAT5 MFI')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
colMedians <- function(x) apply(x, 2, median)
for (cell.type in CELL.TYPES) {
    lines(0:3, colMedians( magnetic.manual.gates.response[grep(cell.type,magnetic.manual.gates.response$cell.type),paste('PSTAT5',1:4,sep='.')] ), col=i, lwd=2)
    i <- i+1
}
legend('topleft', CELL.TYPES, text.col=1:4, bty='n')
dev.off()

# pSTAT5 MFI density per cell type 
par(mfrow=c(2,2))
i <- 1
for (cell.type in CELL.TYPES) {
    X <- magnetic2.manual.gates.response[grep(cell.type,magnetic2.manual.gates.response$cell.type),]
    plot(density(X[,'PSTAT5.2']), col='white', xlim=c(-.5,1.2), xlab='pSTAT5', main=cell.type)
    sapply(2:4, function(j) lines(density(X[,paste('PSTAT5',j,sep='.')]), col=i, lwd=j))
    i <- i+1
}




# repeatability of pSTAT5 MFI per cell subtype per gating approach
for (cell.type in c('Memory Eff', 'Memory Treg', 'Naive Eff', 'Naive Treg')) {
rep.magnetic<-merge(REPEATS, magnetic.manual.gates.response[grep(cell.type, magnetic.manual.gates.response$cell.type),], all.x=TRUE)
#rep.magnetic<-merge(REPEATS, magnetic2.manual.gates.response[grep(cell.type, magnetic.manual.gates.response$cell.type),], all.x=TRUE)
rep.fixed<-merge(REPEATS, fixed.manual.gates.response[grep(cell.type, magnetic.manual.gates.response$cell.type),], all.x=TRUE)
#par(oma=c(0,0,2,0))
#pdf(sprintf('~nikolas/dunwich/Projects/IL2/Paper/IL2stim/article/figures/repeatability-pstat5-mfi-%s.pdf',gsub(' ','-',cell.type)),width=10,height=10)
pdf(sprintf('~nikolas/GoogleDrive/PhD/Thesis/IL2/figures/repeatability-pstat5-mfi-%s.pdf',gsub(' ','-',cell.type)),width=10,height=10)
par(mfrow=c(2,2))
for (i in 2:4) {
pstat5 <- paste('PSTAT5',i,sep='.')
xlim <- range(c(rep.magnetic[,pstat5],rep.fixed[,pstat5]))
plot(rep.fixed[c(TRUE,FALSE),pstat5],rep.fixed[c(FALSE,TRUE),pstat5], xlim=xlim, ylim=xlim, pch=rep.fixed[c(TRUE,FALSE),'pch'],
     xlab=paste('day 1:', 'pSTAT5 MFI at dose',DOSES[[i]]), ylab=paste('day 2:','pSTAT5 MFI at dose',DOSES[[i]]), col='black')
points( rep.magnetic[c(TRUE,FALSE),pstat5], rep.magnetic[c(FALSE,TRUE),pstat5], pch=rep.magnetic[c(TRUE,FALSE),'pch'], col='red')
rp <- vector('expression',2)
rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(rep.fixed[c(TRUE,FALSE),pstat5],rep.fixed[c(FALSE,TRUE),pstat5])**2,3))**2,3))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(rep.magnetic[c(TRUE,FALSE),pstat5],rep.magnetic[c(FALSE,TRUE),pstat5])**2,3))**2,3))[2]
legend('topleft', legend=rp, text.col=c('black','red'), bty='n')
abline(b=1, a=0)
}
pstat5 <- paste('PSTAT5',1:4,sep='.')
fixed.day1 <- rowSums(rep.fixed[c(TRUE,FALSE),pstat5])
fixed.day2 <- rowSums(rep.fixed[c(FALSE,TRUE),pstat5])
magnetic.day1 <- rowSums(rep.magnetic[c(TRUE,FALSE),pstat5])
magnetic.day2 <- rowSums(rep.magnetic[c(FALSE,TRUE),pstat5])
xlim <- range(c(fixed.day1, fixed.day2, magnetic.day1, magnetic.day2))
plot(fixed.day1, fixed.day2, xlim=xlim, ylim=xlim, pch=rep.fixed[c(TRUE,FALSE),'pch'], xlab='day 1: pSTAT5 AUC', ylab='day 2: pSTAT5 AUC', col='black')
points( magnetic.day1, magnetic.day2, pch=rep.magnetic[c(TRUE,FALSE),'pch'], col='red')
rp <- vector('expression',2)
rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(fixed.day1,fixed.day2)**2,3))**2,3))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(magnetic.day1,magnetic.day2)**2,3))**2,3))[2]
legend('topleft', legend=rp, text.col=c('black','red'), bty='n')
abline(b=1, a=0)
#title(cell.type, outer=TRUE)
dev.off()
}



# repeatability of % pSTAT5 positive per cell subtype per gating approach
for (cell.type in c('Memory Eff', 'Memory Treg', 'Naive Eff', 'Naive Treg')) {
rep.magnetic<-merge(REPEATS, magnetic.manual.gates.response[grep(cell.type, magnetic.manual.gates.response$cell.type),], all.x=TRUE)
#rep.magnetic<-merge(REPEATS, magnetic2.manual.gates.response[grep(cell.type, magnetic.manual.gates.response$cell.type),], all.x=TRUE)
rep.fixed<-merge(REPEATS, fixed.manual.gates.response[grep(cell.type, magnetic.manual.gates.response$cell.type),], all.x=TRUE)
#par(oma=c(0,0,2,0))
pdf(sprintf('~nikolas/dunwich/Projects/IL2/Paper/IL2stim/article/figures/repeatability-pstat5-percent-positive-%s.pdf',gsub(' ','-',cell.type)),width=10,height=10)
par(mfrow=c(2,2))
for (i in 2:4) {
#pstat5 <- paste('PSTAT5',i,sep='.')
pstat5 <- paste('pct', 'PSTAT5',i,sep='.')
xlim <- range(c(rep.magnetic[,pstat5],rep.fixed[,pstat5]))
plot(rep.fixed[c(TRUE,FALSE),pstat5],rep.fixed[c(FALSE,TRUE),pstat5], xlim=xlim, ylim=xlim, pch=rep.fixed[c(TRUE,FALSE),'pch'],
     xlab=paste('day 1:', 'pSTAT5 MFI at dose',DOSES[[i]]), ylab=paste('day 2:','pSTAT5 MFI at dose',DOSES[[i]]), col='black')
points( rep.magnetic[c(TRUE,FALSE),pstat5], rep.magnetic[c(FALSE,TRUE),pstat5], pch=rep.magnetic[c(TRUE,FALSE),'pch'], col='red')
rp <- vector('expression',2)
rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(rep.fixed[c(TRUE,FALSE),pstat5],rep.fixed[c(FALSE,TRUE),pstat5])**2,3))**2,3))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(rep.magnetic[c(TRUE,FALSE),pstat5],rep.magnetic[c(FALSE,TRUE),pstat5])**2,3))**2,3))[2]
legend('topleft', legend=rp, text.col=c('black','red'), bty='n')
abline(b=1, a=0)
}
pstat5 <- paste('PSTAT5',1:4,sep='.')
fixed.day1 <- rowSums(rep.fixed[c(TRUE,FALSE),pstat5])
fixed.day2 <- rowSums(rep.fixed[c(FALSE,TRUE),pstat5])
magnetic.day1 <- rowSums(rep.magnetic[c(TRUE,FALSE),pstat5])
magnetic.day2 <- rowSums(rep.magnetic[c(FALSE,TRUE),pstat5])
xlim <- range(c(fixed.day1, fixed.day2, magnetic.day1, magnetic.day2))
plot(fixed.day1, fixed.day2, xlim=xlim, ylim=xlim, pch=rep.fixed[c(TRUE,FALSE),'pch'], xlab='day 1: pSTAT5 AUC', ylab='day 2: pSTAT5 AUC', col='black')
points( magnetic.day1, magnetic.day2, pch=rep.magnetic[c(TRUE,FALSE),'pch'], col='red')
rp <- vector('expression',2)
rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(fixed.day1,fixed.day2)**2,3))**2,3))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(magnetic.day1,magnetic.day2)**2,3))**2,3))[2]
legend('topleft', legend=rp, text.col=c('black','red'), bty='n')
abline(b=1, a=0)
#title(cell.type, outer=TRUE)
dev.off()
}
saturated <- function(x) any(x[2:3] > (1:2)*x[4]/3)

table( as.numeric(apply(pstat5.diff, 1, function(x) saturated(x))) )


### lymphocytes are the most responsive cluster
zoom <- 5
pdf('~nikolas/Thesis/figures/pstat5-response-decision-tree.pdf', height=1*zoom, width=3*zoom)
figure.labels <- iter(paste(letters,')',sep=''))
par(mfrow=c(1,3))
print(load('~nikolas/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/FCS/pstat5-join/CB00406Q_2012-06-12.RData'))
fcs.data <- baseline.relative.pstat5(fcs.data)
#load('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/transforms.RData')
load('~/dunwich/Projects/IL2/transforms.RData')
#for (i in grep('PSTAT5',names(transforms))) transforms[[i]]<-logicleTransform(w=.5)
fcs.data <- applyTransforms(fcs.data,transforms)
#t <- tree::tree( diff.PSTAT5.4 ~ FSCA + SSCA, data=data.frame(fcs.data) )
t <- tree::tree( diff.PSTAT5.4 ~ FSCA + SSCA, data=data.frame(fcs.data) )
#t <- tree::tree( diff.PSTAT5.4 ~ FSCA + SSCA + CD4 + CD3 + CD25 + CD8 + CD56 + FOXP3 + CD45RA, data=data.frame(fcs.data) )
t <- prune.tree( t, best=3 )
w <- as.factor(t$where)
levels(w) <- 1:length(unique(w))
# a)
plot(t)
title(nextElem(figure.labels), adj=0)
text(t)
# b)
y.max <- max(sapply(sort(unique(w)), function(i) normalised.density(fcs.data[which(w==i),'diff.PSTAT5.4'])$y))
plot(normalised.density(fcs.data[,'diff.PSTAT5.4']), main='pSTAT5 response at 1000U', ylim=c(0,y.max))
title(nextElem(figure.labels), adj=0)
sapply(sort(unique(w)), function(i) lines(normalised.density(fcs.data[which(w==i),'diff.PSTAT5.4']),col=i,lwd=2))
# c)
smoothPlot( fcs.data[,c('FSCA','SSCA')], classification=w, ellipses=FALSE, chull.lwd=4 )
title(nextElem(figure.labels), adj=0)
dev.off()

### two different samples yield very different decision trees
zoom <- 7
pdf('~nikolas/Thesis/figures/two-sample-decision-tree.pdf', height=1*zoom, width=2*zoom)
figure.labels <- iter(paste(letters,')',sep=''))
par(mfrow=c(1,3))
load('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/transforms.RData')
for (i in grep('PSTAT5',names(transforms))) transforms[[i]]<-logicleTransform(w=.5)
print(load('~nikolas/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/FCS/pstat5-join/CB00406Q_2012-06-12.RData'))
fcs.data <- applyTransforms(fcs.data,transforms)
fcs.data <- baseline.relative.pstat5(fcs.data)
fcs.data1 <- fcs.data
print(load('~nikolas/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/FCS/pstat5-join/CB00086S_2012-09-18.RData'))
fcs.data <- applyTransforms(fcs.data,transforms)
fcs.data <- baseline.relative.pstat5(fcs.data)
fcs.data2 <- fcs.data
par(mfrow=c(1,2))
t <- tree::tree( diff.PSTAT5.4 ~ FSCA + SSCA + CD4 + CD3 + CD8 + CD56 + CD45RA + CD25 + FOXP3, data=data.frame(fcs.data1) )
# a)
plot(t)
title(nextElem(figure.labels), adj=0)
text(t)
fcs.data <- applyTransforms(fcs.data,transforms)
t <- tree::tree( diff.PSTAT5.4 ~ FSCA + SSCA + CD4 + CD3 + CD8 + CD56 + CD45RA + CD25 + FOXP3, data=data.frame(fcs.data2) )
# b)
plot(t)
title(nextElem(figure.labels), adj=0)
text(t)
dev.off()

# transformed vs untransformed data
library(tree)
load('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/transforms.RData')
for (i in grep('PSTAT5',names(transforms))) transforms[[i]]<-logicleTransform(w=.5)
print(load('~nikolas/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/FCS/pstat5-join/CB00406Q_2012-06-12.RData'))
fcs.data <- baseline.relative.pstat5(fcs.data)
fcs.data2 <- applyTransforms(fcs.data,transforms)
par(mfrow=c(1,2))
t <- tree::tree( diff.PSTAT5.4 ~ FSCA + SSCA + CD4 + CD3 + CD8 + CD56 + CD45RA + CD25 + FOXP3, data=data.frame(fcs.data) )
plot(t)
text(t)
fcs.data <- applyTransforms(fcs.data,transforms)
t <- tree::tree( diff.PSTAT5.4 ~ FSCA + SSCA + CD4 + CD3 + CD8 + CD56 + CD45RA + CD25 + FOXP3, data=data.frame(fcs.data2) )
plot(t)
text(t)



zoom <- 5
pdf('~nikolas/Thesis/figures/pstat5-response-decision-tree.pdf', height=1*zoom, width=3*zoom)
figure.labels <- iter(paste(letters,')',sep=''))
par(mfrow=c(1,3))
# a)
plot(t)
title(nextElem(figure.labels), adj=0)
text(t)
# b)
y.max <- max(sapply(sort(unique(w)), function(i) normalised.density(fcs.data[which(w==i),'diff.PSTAT5.4'])$y))
plot(normalised.density(fcs.data[,'diff.PSTAT5.4']), main='pSTAT5 response at 1000U', ylim=c(0,y.max))
title(nextElem(figure.labels), adj=0)
sapply(sort(unique(w)), function(i) lines(normalised.density(fcs.data[which(w==i),'diff.PSTAT5.4']),col=i,lwd=2))
# c)
smoothPlot( fcs.data[,c('FSCA','SSCA')], classification=w, ellipses=FALSE, chull.lwd=4 )
title(nextElem(figure.labels), adj=0)
dev.off()




# rpart on CD4 lymphocytes
cd4.lymph <- as.data.frame(fcs.data[which(as.logical(CLR[,'CD4'])),])
pdf('~nikolas/dunwich/Projects/IL2/Paper/IL2stim/article/figures/pstat5-response-tree-cd4lymphocytes.pdf')
figure.labels <- iter(paste(letters,')',sep=''))
#par(mfrow = c(3, 3), mai=c(0,0,0,0), oma=c(2,3,2,.5))
par(mfrow = c(3, 3))
#PSTAT5.4
t <- tree::tree( diff.PSTAT5.4 ~ CD45RA + CD25 + FOXP3, data=cd4.lymph )
t <- prune.tree( t, best=4 )
w <- as.factor(t$where)
levels(w) <- 1:length(unique(w)) 
# a)
plot(t)
title(nextElem(figure.labels), adj=0)
text(t)
# b)
plot(normalised.density(cd4.lymph[,'diff.PSTAT5.4']), main='pSTAT5 response at 1000U', xlim=c(-.5,2))
title(nextElem(figure.labels), adj=0)
sapply(sort(unique(w)), function(i) lines(normalised.density(cd4.lymph[which(w==i),'diff.PSTAT5.4']),col=i,lwd=2))
# c)
smoothPlot( cd4.lymph[,c('CD45RA','CD25')], classification=w, ellipses=FALSE, chull.lwd=2, outliers=TRUE )
title(nextElem(figure.labels), adj=0)
#PSTAT5.3
t <- tree::tree( diff.PSTAT5.3 ~ CD45RA + CD25 + FOXP3, data=cd4.lymph )
t <- prune.tree( t, best=4 )
w <- as.factor(t$where)
levels(w) <- 1:length(unique(w))
# d)
plot(t)
title(nextElem(figure.labels), adj=0)
text(t)
# e)
plot(normalised.density(cd4.lymph[,'diff.PSTAT5.3']), main='pSTAT5 response at 10U', xlim=c(-.5,2))
title(nextElem(figure.labels), adj=0)
sapply(sort(unique(w)), function(i) lines(normalised.density(cd4.lymph[which(w==i),'diff.PSTAT5.3']),col=i,lwd=2))
# f)
smoothPlot( cd4.lymph[,c('CD45RA','CD25')], classification=w, ellipses=FALSE, chull.lwd=2, outliers=TRUE )
title(nextElem(figure.labels), adj=0)
#PSTAT5.2
t <- tree::tree( diff.PSTAT5.2 ~ CD45RA + CD25 + FOXP3, data=cd4.lymph )
t <- prune.tree( t, best=4 )
w <- as.factor(t$where)
levels(w) <- 1:length(unique(w))
# g)
plot(t)
title(nextElem(figure.labels), adj=0)
text(t)
# h)
plot(normalised.density(cd4.lymph[,'diff.PSTAT5.2']), main='pSTAT5 response at 0.1U', xlim=c(-.5,2))
title(nextElem(figure.labels), adj=0)
sapply(sort(unique(w)), function(i) lines(normalised.density(cd4.lymph[which(w==i),'diff.PSTAT5.2']),col=i,lwd=2))
# i)
smoothPlot( cd4.lymph[,c('FOXP3','CD25')], classification=w, ellipses=FALSE, chull.lwd=2, outliers=TRUE )
title(nextElem(figure.labels), adj=0)
dev.off()


pdf('~nikolas/dunwich/Projects/IL2/Paper/IL2stim/article/figures/pstat5-response-manual-gates-cd4lymphocytes.pdf')
par(mfrow=c(1,1))
plotClusters(cd4.lymph[,c('CD45RA','CD25','FOXP3')], posteriors=CLR[which(CLR[,'CD4']==1),c('Memory Eff','Memory Treg','Naive Eff','Naive Treg')], chulls=FALSE, outliers=TRUE, lower=smoothPlot)
dev.off()


plotClusters(d[,c('SSCA','FSCA','CD4','PSTAT5.4')], classification=t$where)

fcs.data <- cbind(fcs.data,baseline.relative.pstat5(fcs.data))

f <- function(fcs.data) {
    t <- tree::tree( diff.PSTAT5.4 ~ SSCW + FSCW + SSCH + FSCH + SSCA + FSCA + CD4 + CD25 + FOXP3 + CD45RA, data.frame(fcs.data) )
    par(mfrow=c(2,5))
    for (i in unique(t$where)) {
        plot(normalised.density(fcs.data[,'diff.PSTAT5.4']), main=i, lty=2, xlim=c(0,3))
        sapply(1:4, function(j) lines(normalised.density(fcs.data[which(i==t$where),paste('diff.PSTAT5',j,sep='.')]),col=i,lwd=j))
    }
    return(t)
}

t <- f(fcs.data)
peak.shifts <- do.call('rbind', by(fcs.data, t$where, function(x) colMedians(x[,grep('diff',colnames(x))])))
max.peak.shift <- as.numeric(rownames(peak.shifts[order(peak.shifts[,'diff.PSTAT5.4'],decreasing=TRUE),]))[1]
plotClusters(fcs.data[,c('FSCA','SSCA','CD4','CD25','CD45RA','FOXP3')], classification=as.numeric(t$where==max.peak.shift))

plot(t)
text(t)

library(rpart)
library(partykit)

rp.lymph <- rpart( diff.PSTAT5.4 ~ SSCW + FSCW + SSCH + FSCH + SSCA + FSCA + CD4 + CD25 + FOXP3 + CD45RA, data.frame(fcs.data[which(as.logical(clr[,2])),]), method='anova', model=TRUE )
plot(as.party(rp.lymph))

rp.cd4 <- rpart( diff.PSTAT5.4 ~ SSCW + FSCW + SSCH + FSCH + SSCA + FSCA + CD4 + CD25 + FOXP3 + CD45RA, data.frame(fcs.data[which(as.logical(clr[,3])),]), method='anova', model=TRUE )
plot(as.party(rp.cd4))

rp.naive <- rpart( diff.PSTAT5.4 ~ SSCW + FSCW + SSCH + FSCH + SSCA + FSCA + CD4 + CD25 + FOXP3 + CD45RA, data.frame(fcs.data[which(as.logical(clr[,7])),]), method='anova', model=TRUE )
rp.naive <- rpart( diff.PSTAT5.4 ~ CD25 + FOXP3 , data.frame(fcs.data[which(as.logical(clr[,7])),]), method='anova', model=TRUE )
plot(as.party(rp.naive))
plotClusters(data.frame(fcs.data[which(as.logical(clr[,7])),c('FOXP3','CD25','CD45RA','diff.PSTAT5.4')]), classification=rp.naive$where, outliers=TRUE)

rp.memory <- rpart( diff.PSTAT5.4 ~ CD25 + FOXP3 + CD45RA , data.frame(fcs.data[which(as.logical(clr[,4])),]), method='anova', model=TRUE )
plot(as.party(rp.memory))
plotClusters(data.frame(fcs.data[which(as.logical(clr[,4])),c('FOXP3','CD25','CD45RA','diff.PSTAT5.4')]), classification=rp.memory$where, outliers=TRUE)

# only works as a classification tree
#library(oblique.tree)
#ot.naive <- oblique.tree( formula=diff.PSTAT5.4 ~ SSCW + FSCW + SSCH + FSCH + SSCA + FSCA + CD4 + CD25 + FOXP3 + CD45RA, data=data.frame(fcs.data[which(as.logical(clr[,7])),]), oblique.splits='on' )





plotClusters(fcs.data[,c('FSCA','SSCA','CD4','CD25','CD45RA','FOXP3')], classification=as.numeric(rp$where==17))

t <- f(d)
plotClusters(d[,c('FSCA','SSCA','CD4','CD25','CD45RA','FOXP3')], classification=as.numeric(t$where==7))


lines(normalised.density(d[,'PSTAT5.3']))

t <- tree::tree( PSTAT5.4 ~ CD25 + FOXP3, data.frame(d) )
plotClusters(cbind(d[,c('CD25','FOXP3')],pstat5.diff[,4]), classification=t$where)


library(randomForest)
t <- randomForest( PSTAT5.4 ~ CD25 + FOXP3, data.frame(d) )

plotClusters(cbind(d[,c('CD25','FOXP3')],pstat5.diff[,4]), classification=t$where)

#
source('~nikolas/bin/FCS/fcs.R')
setwd('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/')

file.name <- 'CB01340F_2012-11-21.RData'
file.name <- 'CB01510Q_2012-11-21.RData'
setwd('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/')
load(sprintf('pstat5-join/All/%s',file.name))
load(sprintf('magnetic-manual-gates/CLR/%s',file.name))

fcs.data[which(as.logical(CLR[,'CD4'])),]->d
baseline.relative.pstat5(d)->d

channels <- c("FSCW", "SSCW", "FSCH", "SSCH", "FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3")
form <- formula(sprintf('diff.PSTAT5.4 ~ %s', paste(channels, collapse=' + ')))

t <- tree(form, data.frame(d),split='deviance')
t <- prune.tree(t,best=3)

f <- as.data.frame(t$frame)
plotClusters(d[,c('diff.PSTAT5.4', as.character(unique(f$var[which(f$splits[,'cutleft']!='')])))], classification=t$where, ellipses=FALSE, chull.lwd=2)

plotClusters(d[,c('diff.PSTAT5.4','CD45RA','FOXP3','CD25')], classification=t$where, ellipses=FALSE, chull.lwd=2)

plotClusters(d[,c('diff.PSTAT5.4','CD45RA','FOXP3','CD25')],classification=prune.tree(t,best=3)$where, ellipses=FALSE, chull.lwd=2)

fcs.data <- as.data.frame(fcs.data)
fcs.data <- baseline.relative.pstat5(fcs.data)

m <- lm(form, fcs.data)
coeffs <- data.frame(t(coefficients(m)))
for (cell.subset in CLR.CELL.TYPES) {
    print(cell.subset)
    m <- lm(form, fcs.data[which(as.logical(CLR[,cell.subset])),])
    coeffs <- rbind(coeffs, data.frame(t(coefficients(m))))
} 
rownames(coeffs) <- c('All', CLR.CELL.TYPES)

plot(normalised.density(fcs.data[which(!as.logical(CLR[,'CD4'])),]$diff.PSTAT5.4),col='red')
lines(normalised.density(fcs.data$diff.PSTAT5.4))


lines(normalised.density(fcs.data[which(as.logical(CLR[,'CD4'])),]$diff.PSTAT5.4),col='blue')
lines(normalised.density(fcs.data[which(as.logical(CLR[,'Memory Eff'])),]$diff.PSTAT5.4),col='red', lwd=2)


plotClusters(fcs.data[,c('FSCA','SSCA','CD4','CD25','CD45RA')], classification=as.numeric(!as.logical(CLR[,'CD4']) & fcs.data$diff.PSTAT5.4 > 2))


d <- fcs.data[which(!as.logical(CLR[,'CD4'])),]


t <- tree::tree(form, d)
f <- as.data.frame(t$frame)
plotClusters(d[,c('diff.PSTAT5.4',as.character(unique(f$var[which(f$splits[,'cutleft']!='')])))], classification=t$where, ellipses=FALSE, chull.lwd=2)
#
t <- tree::tree(formula(sprintf('diff.PSTAT5.2 ~ %s', paste(channels, collapse=' + '))), d)
f <- as.data.frame(t$frame)
plotClusters(d[,c('diff.PSTAT5.2', as.character(unique(f$var[which(f$splits[,'cutleft']!='')])))], classification=t$where, ellipses=FALSE, chull.lwd=2)




t <- tree::tree(CD4 ~ diff.PSTAT5.4 + diff.PSTAT5.3 + diff.PSTAT5.2, fcs.data)
plotClusters(fcs.data[,c('CD4','diff.PSTAT5.4','diff.PSTAT5.3','diff.PSTAT5.2')], classification=t$where, ellipses=FALSE, chull.lwd=2)
plotClusters(fcs.data[,c('diff.PSTAT5.4','diff.PSTAT5.3','diff.PSTAT5.2')])

X <- fcs.data[,c('diff.PSTAT5.4','diff.PSTAT5.3','diff.PSTAT5.2')]
library(mclust)
res <- Mclust(X, G=2, modelNames=c('VVV'), initialization=list(subset=1:1000))
plotClusters(X, classification=res$classification)

table(res$classification)


plotClusters(fcs.data[,c('SSCA','FSCA','CD4')])
plotClusters(fcs.data[which(res$classification==1),c('SSCA','FSCA','CD4')])
plotClusters(fcs.data[which(res$classification==2),c('SSCA','FSCA','CD4')])


d <- fcs.data[which(res$classification==2),c('SSCA','FSCA','CD4')]

X <- X[which(res$classification==2),]
res <- Mclust(X, G=2, modelNames=c('VVV'), initialization=list(subset=1:1000))
plotClusters(X, classification=res$classification)

X <- X[which(res$classification==2),]
res <- Mclust(X, G=2, modelNames=c('VVV'), initialization=list(subset=1:1000))
plotClusters(X, classification=res$classification)

X <- X[which(res$classification==1),]
res <- Mclust(X, G=2, modelNames=c('VVV'), initialization=list(subset=1:1000))
plotClusters(X, classification=res$classification)


# divide cells as low responders and high responders on pSTAT5 response (i.e baseline subtracted) at 1000U 
# within responders further divide on low/high on pSTAT5 response at 10U
# repeat on pSTAT5 response at 0.1U
# cluster ones which are consistently high, these are the most sensitive cell populations

library(mixtools)

d <- fcs.data


for (i in 4:2) {
diff.pstat5 <- paste('diff','PSTAT5',i,sep='.')
x <- d[,diff.pstat5]
res <- kmeans(x, centers=c(0,1))
mu <- as.numeric(res$centers)
lambda <- as.numeric(prop.table(table(res$cluster)))
sigsqrd <- as.numeric(tapply(x,res$cluster,var))
m <- mixtools::normalmixEM2comp(x,mu=mu,lambda=lambda,sigsqrd=sigsqrd,verb=TRUE)
d <- d[which( (d[,diff.pstat5] > 0) & (m$posterior[,2]>.99)),]
#plot(m,which=2)
#lines(density(d[,diff.pstat5]))
}


# divide cells as low responders and high responders on pSTAT5 response (i.e baseline subtracted) at 1000U 
# within responders further divide on low/high on pSTAT5 response at 10U
# repeat on pSTAT5 response at 0.1U
# cluster ones which are consistently high, these are the most sensitive cell populations 

split.response <- function(d, i) {
    if (i < 2) return(list(d=d))
    if (nrow(d)<10) return(list(d=d))
    print(diff.pstat5 <- paste('diff','PSTAT5',i,sep='.'))
    x <- d[,diff.pstat5]
    res <- kmeans(x, centers=c(0,1))
    mu <- as.numeric(res$centers)
    lambda <- as.numeric(prop.table(table(res$cluster)))
    sigsqrd <- as.numeric(tapply(x,res$cluster,var))
    m <- mixtools::normalmixEM2comp(x,mu=mu,lambda=lambda,sigsqrd=sigsqrd)
    i1 <- which( (d[,diff.pstat5] > 0) & (m$posterior[,1]>.99))
    d1 <- d[i1,]
    i2 <- which( (d[,diff.pstat5] > 0) & (m$posterior[,2]>.99))
    d2 <- d[i2,]
    return(list(i1=i1, d1=split.response(d1, i-1), i2=i2, d2=split.response(d2, i-1), d=d, m=m))
}

X <- split.response(fcs.data[which(!as.logical(CLR[,'CD4'])),], 4) 
X <- split.response(fcs.data[which(as.logical(CLR[,'CD4'])),], 4)
X <- split.response(fcs.data, 4)

hi.x <- X$d2$d2$d2$d
hc <- hclust(dist(hi.x[,c('FSCA','SSCA','CD4','CD45RA','FOXP3','CD25')]))
pairs(hi.x[,c('FSCA','SSCA','CD4','CD45RA','FOXP3','CD25')], col=cutree(hc,k=8), pch=20)
points(hi.x[,c('SSCH','SSCW')], col=cutree(hc,k=6))
plotClusters(fcs.data[,c('FSCA','SSCA','CD4','CD45RA','FOXP3','CD25')], outliers=TRUE, plot.points=hi.x, classification=cutree(hc,k=8), upper=contourPlot, ellipses=FALSE, chulls=FALSE)

b <- fcs.data[X$d2$i2,][X$d2$d2$i2,][X$d2$d2$d2$i2,]

## the recursive partitioning of pSTAT5 starts at the highest dose then moves
## down to lower doses to find the subset which are the first responders
zoom <- 2
pdf('~nikolas/dunwich/Projects/IL2/Paper/IL2stim/article/figures/pstat5-rpart.pdf',width=5*zoom,height=3*zoom)
par(mfrow=c(3,5))
frame()
frame()
plot(X$m,which=2,main2='',xlab2='pSTAT5 response at 1000U')
frame()
frame()
#
frame()
plot(X$d1$m,which=2,main2='low',xlab2='pSTAT5 response at 10U')
frame()
plot(X$d2$m,which=2,main2='high',xlab2='pSTAT5 response at 10U')
frame()
#
plot(X$d1$d1$m,which=2,main2='low',xlab2='pSTAT5 response at 0.1U')
plot(X$d1$d2$m,which=2,main2='high',xlab2='pSTAT5 response at 0.1U')
frame()
plot(X$d2$d1$m,which=2,main2='low',xlab2='pSTAT5 response at 0.1U')
plot(X$d2$d2$m,which=2,main2='high',xlab2='pSTAT5 response at 0.1U')
dev.off()

plotClusters(fcs.data[,c('FSCA','SSCA','CD4','CD45RA','FOXP3','CD25')], outliers=TRUE, plot.points=X$d2$d2$d2$d)

#nn <- nn2(fcs.data,query=hi.x,k=1)

first.resp <- data.frame()
for (f in list.files('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/first-response-cells',pattern='.*.RData', full.names=TRUE)) {
    print(load(f))
    print(nrow(hi.x))
    #first.resp <- rbind(first.resp, cbind(hi.x,gsub('.RData','',basename(f))))
    first.resp <- rbind(first.resp, hi.x)
}

d <- first.resp[,-grep('PSTAT5',colnames(first.resp))]
library(fastcluster)

hc <- hclust(dist(d))

princomp(d)->pc


plotClusters(first.resp[,c('FSCA','SSCA','CD4','CD45RA','FOXP3','CD25')])

rpart::rpart(diff.PSTAT5.2 ~ FSCA + SSCA + CD4 + CD25 + CD45RA + FOXP3, data=first.resp)$where->w


plot(density(d$diff.PSTAT5.4))

plot(density(x <- d$diff.PSTAT5.3))
res <- kmeans(x, centers=c(0,1))
mu <- as.numeric(res$centers)
lambda <- as.numeric(prop.table(table(res$cluster)))
sigsqrd <- as.numeric(tapply(x,res$cluster,var))
m <- mixtools::normalmixEM2comp(x,mu=mu,lambda=lambda,sigsqrd=sigsqrd,verb=TRUE)
plot(m,which=2)
d <- d[which( (x > 0) & (m$posterior[,2]>.99) ),]

plot(density(d$diff.PSTAT5.2))
lines(density(fcs.data$diff.PSTAT5.2))

plot(density(x <- d$diff.PSTAT5.2))
res <- kmeans(x, centers=c(0,1))
mu <- as.numeric(res$centers)
lambda <- as.numeric(prop.table(table(res$cluster)))
sigsqrd <- as.numeric(tapply(x,res$cluster,var))
m <- mixtools::normalmixEM2comp(x,mu=mu,lambda=lambda,sigsqrd=sigsqrd,verb=TRUE)
plot(m,which=2)
d <- d[which( (x > 0) & (m$posterior[,2]>.99) ),]


plot(density(x <- d$diff.PSTAT5.2))
res <- kmeans(x, centers=c(0,1))
mu <- as.numeric(res$centers)
lambda <- as.numeric(prop.table(table(res$cluster)))
sigsqrd <- as.numeric(tapply(x,res$cluster,var))
m <- mixtools::normalmixEM2comp(x,mu=mu,lambda=lambda,sigsqrd=sigsqrd,verb=TRUE)
plot(m,which=2)

plot(density(fcs.data[,'CD4']))
lines(density(d[i<-which(m$posterior[,2]>.99),'CD4']))

plot(density(fcs.data[which(as.logical(CLR[,'CD4'])),'CD45RA']))
lines(density(d[i<-which(m$posterior[,2]>.99),'CD45RA']))

plot(density(fcs.data[which(as.logical(CLR[,'Memory'])),'CD25']))
lines(density(d[i<-which(m$posterior[,2]>.99),'CD25']), col='red')

plot(density(d[i,'diff.PSTAT5.2']),xlim=c(-.5,3))
lines(density(d[i,'diff.PSTAT5.3']))
lines(density(d[i,'diff.PSTAT5.4']))


plot(density(d$diff.PSTAT5.3), xlim=c(-.5,3))
lines(density(d$diff.PSTAT5.3))
lines(density(d$diff.PSTAT5.2))


ip <- read.csv('~nikolas/Projects/IL2/IL2RA-genotype/ip-IL2RA-2013-12-23.csv')
d <- data.frame()
individuals <- c()
gate.dir <- sprintf('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates/CLR/')
for (f in list.files(pattern='.*.RData')) {
    load(f)
    individual <- (gsub('_.*.RData','',f))
    t1d <- ip[grep(individual,as.character(ip$uniqueID)),'t1d']-1
    cat(individual, t1d,'\n')
    load(file.path(gate.dir, f))
    d <- rbind(d,cbind(fcs.data[which(as.logical(CLR[,'Naive'])),],t1d=t1d))
    individuals <- c(individual, individuals)
}
d$t1d <- as.factor(d$t1d)
d <- baseline.relative.pstat5(d)
form <- formula(paste('t1d','~', paste(grep('t1d',colnames(d),value=T,invert=T),collapse=' + ')))
library(tree)
library(rpart)
t <- rpart::rpart(form, data=d)
library(partykit)
plot(as.party(t))


X <- data.frame(do.call('rbind',strsplit(gsub('.RData','',list.files(pattern='.*.RData')), '_')))
colnames(X) <- c('individual', 'date')
colnames(ip) <- gsub('uniqueID','individual',colnames(ip))
X <- (merge(X,ip))

#control
i1 <- 'CB00396E_2012-09-25.RData'
load(i1)
f1 <- fcs.data
load(file.path(gate.dir, i1))
g1 <- CLR

#case
i2 <- 'CB01484M_2012-09-25.RData'
load(i2)
f2 <- fcs.data
load(file.path(gate.dir, i2))
g2 <- CLR

mean(f1[which(as.logical(g1[,'Memory'])),'CD25'])
mean(f2[which(as.logical(g2[,'Memory'])),'CD25'])


library(randomForest)
cd4.lymph <- as.data.frame(fcs.data[which(as.logical(CLR[,'CD4'])),])
form <- formula(paste('diff.PSTAT5.4','~',paste(grep('PSTAT5',colnames(fcs.data),value=TRUE,invert=TRUE),collapse='+')))
RF <- randomForest(form, data=cd4.lymph)
dim(cd4.lymph)
head(getTree(RF,1))

library(rpart)
RPART <- rpart::rpart(form,data=cd4.lymph) 

library(party)
CTREE <- party::ctree(form,data=cd4.lymph)




####
load('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/CB00010K_2012-11-13.RData')
fcs.data <- fcs.data[percentile.filter(fcs.data),]
fcs.data <- baseline.relative.pstat5(fcs.data)

fcs.data <- fcs.data[,-grep('^PSTAT5', colnames(fcs.data))]
fcs.data <- fcs.data[,-grep('diff.PSTAT5.1',colnames(fcs.data))]

stats <- data.frame()
for (k in 2:10) {
    system.time(res <- kmeans(fcs.data, centers=k))
    t.var <- apply(fcs.data,2,var)
    w.var <- apply(fcs.data,2,function(x) tapply(x, res$cluster, var))
    stats <- rbind(stats, t(colSums(w.var)/t.var))
    print(k)
    print(sum(sort(colSums(w.var)/t.var)))
}

plot(NULL, xlim=c(1,11), ylim=range(stats))
apply(stats, 2, function(x) lines(2:10,x))



par(mfrow=c(3,5)) 
it<-iter(colnames(fcs.data))
apply(fcs.data,2,function(x) {
      plot(normalised.density((x)),main=nextElem(it))
      for (i in unique(sort(res$cluster))) lines(normalised.density(x[which(res$cluster==i)]),col=i, lwd=2)
     })



d <- normalised.density(fcs.data, from=min(fcs.data[,'CD4']), to=max(fcs.data[,'CD4']))
i <- 1
d1 <- normalised.density(fcs.data[which(res$cluster==i),'CD4'], from=min(fcs.data[,'CD4']), to=max(fcs.data[,'CD4']))
i <- 2
d2 <- normalised.density(fcs.data[which(res$cluster==i),'CD4'], from=min(fcs.data[,'CD4']), to=max(fcs.data[,'CD4']))

chan <- 'CD4'
dens <- sapply( sort(unique(res$cluster)), function(i) normalised.density(fcs.data[which(res$cluster==i),chan], from=min(fcs.data[,chan]), to=max(fcs.data[,chan]))$y )

# KL distance (not a symmetric function)
KL <- function(y1,y2) {
    kl <- log(y1/y2)*y1
    kl[which(!is.finite(kl))] <- 0
    sum(kl)
}


####
source('~nikolas/bin/FCS/hdr.R')
X <- fcs.data[,-grep('PSTAT5',colnames(fcs.data))]
system.time(dens <- knn.dens(X,k=1000))

library(FNN)
fnn.dist <- FNN::knn.dist(x,k)

d <- fcs.data

for (i in 1:100) {
    print(i)
    idx <- FNN::knnx.index(d, query=t(d[sample(1:nrow(d),size=1),]), k=13200)
    d <- d[-idx,]
}

#distance matrix
distance.matrix <- dist(d)

#fast hclust single linkage
library(fastcluster)
hc <- fastcluster::hclust(distance.matrix, method='single')

#minimum spanning tree
library(vegan)
st <- spantree(distance.matrix)

library(nnclust)
system.time(st <- mst(d))

st <- mst(fcs.data)


#pstat5-peak-normalisation


## new sensitive cell subset?

f <- 'CB01496A_2012-11-19.RData'
f <- 'CB00366X_2012-11-07.RData'
f <- 'CB00366X_2012-11-07.RData'

panel <- 'CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5'
f <- 'CB00406Q_2012-06-12.RData'
f <- 'CB00573X_2012-06-12.RData'

panel <- 'PSTAT5-CD25-CD45RA-CD4-FOXP3' 
f <- 'CB01494Y_2012-10-09.RData'


print(length(files <- list.files(pattern='.*.RData')))
X<-do.call('rbind',lapply(files,function(f){print(load(f));return(hi.x)}))
dim(X)
smoothPlot(X[,c('CD4','CD25')], outliers=TRUE)
res <- mclust::Mclust(X[,c('CD4','CD25')], G=6, modelNames='VVV', initialization=list(subset=1:2000))
plotClustRes(X[,c('CD4','CD25')], res, outliers=TRUE)
xx <- hi.x[which(predict( res, hi.x[,c('CD4','CD25')] )$classification==6),]

load('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/transform-w1.RData')
#
#pdf('~nikolas/GoogleDrive/PhD/Thesis/IL2/figures/new-cell-subset.pdf',width=2*5,height=2*5)
#pdf('~nikolas/Thesis/figures/new-cell-subset.pdf',width=2*5,height=2*5)
#pdf('~nikolas/Thesis/figures/new-cell-subset2.pdf',width=2*5,height=2*5)
par(mfrow=c(2,2))
setwd(file.path('~/dunwich/Projects/IL2/',panel,'first-response-cells')) 
print(length(files <- list.files(pattern='.*.RData')))
X<-do.call('rbind',lapply(files,function(f){print(load(f));return(hi.x)}))
figure.labels <- iter(paste(letters,')',sep=''))
# pooled samples
smoothPlot(X[,c('FSCA','SSCA')])
#smoothPlot(X[,c('CD4','CD25')])
title(nextElem(figure.labels), adj=0)
#G <- locator(type='l',col='purple',lwd=4) 
lines(G$x,G$y,col='purple',lwd=4)
print(load(f))
#in.poly <- point.in.polygon( hi.x[,'CD4'] , hi.x[,'CD25'], G$x, G$y)
in.poly <- point.in.polygon( hi.x[,'FSCA'] , hi.x[,'SSCA'], G$x, G$y)
xx<-hi.x[as.logical(in.poly),]
# single sample
#plot(hi.x[,c('CD4','CD25')],pch=20)
plot(hi.x[,c('FSCA','SSCA')],pch=20)
title(nextElem(figure.labels), adj=0)
#points(xx[,c('CD4','CD25')],col='purple',pch=20)
points(xx[,c('FSCA','SSCA')],col='purple',pch=20)
# not lymphocytes
print(load(file.path('~/dunwich/Projects/IL2/',panel,'RData','pstat5-join',basename(f))))
fcs.data <- applyTransforms(fcs.data, transforms)
fcs.data <- baseline.relative.pstat5(fcs.data)
print(load(file.path('~/dunwich/Projects/IL2/',panel,'CLR',basename(f))))
smoothPlot(fcs.data[,c('FSCA','SSCA')], outliers=FALSE)
points(xx[,c('FSCA','SSCA')], col='purple', pch=20)
title(nextElem(figure.labels), adj=0)
#plot.ellipses(list(e.lymphocytes)) 
# single sample dose response
ylim <- range(sapply(CELL.TYPES, function(cell.type) colMedians(fcs.data[as.logical(CLR[,cell.type]),paste('diff.PSTAT5',1:4,sep='.')]) ))
plot(NULL, xlim=c(0,3), ylim=ylim, xaxt='n', xlab='dose', ylab='pSTAT5 MFI')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
for (cell.type in CELL.TYPES) {
    mfi <- fcs.data[which(as.logical(CLR[,cell.type])),paste('diff.PSTAT5',1:4,sep='.')]
    lines(0:3, colMedians( mfi ), col=i, lwd=3)
    #polygon(
            #c(0:3,3:0),
            #c( colQuantile( mfi, prob=.75 ),
            #rev(colQuantile( mfi, prob=.25)) ),
            #col=do.call('rgb', c(as.list(t(col2rgb(i)/255)),alpha=.1)),
            #border=i,
            #lwd=.5)
    i <- i+1
} 
a <- xx[,paste('diff.PSTAT5',1:4,sep='.')]
lines(0:3, colMedians( a ), col='purple', lwd=2)
polygon(
        c(0:3,3:0),
        c( colQuantile( a, prob=.75 ),
        rev(colQuantile( a, prob=.25)) ),
        col=do.call('rgb', c(as.list(t(col2rgb('purple')/255)),alpha=.1)),
        border='purple',
        lwd=.5)
legend('topleft', CELL.TYPES, text.col=1:4, bty='n')
#dev.off()


smoothPlot(X[,c('FSCA','SSCA')])
G <- locator(type='l',col='purple',lwd=4) 
in.poly <- point.in.polygon(X[,'FSCA'],X[,'SSCA'],G$x,G$y)
X2 <- X[as.logical(in.poly),]
smoothPlot(X2[,c('CD4','CD25')])
G2 <- locator(type='l',col='purple',lwd=4) 
in.poly <- point.in.polygon(X2[,'CD4'],X2[,'CD25'],G2$x,G2$y)
X3 <- X2[as.logical(in.poly),]



channels <- c('CD3','CD56')
channels <- c('SSCH','SSCW')
channels <- c('CD45RA','CD8')
channels <- c('CD45RA','CD3')
channels <- c('CD8','CD3')
smoothPlot(fcs.data[,channels])
points(xx[,channels], col='purple', pch=20)



points(applyTransforms(hi.x[,c('CD4','CD25')], transforms),pch=20,col='purple')
G <- locator(type='l') 
in.poly <- point.in.polygon( transforms[['CD4']](hi.x[,'CD4']), transforms[['CD25']](hi.x[,'CD25']), G$x, G$y)

smoothPlot(fcs.data[,c('FSCA','SSCA')])
points(hi.x[as.logical(in.poly),c('FSCA','SSCA')], pch=20, col='purple')


d <- read.FCS('/chiswick/data/store/facs/Tony-FCS/CD14-CD19-CD3-CD45RA-CD56-HLA-DR-PSTAT5/CB01496A_0U_2013-06-07.fcs', channels=c('FSCA','SSCA','CD14','CD19','CD3','CD45RA','CD56','HLADR'),TRANS=logicleTransform(w=1))

#merge on c('FSCA','SSCA','CD45RA',)

# CD56 CD8 CD3
d2 <- read.FCS('/chiswick/data/store/facs/Tony-FCS/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/CB01496A_0U_2013-09-10.fcs', channels=c('CD25','CD3','CD4','CD45RA','CD56','CD8','FOXP3','PSTAT5'), TRANS=logicleTransform(w=1))

d2 <- read.FCS('/chiswick/data/store/facs/Tony-FCS/CD25-CD3-CD4-CD45RA-CD56-CD8-PSTAT5/CB01496A_0U_2013-09-10.fcs', TRANS=logicleTransform(w=1), channels=c('CD25','CD3','CD4','CD45RA','CD56','CD8','PSTAT5'))

   
# does not exist
#d2 <- read.FCS('/chiswick/data/store/facs/Tony-FCS/CD25-CD3-CD4-CD45RA-CD56-FOXP3-PSTAT5/CB01496A_0U_2013-09-10.fcs', channels=c('FSCA','SSCA','CD25','CD3','CD4','CD45RA','CD56','CD8','FOXP3','PSTAT5'), TRANS=logicleTransform(w=1))


print(load('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/CB01496A_2012-11-19.RData'))
X2 <- apply(fcs.data,2,logicleTransform(w=1))
join.chan <- grep('PSTAT5', colnames(X2), invert=TRUE, value=TRUE)
new.markers <- c('CD8','CD56','CD3')
X1 <- read.FCS('/chiswick/data/store/facs/Tony-FCS/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/CB01496A_0U_2013-09-10.fcs', channels=c(join.chan,new.markers), TRANS=logicleTransform(w=1))

nn <- RANN::nn2(X2[,join.chan],query=X1[,join.chan],k=1)

dim(X <- cbind(X1, X2[nn$nn.idx,grep('PSTAT5',colnames(X2))]))

fcs.data <- baseline.relative.pstat5(X)

save(fcs.data,file='~nikolas/x.RData')

