library(plyr)
source('~nikolas/Projects/IL2/bin/common.R')

DOSES <- gsub('\\.', '', DOSES)

#barplot representing number of cases and controls analysed per day
d <- structure(c(4L, 1L, 1L, 2L, 2L, 2L, 1L, 1L, 3L, 2L, 1L, 1L, 1L, 
1L, 2L, 2L, 1L, 1L, 1L, 2L, 1L, 2L, 4L, 1L, 2L, 1L, 1L, 3L, 1L, 
1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 3L, 1L, 2L, 1L, 2L, 2L, 2L, 1L, 
2L, 2L, 1L, 1L, 0L, 1L, 0L, 1L, 2L, 1L, 1L, 2L, 1L, 1L, 2L, 1L, 
1L, 1L, 0L, 2L, 2L, 1L, 2L, 3L, 1L, 1L, 1L, 1L, 1L), .Dim = c(37L, 
2L), .Dimnames = list(c("2012-05-28", "2012-05-29", "2012-06-12", 
"2012-06-25", "2012-06-29", "2012-07-02", "2012-07-23", "2012-07-26", 
"2012-07-31", "2012-08-03", "2012-08-23", "2012-08-30", "2012-09-04", 
"2012-09-11", "2012-09-14", "2012-09-18", "2012-09-20", "2012-09-25", 
"2012-09-26", "2012-10-02", "2012-10-03", "2012-10-05", "2012-10-09", 
"2012-10-16", "2012-10-18", "2012-10-23", "2012-10-25", "2012-11-07", 
"2012-11-13", "2012-11-19", "2012-11-21", "2012-11-29", "2013-01-22", 
"2013-01-29", "2013-03-07", "2013-03-11", "2013-03-27"), c("case", 
"control")))
d <- t(d)
  
pdf('~nikolas/Thesis/figures/IL2-sample-time.pdf',width=10,height=5)
par(las=2)
par(mar=c(8,5,1,1))
barplot(d,col=c('red','darkblue'),legend=rownames(d))
dev.off()

## Which metric to use for pSTAT5 response?
# choose metric which gives lowest within-individual variance?
# area under curve
# percent positive

# not joined
individual <- 'CB01495Z'
pstat5dens <- data.frame()
pch <- REPEATS[which(REPEATS$individual==individual),'pch']
day <- REPEATS[which(REPEATS$individual==individual),'day1']
#
for (gate in CELL.TYPES) {
   print(gate)
   for (dose in DOSES) {
      print(dose)
      print((load(sprintf('/chiswick/data/store/facs/Tony-RData/PSTAT5-CD25-CD45RA-CD4-FOXP3/%s_%s_%s.RData',individual,dose,day))))
      print(dim(fcs.data))
      print(load(file.path(base.dir,'magnetic-manual-gates2','CLR', sprintf('%s.RData',paste(individual, dose, day, sep='_')))))
      print(length(g <- which(as.logical(CLR[,gate]))))
      X <- logicleTransform(w=1)(fcs.data[g,'PSTAT5'])
      y <- normalised.density(X,from=-1,to=3,n=512)$y
      x <- normalised.density(X,from=-1,to=3,n=512)$x
      d <- data.frame(individual=individual,day=day,cell.type=gate,dose=dose,t(y))
      pstat5dens <- rbind(pstat5dens,d)
    }
}
pdf('~nikolas/Thesis/figures/dose-effect-pstat5-cellsubsets-density.pdf')
par(mfrow=c(2,2))
figure.labels <- iter(paste(letters,')',sep=''))
for (gate in CELL.TYPES) {
  X <- pstat5dens[which(pstat5dens$cell.type==gate),]
  x.range <- c(0.5, 3)
  y.range <- range(X[,grep('^X',colnames(X))])
  plot(NULL, xlim=x.range, ylim=y.range, main=gate, xlab='pSTAT5', ylab='')
  title(nextElem(figure.labels), adj=0)
  j <- 0
  for (day in unique(X$day)) {
    print(day)
    x1 <- as.numeric(X[which(X$day==day & X$dose=='0U'),grep('^X',colnames(X))])
    j <- j+1
    abline(v=x[max(which(cumsum(x1)/sum(x1)<.99))],col=j,lwd=.5)
    for (i in 1:4) {
    #for (i in 1:1) {
    #for (i in c(1,4)) {
      print(dose <- DOSES[i])
      lines(x, X[which(X$day==day & X$dose==dose),grep('^X',colnames(X))],lwd=i,col=blues4[i])
    }
  }
  cat('\n')
  legend('topright',DOSES, lwd=1:4, col=blues4[1:4])
}
dev.off()


base.dir <- '~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/RData'
fcs.data.day1 <- file.path(base.dir,'pstat5-join/', sprintf('%s.RData',paste(REPEATS$individual, REPEATS$day1, sep='_')))
fcs.data.day2 <- file.path(base.dir,'pstat5-join/', sprintf('%s.RData',paste(REPEATS$individual, REPEATS$day2, sep='_')))

clr.day1 <- file.path(base.dir,'magnetic-manual-gates2/CLR', sprintf('%s.RData',paste(REPEATS$individual, REPEATS$day1, sep='_')))
clr.day1 <- file.path(base.dir,'magnetic-manual-gates2/CLR', sprintf('%s.RData',paste(REPEATS$individual, REPEATS$day1, sep='_')))
clr.day2 <- file.path(base.dir,'magnetic-manual-gates2/CLR', sprintf('%s.RData',paste(REPEATS$individual, REPEATS$day2, sep='_')))

# joined
for ( individual in REPEATS$individual ) {
    pstat5dens <- data.frame()
    pch <- REPEATS[which(REPEATS$individual==individual),'pch']
    days <- REPEATS[which(REPEATS$individual==individual),c('day1', 'day2')]
    for (day in days) {
      #
      print(load(file.path(base.dir,'pstat5-join/All/', sprintf('%s.RData',paste(individual, day, sep='_')))))
      fcs.data <- baseline.relative.pstat5(fcs.data,REPLACE=FALSE)
      print(load(file.path(base.dir,'magnetic-manual-gates2/CLR', sprintf('%s.RData',paste(individual, day, sep='_')))))
      for (gate in CELL.TYPES) {
        X <- fcs.data[which(as.logical(CLR[,gate])),]
        for (dose in paste('PSTAT5',1:4,sep='.')) {
          y <- normalised.density(logicleTransform(w=1)(X[,dose]),from=-1,to=3,n=512)$y
          x <- normalised.density(logicleTransform(w=1)(X[,dose]),from=-1,to=3,n=512)$x
          d <- data.frame(individual=individual,day=day,cell.type=gate,dose=dose,t(y))
          pstat5dens <- rbind(pstat5dens,d)
        }
      }
    }
    # the density curves look very different from one day to the next
    print(sprintf('~nikolas/Thesis/figures/nn-dose-effect-pstat5-cellsubsets-density-repeatability-%s.pdf',pch))
    pdf(sprintf('~nikolas/Thesis/figures/nn-dose-effect-pstat5-cellsubsets-density-repeatability-%s.pdf',pch))
    par(mfrow=c(2,2))
    figure.labels <- iter(paste(letters,')',sep=''))
    for (gate in CELL.TYPES) {
      X <- pstat5dens[which(pstat5dens$cell.type==gate),]
      x.range <- c(0.5, 3)
      y.range <- range(X[,grep('^X',colnames(X))])
      plot(NULL, xlim=x.range, ylim=y.range, main=gate, xlab='pSTAT5', ylab='')
      title(nextElem(figure.labels), adj=0)
      j <- 0
      for (day in unique(X$day)) {
        print(day)
        x1 <- as.numeric(X[which(X$day==day & X$dose=='PSTAT5.1'),grep('^X',colnames(X))])
        j <- j+1
        abline(v=x[max(which(cumsum(x1)/sum(x1)<.99))],col=j,lwd=.5)
        for (dose in 1:4) {
          print(dose)
          lines(x, X[which(X$day==day & X$dose==paste('PSTAT5',dose,sep='.')),grep('^X',colnames(X))],lwd=dose,col=j)
        }
      }
      #legend('topright',DOSES, lwd=1:4) #, col=blues4[1:4])
    }
    dev.off()
}

# not joined
for ( individual in REPEATS$individual ) {
    pstat5dens <- data.frame()
    pch <- REPEATS[which(REPEATS$individual==individual),'pch']
    days <- REPEATS[which(REPEATS$individual==individual),c('day1', 'day2')]
    for (day in days) {
      #
      for (gate in CELL.TYPES) {
        print(gate)
        for (dose in DOSES) {
          print(dose)
          print((load(sprintf('/chiswick/data/store/facs/Tony-RData/PSTAT5-CD25-CD45RA-CD4-FOXP3/%s_%s_%s.RData',individual,dose,day))))
          print(dim(fcs.data))
          print(load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3','magnetic-manual-gates2/CLR', sprintf('%s.RData',paste(individual, dose, day, sep='_')))))
          print(length(g <- which(as.logical(CLR[,gate]))))
          X <- logicleTransform(w=1)(fcs.data[g,'PSTAT5'])
          y <- normalised.density(X,from=-1,to=3,n=512)$y
          x <- normalised.density(X,from=-1,to=3,n=512)$x
          d <- data.frame(individual=individual,day=day,cell.type=gate,dose=dose,t(y))
          pstat5dens <- rbind(pstat5dens,d)
        }
      }
    }
    # the density curves look very different from one day to the next
    print(sprintf('~nikolas/Thesis/figures/dose-effect-pstat5-cellsubsets-density-repeatability-%s.pdf',pch))
    pdf(sprintf('~nikolas/Thesis/figures/dose-effect-pstat5-cellsubsets-density-repeatability-%s.pdf',pch))
    #print(sprintf('~nikolas/Thesis/figures/resting-pstat5-cellsubsets-density-repeatability-%s.pdf',pch))
    #pdf(sprintf('~nikolas/Thesis/figures/resting-pstat5-cellsubsets-density-repeatability-%s.pdf',pch))
    par(mfrow=c(2,2))
    figure.labels <- iter(paste(letters,')',sep=''))
    for (gate in CELL.TYPES) {
      X <- pstat5dens[which(pstat5dens$cell.type==gate),]
      x.range <- c(0.5, 3)
      y.range <- range(X[,grep('^X',colnames(X))])
      plot(NULL, xlim=x.range, ylim=y.range, main=gate, xlab='pSTAT5', ylab='')
      title(nextElem(figure.labels), adj=0)
      j <- 0
      for (day in unique(X$day)) {
        print(day)
        x1 <- as.numeric(X[which(X$day==day & X$dose=='0U'),grep('^X',colnames(X))])
        j <- j+1
        abline(v=x[max(which(cumsum(x1)/sum(x1)<.99))],col=j,lwd=.5)
        #for (i in 1:4) {
        #for (i in 1:1) {
        for (i in c(1,4)) {
          print(dose <- DOSES[i])
          lines(x, X[which(X$day==day & X$dose==dose),grep('^X',colnames(X))],lwd=i,col=j)
        }
      }
      #legend('topright',DOSES, lwd=1:4) #, col=blues4[1:4])
    }
    dev.off()
}


# not joined
pstat5.mfi <- data.frame()
for ( individual in REPEATS$individual) {
    pch <- REPEATS[which(REPEATS$individual==individual),'pch']
    days <- REPEATS[which(REPEATS$individual==individual),c('day1', 'day2')]
    for (day in days) {
      #
      for (gate in CELL.TYPES) {
        print(gate)
        for (dose in gsub('\\.','',DOSES)) {
          print(dose)
          print((load(sprintf('/chiswick/data/store/facs/Tony-RData/PSTAT5-CD25-CD45RA-CD4-FOXP3/%s_%s_%s.RData',individual,dose,day))))
          print(dim(fcs.data))
          print(load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3','magnetic-manual-gates2/CLR', sprintf('%s.RData',paste(individual, dose, day, sep='_')))))
          print(length(g <- which(as.logical(CLR[,gate]))))
          mfi <- median(fcs.data[g,'PSTAT5'])
          d <- data.frame(individual=individual,day=day,cell.type=gate,dose=dose,mfi)
          pstat5.mfi <- rbind(pstat5.mfi,d)
        }
      }
    }
}




pstat5 <- lapply(gsub('\\.','',DOSES), function(dose) {
      print((load(sprintf('/chiswick/data/store/facs/Tony-RData/PSTAT5-CD25-CD45RA-CD4-FOXP3/%s_%s_%s.RData',individual,dose,days[[1]]))))
      print(dim(fcs.data))
      print(load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3','magnetic-manual-gates2/CLR', sprintf('%s.RData',paste(individual, dose, days[[1]], sep='_')))))
      print(length(g <- which(as.logical(CLR[,'Naive Treg']))))
      return(logicleTransform(w=1)(fcs.data[g,'PSTAT5']))
})

load(sprintf('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/%s_%s.RData',individual,days[[1]]))
print(load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3','magnetic-manual-gates2/CLR', sprintf('%s.RData',paste(individual, days[[1]], sep='_')))))
fcs.data <- fcs.data[which(as.logical(CLR[,'Naive Treg'])),]

par(mfrow=c(1,2))
plot(normalised.density(pstat5[[i]]),col='white')
for (i in 1:4) lines(normalised.density(pstat5[[i]]),col=blues4[[i]],lwd=i)
plot(normalised.density(logicleTransform(w=1)(fcs.data[,paste('PSTAT5',4,sep='.')])),col='white')
for (i in 1:4) lines(normalised.density(logicleTransform(w=1)(fcs.data[,paste('PSTAT5',i,sep='.')])),col=blues4[[i]],lwd=i)


# the density curves look very different from one day to the next and correcting for resting sample does not improve things much
pdf('~nikolas/Thesis/figures/dose-effect-pstat5-cellsubsets-density-repeatability.pdf')
par(mfrow=c(2,2))
figure.labels <- iter(paste(letters,')',sep=''))
for (gate in CELL.TYPES) {
  X <- pstat5dens[which(pstat5dens$cell.type==gate),]
  x.range <- c(0.5, 3)
  y.range <- range(X[,grep('^X',colnames(X))])
  plot(NULL, xlim=x.range, ylim=y.range, main=gate, xlab='pSTAT5', ylab='')
  title(nextElem(figure.labels), adj=0)
  j <- 0
  for (day in unique(X$day)) {
    print(day)
    x1 <- as.numeric(X[which(X$day==day & X$dose=='PSTAT5.1'),grep('^X',colnames(X))])
    j <- j+1
    abline(v=x[max(which(cumsum(x1)/sum(x1)<.99))],col=j,lwd=.5)
    for (dose in 1:4) {
      print(dose)
      lines(x, X[which(X$day==day & X$dose==paste('PSTAT5',dose,sep='.')),grep('^X',colnames(X))],lwd=dose,col=j)
    }
  }
  #legend('topright',DOSES, lwd=1:4) #, col=blues4[1:4])
}
dev.off()


#
ddply( pstat5mfi,c('individual','cell.type'), function(x) {
    response.diff <- sum(abs(x[1,grep('PSTAT5',colnames(x))]-x[2,grep('PSTAT5',colnames(x))]))
    return(data.frame(individual=x$individual,cell.type=x$cell.type,response.diff=response.diff)[1,])
}) 
ddply( baseline.pstat5mfi,c('individual','cell.type'), function(x) {
    response.diff <- sum(abs(x[1,grep('PSTAT5',colnames(x))]-x[2,grep('PSTAT5',colnames(x))]))
    return(data.frame(individual=x$individual,cell.type=x$cell.type,response.diff=response.diff)[1,])
})
#
ddply( pstat5pos,c('individual','cell.type'), function(x) {
    response.diff <- sum(abs(x[1,grep('PSTAT5',colnames(x))]-x[2,grep('PSTAT5',colnames(x))]))
    return(data.frame(individual=x$individual,cell.type=x$cell.type,response.diff=response.diff)[1,])
}) 
ddply( baseline.pstat5pos,c('individual','cell.type'), function(x) {
    response.diff <- sum(abs(x[1,grep('PSTAT5',colnames(x))]-x[2,grep('PSTAT5',colnames(x))]))
    return(data.frame(individual=x$individual,cell.type=x$cell.type,response.diff=response.diff)[1,])
})



X <- fcs.data[,c('FSCA','SSCA')]
smoothPlot(X)
G <- locator(type='l')

library(sp)
in.poly <- point.in.polygon( X[,1] , X[,2], G$x, G$y)
e <- classification.to.ellipse(X,in.poly)
plot.ellipses(list(e))

dim(X[in.poly,])


points(X[in.poly,],pch='.')
ch <- chull(X[in.poly,])
lines(X[ch,],col='black')
points(X[ch,],pch='x')


f <- function(fcs.data, channels, main, outliers=FALSE, plot.gates=NULL) {
    xquant <- quantile(fcs.data[,channels[[1]]],probs=seq(0,1,.01))
    yquant <- quantile(fcs.data[,channels[[2]]],probs=seq(0,1,.01))
    print(xlim <- c(xquant[['1%']],xquant[['99%']]))
    print(ylim <- c(yquant[['1%']],yquant[['99%']]))
    if (!outliers)
        smoothScatter( fcs.data[,channels], xlab=channels[[1]], ylab=channels[[2]], main=main, colramp=colorRampPalette(c('white','blue','green','yellow','orange','red')), nrpoints=0, cex.lab=2, cex.main=2, ylim=ylim, xlim=xlim )
    else
        smoothScatter( fcs.data[,channels], xlab=channels[[1]], ylab=channels[[2]], main=main, colramp=colorRampPalette(c('white','blue','green','yellow','orange','red')), nrpoints=0, cex.lab=2, cex.main=2 )
    if (!is.null(plot.gates))
    for (e in plot.gates)
        lines(ellipse::ellipse(e$Sigma,centre=e$Mu),col='red',lwd=2,lty=1)
}


par(mfrow=c(3,2), cex.lab=2, cex.main=2, las=1, mar=c(6,6,2,1))
x <- fcs.data
figure.labels <- iter(paste(letters,')',sep='')) 
# Lymphocytes
channels <- c('FSCA','SSCA')
#e.lymphocytes <- mahalanobis.magnetic.ellipse(x,e.lymphocytes)
f(x, channels, main='Lymphocytes', plot.gates=list(e.lymphocytes))
title(nextElem(figure.labels), adj=0)
lymphocytes.filter <- filter.mahalanobis.ellipse(x,e.lymphocytes)
x <- x[ which(lymphocytes.filter) , ] 
plot.gate.chull(fcs.data[,channels], CLR[,'Lymphocytes'])


# Single cells
channels <- c('SSCH','SSCW')
e.singlecells <- mahalanobis.magnetic.ellipse(x,e.singlecells)
f(x, channels, main='Single cells', outliers=FALSE, plot.gates=list(e.singlecells))
title(nextElem(figure.labels), adj=0)
singlecells.filter <- filter.mahalanobis.ellipse(x,e.singlecells)
x <- x[ which(singlecells.filter) , ] 
plot.gate.chull(fcs.data[,channels], CLR[,'Single Cells'])

# CD4
channels <- c('CD4','SSCA')
e.cd4 <- mahalanobis.magnetic.ellipse(x,e.cd4)
f(x, channels, main='CD4+', outliers=FALSE, plot.gates=list(e.cd4))
title(nextElem(figure.labels), adj=0)
cd4.filter <- filter.mahalanobis.ellipse(x,e.cd4)
x <- x[ which(cd4.filter) , ] 
plot.gate.chull(fcs.data[,channels], CLR[,'CD4'])

# Memory / Naive
channels <- c('CD45RA','SSCA')
new.ellipses <- em.magnetic.ellipses(x, list(e.naive, e.memory))
e.naive <- new.ellipses[[1]]
e.memory <- new.ellipses[[2]]
f(x, channels, main='Memory / Naive', outliers=TRUE, plot.gates=list(e.memory, e.naive))
title(nextElem(figure.labels), adj=0)
naive.filter <- filter.mahalanobis.ellipse(x, e.naive)
memory.filter <- filter.mahalanobis.ellipse(x, e.memory)
x.naive <- x[ which(naive.filter) , ]
x.memory <- x[ which(memory.filter) , ] 
plot.gate.chull(fcs.data[,channels], CLR[,'Memory'])
plot.gate.ellipse(fcs.data[,channels], CLR[,'Memory'])
plot.gate.chull(fcs.data[,channels], CLR[,'Naive'])
plot.gate.ellipse(fcs.data[,channels], CLR[,'Naive'])

# Naive Eff / Treg
channels <- c('CD25','FOXP3')
f(x.naive, channels, main='Naive Eff / TReg', outliers=TRUE)
#lines(ellipse::ellipse(e.naive.conv$Sigma,centre=e.naive.conv$Mu),col='darkgreen',lwd=3)
#lines(ellipse::ellipse(e.naive.tregs$Sigma,centre=e.naive.tregs$Mu),col='blue',lwd=3)
plot.gate.chull(fcs.data[,channels], CLR[,'Naive Eff'])
plot.gate.ellipse(fcs.data[,channels], CLR[,'Naive Eff'])
plot.gate.chull(fcs.data[,channels], CLR[,'Naive Eff'])
plot.gate.ellipse(fcs.data[,channels], CLR[,'Naive Treg'])
title(nextElem(figure.labels), adj=0)
naive.eff.filter <- filter.mahalanobis.ellipse(x.naive, e.naive.conv)
naive.tregs.filter <- filter.mahalanobis.ellipse(x.naive, e.naive.tregs) 


# Memory Eff / Treg
channels <- c('CD25','FOXP3')
f(x.memory, channels, main='Memory Eff / TReg', outliers=TRUE)
#lines(ellipse::ellipse(e.memory.conv$Sigma,centre=e.memory.conv$Mu),col='black',lwd=3)
#lines(ellipse::ellipse(e.memory.tregs$Sigma,centre=e.memory.tregs$Mu),col='red',lwd=3)
plot.gate.chull(fcs.data[,channels], CLR[,'Memory Eff'])
plot.gate.ellipse(fcs.data[,channels], CLR[,'Memory Eff'])
plot.gate.chull(fcs.data[,channels], CLR[,'Memory Eff'])
plot.gate.ellipse(fcs.data[,channels], CLR[,'Memory Treg'])
title(nextElem(figure.labels), adj=0)
memory.eff.filter <- filter.mahalanobis.ellipse(x.memory, e.memory.conv)
memory.tregs.filter <- filter.mahalanobis.ellipse(x.memory, e.memory.tregs)



#flowWorkspace
library(flowWorkspace)
ws<-openWorkspace('/dunwich/scratch/nikolas/Projects/IL2/IL2_s_T1D_Treg_pSTAT5_copy.wsp')
G<-parseWorkspace(ws,name=1,path=ws@path,isNcdf=FALSE,cleanup=FALSE,keep.indices=TRUE)
pop.stats <- lapply(G,getPopStats)

#11 memory eff %pSTAT5+
head(d <- data.frame(id=names(pop.stats), cell.type='Memory Eff', do.call('rbind', lapply(pop.stats, function(x) x[11,]))))
#9 memory Treg %pSTAT5+
d <- rbind(d, data.frame(id=names(pop.stats), cell.type='Memory Treg', do.call('rbind', lapply(pop.stats, function(x) x[9,]))))
#17 naive eff %pSTAT5+
d <- rbind(d, data.frame(id=names(pop.stats), cell.type='Naive Eff', do.call('rbind', lapply(pop.stats, function(x) x[17,]))))
#24 naive Treg %pSTAT5+
d <- rbind(d, data.frame(id=names(pop.stats), cell.type='Naive Treg', do.call('rbind', lapply(pop.stats, function(x) x[24,]))))



load('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/CB00055J_2012-10-03.RData')
print(load('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR/CB00055J_2012-10-03.RData'))



BASE.DIR <- '~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3'
BASE.DIR <- '~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/'

load('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/transform-w1.RData')

par(mfrow=c(2,5))
for (f in list.files(file.path(BASE.DIR,'FCS','pstat5-join'), pattern='.*.RData', full.names=TRUE)) { 
    print(load(f))
    fcs.data <- baseline.relative.pstat5(fcs.data)
    t <- tree(diff.PSTAT5.4 ~ FSCA + SSCA, data=data.frame(fcs.data))
    t <- prune.tree( t, best=3 )
    w <- as.factor(t$where)
    levels(w) <- 1:length(unique(w))
    smoothPlot(fcs.data[,c('FSCA','SSCA')], classification=w, ellipses=FALSE, chull.lwd=4)
}




base.dir <- '~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/' 

CORE.FMARKERS <- c('CD25','CD3','CD4','CD45RA','CD56','CD8','FOXP3')
CORE.MARKERS <- c('SSCA','FSCA','CD25','CD3','CD4','CD45RA','CD56','CD8','FOXP3')
 
load("~/dunwich/Projects/IL2/transforms.RData")
#scale scatter to be on a similar scale as fluorescence
transforms[['SSCA']] <- function(x) 10*(x-min(x))/(max(x)-min(x))
transforms[['FSCA']] <- function(x) 10*(x-min(x))/(max(x)-min(x))

f <- function(individual,date) {
  FILES <- sprintf('%s.RData', paste(individual, DOSES_, date, sep='_'))
  load(file.path(base.dir,'RData',FILES[[1]]))
  fcs.data <- applyTransforms(fcs.data,transforms)
  smoothPlot(fcs.data[,c('CD3','CD56')],outliers=TRUE)
}
 
figure.labels <- iter(paste(letters,')',sep=''))
par(mfrow=c(1,2))
#good
individual <- 'CB00086S'
date <- '2012-09-18' 
f(individual,date)
title(nextElem(figure.labels), adj=0)
#bad
individual <- 'CB00406Q'
date <- '2012-06-12'
f(individual,date)
title(nextElem(figure.labels), adj=0)

pdf('~/Thesis/figures/lymphocytes-all-markers.pdf',width=10,height=10)
plotClusters(fcs.data[as.logical(CLR[,'Lymphocytes']),CORE.FMARKERS],posteriors=CLR[as.logical(CLR[,'Lymphocytes']),CELL.TYPES],outliers=TRUE,ellipses=FALSE,chull.lwd=2,clusters.col=CELL.TYPES.COL)
dev.off()







