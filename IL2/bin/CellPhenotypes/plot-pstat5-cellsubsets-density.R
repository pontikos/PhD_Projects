#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))

option_list <- list( 
make_option(c("--individual"), default=NULL, help = "individual"),
make_option(c("--date"), default=NULL, help = "date")
)
OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

individual <- opt$individual
day <- opt$date

source('~nikolas/Projects/IL2/bin/common.R')
DOSES <- gsub('\\.', '', DOSES)

pstat5dens <- data.frame()
#
for (gate in CELL.TYPES) {
   print(gate)
   for (dose in DOSES) {
      print(dose)
      print((load(sprintf('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/RData/%s_%s_%s.RData',individual,dose,day))))
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

base.dir <- '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/'
#
print(plot.file <- file.path(base.dir,'Plots',sprintf('pstat5_%s_%s.pdf', individual, day)))
pdf(plot.file)
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
      print(dose <- DOSES[i])
      lines(x, X[which(X$day==day & X$dose==dose),grep('^X',colnames(X))],lwd=i,col=blues4[i])
    }
  }
  cat('\n')
  legend('topright',DOSES, lwd=1:4, col=blues4[1:4])
}
dev.off()


#source('~nikolas/bin/FCS/fcs.R') 
#individual <- 'CB01422V'
#day <- '2012-10-23' 
#blues4 <- blues9[5:9] 
## not joined
#par(mfrow=c(2,2))
##
#for (gate in CELL.TYPES) {
    #print(gate)
  #i <- 1
  #for (dose in DOSES) {
      #print(dose)
      #print((load(sprintf('/chiswick/data/store/facs/Tony-RData/PSTAT5-CD25-CD45RA-CD4-FOXP3/%s_%s_%s.RData',individual,dose,day))))
      #print(dim(fcs.data))
      #print(load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3','magnetic-manual-gates2/CLR', sprintf('%s.RData',paste(individual, dose, day, sep='_')))))
      #print(length(g <- which(as.logical(CLR[,gate]))))
      #X <- logicleTransform(w=1)(fcs.data[g,'PSTAT5'])
      #if (dose=='0U') {
          #plot(normalised.density(X), main=gate, xlim=c(0,3), lwd=1, col=blues4[1])
          #abline(v=quantile(X,.99))
      #} else {
          #lines(normalised.density(X), lwd=i, col=blues4[i])
      #}
      #i <- i+1
  #}
#} 
## joined
#par(mfrow=c(2,2))
##
#for (gate in CELL.TYPES) {
    #print(gate)
  #i <- 1
  #for (dose in DOSES) {
      #print(dose)
      #print((load(sprintf('/chiswick/data/store/facs/Tony-RData/PSTAT5-CD25-CD45RA-CD4-FOXP3/%s_%s_%s.RData',individual,dose,day))))
      #print(dim(fcs.data))
      #print(load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3','magnetic-manual-gates2/CLR', sprintf('%s.RData',paste(individual, dose, day, sep='_')))))
      #print(length(g <- which(as.logical(CLR[,gate]))))
      #X <- logicleTransform(w=1)(fcs.data[g,'PSTAT5'])
      #if (dose=='0U') {
          #plot(normalised.density(X), main=gate, xlim=c(0,3), lwd=1, col=blues4[1])
          #abline(v=quantile(X,.99))
      #} else {
          #lines(normalised.density(X), lwd=i, col=blues4[i])
      #}
      #i <- i+1
  #}
#}
