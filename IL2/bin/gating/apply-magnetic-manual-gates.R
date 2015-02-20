#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
source('~nikolas/bin/FCS/fcs.R',chdir=T)
source('~nikolas/Projects/IL2/bin/gating/manual-gates.R')

option_list <- list( 
make_option(c("--in.file"), default=NULL, help = ".RData file")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

file.name <- opt$in.file
print(basename(file.name))
base.dir <- '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates/'


#file.name <- '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/CB00406Q_2013-01-22.RData' 
print(load(file.name))
fcs.data <- applyTransforms(fcs.data, transforms) 

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


plot.gates <- function(gates, col='black') {
    for (e in gates) lines(ellipse::ellipse(e$Sigma,centre=e$Mu),col=col,lwd=2,lty=2)
    points(do.call('rbind', lapply(gates,function(x)x$Mu)), pch='X', cex=2, col=col)
}


pdf(file.path(base.dir,'Plots',sprintf('%s.pdf',gsub('.RData','',basename(file.name)))))
par(mfrow=c(3,2), cex.lab=2, cex.main=2, las=1, mar=c(6,6,2,1))
x <- fcs.data
figure.labels <- iter(paste(letters,')',sep='')) 
# Lymphocytes
channels <- c('FSCA','SSCA')
e.lymphocytes <- mahalanobis.magnetic.ellipse(x,e.lymphocytes)
f(x, channels, main='Lymphocytes', plot.gates=list(e.lymphocytes))
title(nextElem(figure.labels), adj=0)
lymphocytes.filter <- filter.mahalanobis.ellipse(x,e.lymphocytes)
x <- x[ which(lymphocytes.filter) , ] 
# Single cells
channels <- c('SSCH','SSCW')
e.singlecells <- mahalanobis.magnetic.ellipse(x,e.singlecells)
f(x, channels, main='Single cells', outliers=FALSE, plot.gates=list(e.singlecells))
title(nextElem(figure.labels), adj=0)
singlecells.filter <- filter.mahalanobis.ellipse(x,e.singlecells)
x <- x[ which(singlecells.filter) , ] 
# CD4
channels <- c('CD4','SSCA')
e.cd4 <- mahalanobis.magnetic.ellipse(x,e.cd4)
f(x, channels, main='CD4+', outliers=FALSE, plot.gates=list(e.cd4))
title(nextElem(figure.labels), adj=0)
cd4.filter <- filter.mahalanobis.ellipse(x,e.cd4)
x <- x[ which(cd4.filter) , ] 
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
# Naive Eff / Treg
channels <- c('CD25','FOXP3')
new.ellipses <- em.magnetic.ellipses(x.naive, list(e.naive.conv, e.naive.tregs), iterations=100, background=.001)
e.naive.conv <- new.ellipses[[1]]
e.naive.tregs <- new.ellipses[[2]]
f(x.naive, channels, main='Naive Eff / TReg', outliers=TRUE)
lines(ellipse::ellipse(e.naive.conv$Sigma,centre=e.naive.conv$Mu),col='darkgreen',lwd=3)
lines(ellipse::ellipse(e.naive.tregs$Sigma,centre=e.naive.tregs$Mu),col='blue',lwd=3)
title(nextElem(figure.labels), adj=0)
naive.eff.filter <- filter.mahalanobis.ellipse(x.naive, e.naive.conv)
naive.tregs.filter <- filter.mahalanobis.ellipse(x.naive, e.naive.tregs) 
# Memory Eff / Treg
channels <- c('CD25','FOXP3')
new.ellipses <- em.magnetic.ellipses(x.memory, list(e.memory.conv, e.memory.tregs), iterations=10, background=.0005)
e.memory.conv <- new.ellipses[[1]]
e.memory.tregs <- new.ellipses[[2]]
f(x.memory, channels, main='Memory Eff / TReg', outliers=TRUE)
lines(ellipse::ellipse(e.memory.conv$Sigma,centre=e.memory.conv$Mu),col='black',lwd=3)
lines(ellipse::ellipse(e.memory.tregs$Sigma,centre=e.memory.tregs$Mu),col='red',lwd=3)
title(nextElem(figure.labels), adj=0)
memory.eff.filter <- filter.mahalanobis.ellipse(x.memory, e.memory.conv)
memory.tregs.filter <- filter.mahalanobis.ellipse(x.memory, e.memory.tregs)
dev.off()

#memory
#smoothPlot(fcs.data[which(as.logical(CLR[,'Memory'])),c('CD25','FOXP3')], outliers=TRUE)
#points( fcs.data[which(as.logical(CLR[,'Memory Treg'])), c('CD25','FOXP3')], pch=20)
#points( fcs.data[which(as.logical(CLR[,'Memory Eff'])), c('CD25','FOXP3')], pch=20, col='red')
#naive
#smoothPlot(fcs.data[which(as.logical(CLR[,'Naive'])),c('CD25','FOXP3')], outliers=TRUE)
#points( fcs.data[which(as.logical(CLR[,'Naive Treg'])), c('CD25','FOXP3')], pch=20)
#points( fcs.data[which(as.logical(CLR[,'Naive Eff'])), c('CD25','FOXP3')], pch=20, col='red')


CLR <- matrix(0, nrow=nrow(fcs.data), ncol=9)
colnames(CLR) <- c('Lymphocytes', 'Single cells', 'CD4', 'Memory', 'Naive', 'Naive Eff', 'Naive Treg', 'Memory Eff', 'Memory Treg') 
CLR[,'Lymphocytes'] <- as.numeric(lymphocytes.filter)
CLR[which(as.logical(CLR[,'Lymphocytes'])),'Single cells'] <- as.numeric(singlecells.filter)
CLR[which(as.logical(CLR[,'Single cells'])),'CD4'] <- as.numeric(cd4.filter) 
CLR[which(as.logical(CLR[,'CD4'])),'Naive'] <- as.numeric(naive.filter)
CLR[which(as.logical(CLR[,'CD4'])),'Memory'] <- as.numeric(memory.filter)
CLR[which(as.logical(CLR[,'Naive'])),'Naive Eff']  <- as.numeric(naive.eff.filter)
CLR[which(as.logical(CLR[,'Naive'])),'Naive Treg']  <- as.numeric(naive.tregs.filter)
CLR[which(as.logical(CLR[,'Memory'])),'Memory Eff']  <- as.numeric(memory.eff.filter)
CLR[which(as.logical(CLR[,'Memory'])),'Memory Treg']  <- as.numeric(memory.tregs.filter)

save(CLR, file=file.path(base.dir,'CLR',basename(file.name)))



