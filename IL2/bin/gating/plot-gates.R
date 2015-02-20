#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))

# Plots the gates from the FCS and corresponding CLR file following the
# Lymphocyte -> Single cells -> CD4+ ...
# hierarchy.


#for x in RData/pstat5-join/*.RData; do x=`basename $x`; Rscript ~nikolas/Projects/IL2/bin/gating/plot-gates.R --in.file  RData/pstat5-join/$x --plot.file CLR/Plots/${x%.RData}.pdf --gate.file CLR/$x ; done

option_list <- list( 
make_option(c("--in.file"), default=NULL, help = ".RData file"),
make_option(c("--gate.file"), default=NULL, help = ".RData file"),
make_option(c("--plot.file"), default=NULL, help = ".RData file")
)
OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt


#in.file
print(basename(opt$in.file))
load(opt$in.file)
print(colnames(fcs.data))

#gate.file
print(basename(opt$gate.file))
load(opt$gate.file)
print(colnames(CLR))

if (nrow(fcs.data)!=nrow(CLR)) stop('number of rows of in.file and gate.file do not match!')

#plot.file
print(plot.file <- opt$plot.file)

source('~nikolas/bin/FCS/fcs.R',chdir=T)
print(load('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/transform-w1.RData'))
print(names(transforms))
#fcs.data <- applyTransforms(fcs.data, transforms) 
fcs.data <- cbind(fcs.data[,c('FSCA','SSCA')], applyTransforms(fcs.data[,-grep('FSCA|SSCA',colnames(fcs.data))], transforms))

f <- function(fcs.data, channels, main, outliers=FALSE) {
    xquant <- quantile(fcs.data[,channels[[1]]],probs=seq(0,1,.01))
    yquant <- quantile(fcs.data[,channels[[2]]],probs=seq(0,1,.01))
    print(xlim <- c(xquant[['1%']],xquant[['99%']]))
    print(ylim <- c(yquant[['1%']],yquant[['99%']]))
    smoothPlot( fcs.data[,channels], xlab=channels[[1]], ylab=channels[[2]], main=main, outliers=outliers )
}


plot.gate.chull <- function(X, classification, col) {
    print(table(classification==1))
    X1 <- X[which(classification==1),]
    p <- X1[chull(X1),]
    p <- rbind(p, p)
    lines(p, col=col, lwd=2)
}

pdf(plot.file)
par(mfrow=c(3,2), cex.lab=2, cex.main=2, las=1, mar=c(6,6,2,1))
x <- fcs.data
figure.labels <- iter(paste(letters,')',sep='')) 
# Lymphocytes
channels <- c('FSCA','SSCA')
f(x, channels, main='Lymphocytes')
plot.gate.chull(fcs.data[,channels], classification=CLR[,'Lymphocytes'], col='red')
title(nextElem(figure.labels), adj=0)
x <- fcs.data[ which(as.logical(CLR[,'Lymphocytes'])) , ] 
# Single cells
channels <- c('SSCH','SSCW')
f(x, channels, main='Single cells', outliers=FALSE)
plot.gate.chull(fcs.data[,channels], classification=CLR[,'Single cells'], col='red')
title(nextElem(figure.labels), adj=0)
x <- fcs.data[ which(as.logical(CLR[,'Single cells'])) , ] 
# CD4
channels <- c('CD4','SSCA')
f(x, channels, main='CD4+', outliers=FALSE)
plot.gate.chull(fcs.data[,channels], classification=CLR[,'CD4'], col='red')
title(nextElem(figure.labels), adj=0)
x <- fcs.data[ which(as.logical(CLR[,'CD4'])) , ] 
# Memory / Naive
channels <- c('CD45RA','SSCA')
f(x, channels, main='Memory / Naive', outliers=TRUE)
plot.gate.chull(fcs.data[,channels], classification=CLR[,'Memory'], col='black')
plot.gate.chull(fcs.data[,channels], classification=CLR[,'Naive'], col='red')
title(nextElem(figure.labels), adj=0)
x.naive <- fcs.data[ as.logical(CLR[,'Naive']) , ]
x.memory <- fcs.data[ as.logical(CLR[,'Memory']) , ] 
# Naive Eff / Treg
channels <- c('CD25','FOXP3')
f(x.naive, channels, main='Naive Eff / TReg', outliers=TRUE)
plot.gate.chull(fcs.data[,channels], classification=CLR[,'Naive Eff'], col='green')
plot.gate.chull(fcs.data[,channels], classification=CLR[,'Naive Treg'], col='blue')
title(nextElem(figure.labels), adj=0)
# Memory Eff / Treg
channels <- c('CD25','FOXP3')
f(x.memory, channels, main='Memory Eff / TReg', outliers=TRUE)
plot.gate.chull(fcs.data[,channels], classification=CLR[,'Memory Eff'], col='black')
plot.gate.chull(fcs.data[,channels], classification=CLR[,'Memory Treg'], col='red')
title(nextElem(figure.labels), adj=0)
dev.off()

