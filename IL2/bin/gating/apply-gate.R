#!/usr/bin/env Rscript

# Given a flow file and a CLR file, apply the gating in the CLR to the flow file.
# The challenge is that the number of rows in the two files are different.
# Since the dimensions of the flow file (file to be gated) and the CLR file (file containing the gating)
# are different, I use the classification.to.chull function which returns the coordinates of the convex hull points
# (i.e the points lying at the periphery of the gate) in the original file (the file on which the gating was performed).
# Note that the 2 channels on which the gate is defined need
# to be specified.
# The output will be a CLR file.

# Gating dimensions
# 1. Lymphocytes, SSCA, FSCA
# 2. Single cells, SSCH, SSCW
# 3. CD4, CD4, SSCA
# 3.1 Memory, CD45RA, SSCA
# 3.2 Naive, CD45RA, SSCA
# 3.1.1 Memory Eff, CD25, FOXP3
# 3.1.2 Memory Treg, CD25, FOXP3
# 3.2.1 Naive Eff, CD25, FOXP3
# 3.2.2 Naive Treg, CD25, FOXP3

#BASEDIR=~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5; for x in `cat $BASEDIR/individual-date.csv`; do individual=`echo $x | cut -d, -f1`; date=`echo $x | cut -d, -f2`; Rscript ~nikolas/Projects/IL2/bin/gating/apply-gate.R --individual $individual --date $date --base.dir $BASEDIR; done

suppressPackageStartupMessages(library(optparse))
source('~nikolas/bin/FCS/fcs.R',chdir=T)

#option_list <- list( 
#make_option(c("--in.gate"), default=NULL, help = ".RData CLR file"),
#make_option(c("--out.gate"), default=NULL, help = ".RData CLR file")
#)

#individual
#date

option_list <- list( 
make_option(c("--individual"), default=NULL, help = "individual"),
make_option(c("--date"), default=NULL, help = "date"),
make_option(c("--base.dir"), default=NULL, help = "base dir")
)
OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

individual <- opt$individual
date <- opt$date

#gate.dir <- '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR'
gate.dir <- file.path(opt$base.dir, 'CLR')

#in.gate
in.gate <- sprintf('%s.RData',paste(individual, date, sep='_'))
in.gate <- file.path(gate.dir,in.gate)

print(basename(in.gate))
# gate file CLR
print(load(in.gate))
gate.CLR <- CLR
# gate file FCS data
# first look for file here
fcs.dir <- file.path(opt$base.dir, 'RData', 'pstat5-join')
load(file.path(fcs.dir,basename(in.gate)))
gate.fcs.data <- fcs.data
# make sure the num of rows of the CLR and FCS file match
print(nrow(CLR)==nrow(gate.fcs.data))

#
f <- function(fcs.data, channels, main, outliers=FALSE, plot.gates=NULL) {
    X <- fcs.data[,channels]
    xquant <- quantile(fcs.data[,channels[[1]]],probs=seq(0,1,.01))
    yquant <- quantile(fcs.data[,channels[[2]]],probs=seq(0,1,.01))
    print(xlim <- c(xquant[['1%']],xquant[['99%']]))
    print(ylim <- c(yquant[['1%']],yquant[['99%']]))
    if (!outliers)
        smoothScatter( fcs.data[,channels], xlab=channels[[1]], ylab=channels[[2]], main=main, colramp=colorRampPalette(c('white','blue','green','yellow','orange','red')), nrpoints=0, cex.lab=2, cex.main=2, ylim=ylim, xlim=xlim )
    else
        smoothScatter( fcs.data[,channels], xlab=channels[[1]], ylab=channels[[2]], main=main, colramp=colorRampPalette(c('white','blue','green','yellow','orange','red')), nrpoints=0, cex.lab=2, cex.main=2 )
    gates <- list()
    for (gate in plot.gates) {
        ch.points <- classification.to.chull(gate.fcs.data[,channels],CLR[,gate])
        polygon(ch.points)
        gates[[gate]] <- as.logical(point.in.polygon(X[,1],X[,2],ch.points[,1], ch.points[,2]))
    }
    return(gates)
}


f <- function(fcs.data, channels, main, outliers=FALSE, plot.gates=NULL) {
    X <- fcs.data[,channels]
    gates <- list()
    for (gate in plot.gates) {
        print(gate)
        print(dim(gate.fcs.data))
        print(dim(gate.CLR))
        ch.points <- classification.to.chull(gate.fcs.data[,channels],gate.CLR[,gate])
        gates[[gate]] <- as.logical(point.in.polygon(X[,1],X[,2],ch.points[,1], ch.points[,2]))
    }
    return(gates)
}

#out.gate
DOSES <- c( '0U', '01U', '10U', '1000U')


for (dose in DOSES) {

out.gate <- sprintf('%s.RData',paste(individual, dose, date, sep='_'))
out.gate <- file.path(gate.dir,out.gate)

print(basename(out.gate))

# file on which gating is to applied
# first look for file here
fcs.dir <- file.path(opt$base.dir, 'RData')
load(file.path(fcs.dir,basename(out.gate))) 
print(dim(x <- fcs.data))

# Lymphocytes
channels <- c('FSCA','SSCA')
gates <- f(x, channels, main='Lymphocytes', plot.gates='Lymphocytes')
x <- x[gates[['Lymphocytes']],]
# Single cells
channels <- c('SSCH','SSCW')
gates <- c(gates, f(x, channels, main='Single cells', outliers=FALSE, plot.gates='Single cells'))
x <- x[gates[['Single cells']],]
# CD4
channels <- c('CD4','SSCA')
gates <- c(gates, f(x, channels, main='CD4+', outliers=FALSE, plot.gates='CD4'))
x <- x[gates[['CD4']],]
# Memory / Naive
channels <- c('CD45RA','SSCA')
gates <- c(gates, f(x, channels, main='Memory / Naive', outliers=FALSE, plot.gates=c('Memory','Naive')))
x.naive <- x[ gates[['Naive']] , ]
x.memory <- x[ gates[['Memory']] , ] 
# Naive Eff / Treg
channels <- c('CD25','FOXP3')
gates <- c(gates, f(x.naive, channels, main='Naive Eff / TReg', outliers=FALSE, plot.gates=c('Naive Eff','Naive Treg')))
# Memory Eff / Treg
channels <- c('CD25','FOXP3')
gates <- c(gates, f(x.memory, channels, main='Memory Eff / TReg', outliers=FALSE, plot.gates=c('Memory Eff','Memory Treg')))

#
CLR <- matrix(0, nrow=nrow(fcs.data), ncol=9)
colnames(CLR) <- c('Lymphocytes', 'Single cells', 'CD4', 'Memory', 'Naive', 'Naive Eff', 'Naive Treg', 'Memory Eff', 'Memory Treg') 
CLR[,'Lymphocytes'] <- as.numeric(gates[['Lymphocytes']])
CLR[which(as.logical(CLR[,'Lymphocytes'])),'Single cells'] <- as.numeric(gates[['Single cells']])
CLR[which(as.logical(CLR[,'Single cells'])),'CD4'] <- as.numeric(gates[['CD4']])
CLR[which(as.logical(CLR[,'CD4'])),'Naive'] <- as.numeric(gates[['Naive']])
CLR[which(as.logical(CLR[,'CD4'])),'Memory'] <- as.numeric(gates[['Memory']])
CLR[which(as.logical(CLR[,'Naive'])),'Naive Eff']  <- as.numeric(gates[['Naive Eff']])
CLR[which(as.logical(CLR[,'Naive'])),'Naive Treg']  <- as.numeric(gates[['Naive Treg']])
CLR[which(as.logical(CLR[,'Memory'])),'Memory Eff']  <- as.numeric(gates[['Memory Eff']])
CLR[which(as.logical(CLR[,'Memory'])),'Memory Treg']  <- as.numeric(gates[['Memory Treg']])
save(CLR, file=out.gate)

}

