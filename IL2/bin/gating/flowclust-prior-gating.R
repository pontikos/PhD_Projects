#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tree))
suppressPackageStartupMessages(library(rpart))
#suppressPackageStartupMessages(library(mvpart))
suppressPackageStartupMessages(library(partykit))
suppressPackageStartupMessages(library(diptest))
suppressPackageStartupMessages(library(flowClust))
source('~nikolas/bin/FCS/fcs.R',chdir=T)

option_list <- list( 
make_option(c("--in.file"), default=NULL, help = ".RData file"),
make_option(c("--plot.file"), default=NULL, help = ".RData file"),
make_option(c("--gate.file"), default=NULL, help = ".RData file"),
make_option(c("--prior.file"), default=NULL, help = ".RData file")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

if (!is.null(opt$channels)) {
    channels <- unlist(strsplit(opt$channels, ","))
} else  {
    channels <- NULL
}

print(file.name <- opt$in.file)
load(file.name)

print(gate.file <- opt$gate.file)

# load priors
#load('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/flowclust-priors.RData')
load(opt$prior.file)

fcs.data <- applyTransforms(fcs.data,transforms)

CLR <- matrix(0,nrow(fcs.data),9)
colnames(CLR) <- c("Lymphocytes", "Single cells", "CD4", "Memory", "Naive", "Naive Eff", "Naive Treg", "Memory Eff", "Memory Treg")

level=.9
u.cutoff=.5
B=500
kappa=1

#
save.gate <- function() {
    #save( CLR, file=file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3','flowclust-prior-manual-gates','CLR',basename(file.name)) )
    save( CLR, file=gate.file )
}

### gating plot
print(plot.file <- opt$plot.file)
pdf(plot.file)
par(mfrow=c(3,2), cex.lab=2, cex.main=2, las=1, mar=c(6,6,2,1))
figure.labels <- iter(paste(letters,')',sep=''))
#lymphocytes
res <- flowClust::flowClust(fcs.data,K=prior[['Lymphocytes']]$K,varNames=c('FSCA','SSCA'),B=B,level=level,u.cutoff=u.cutoff,trans=0,lambda=1,control=list(B.lambda=0),usePrior='yes', prior=prior[['Lymphocytes']])
classification <- MAP(fcs.data, res)
smoothPlot( fcs.data[,c('SSCA','FSCA')], classification=classification, ellipses=TRUE, chulls=FALSE)
title(nextElem(figure.labels), adj=0) 
lymphocytes <- as.numeric(classification==2)
CLR[,'Lymphocytes'] <- lymphocytes
fcs.data <- fcs.data[as.logical(lymphocytes),]
save.gate()
#single cells 
res <- flowClust::flowClust(fcs.data,K=prior[['Single cells']]$K,varNames=c('SSCH','SSCW'),B=B,level=level,u.cutoff=u.cutoff,trans=0,lambda=1,control=list(B.lambda=0),usePrior='yes',prior=prior[['Single cells']])
classification <- MAP(fcs.data, res)
single.cells <- as.numeric(classification==2)
CLR[as.logical(lymphocytes),'Single cells'] <- as.numeric(single.cells)
smoothPlot( fcs.data[,c('SSCH','SSCW')], classification=single.cells, ellipses=TRUE, chulls=FALSE)
title(nextElem(figure.labels), adj=0) 
fcs.data <- fcs.data[as.logical(single.cells),]
save.gate()
# CD4+ lymphocytes
#
res <- flowClust::flowClust(fcs.data,K=prior[['CD4']]$K,varNames=c('CD4','SSCA'),B=B,level=level,u.cutoff=u.cutoff,trans=0,lambda=1,control=list(B.lambda=0),usePrior='yes', prior=prior[['CD4']])
classification <- MAP(fcs.data, res)
cat('CD4\n')
classification <- MAP(fcs.data, res)
cd4 <- as.numeric(classification==2)
CLR[as.logical(CLR[,'Single cells']),'CD4'] <- cd4
smoothPlot( fcs.data[,c('CD4','SSCA')], classification=cd4, ellipses=TRUE, chulls=FALSE)
title(nextElem(figure.labels), adj=0) 
fcs.data <- fcs.data[as.logical(cd4),]
save.gate()
# memory / naive
res <- flowClust::flowClust(fcs.data,K=prior[['Memory / Naive']]$K,varNames=c('CD45RA','SSCA'),B=B,level=level,u.cutoff=u.cutoff,trans=0,lambda=1,control=list(B.lambda=0),usePrior='yes',prior=prior[['Memory / Naive']])
classification <- MAP(fcs.data, res)
smoothPlot( fcs.data[,c('CD45RA','SSCA')], classification=classification, ellipses=TRUE, chulls=FALSE)
cat('CD45RA\n')
#memory <- as.numeric(classification==which.min(tapply(fcs.data[,'CD45RA'], classification, mean)))
memory <- as.numeric(classification==2)
CLR[as.logical(CLR[,'CD4']),'Memory'] <- memory
#naive <-  as.numeric(classification==which.max(tapply(fcs.data[,'CD45RA'], classification, mean)))
naive <-  as.numeric(classification==3)
CLR[as.logical(CLR[,'CD4']),'Naive'] <- naive
title(nextElem(figure.labels), adj=0)
naive <- fcs.data[as.logical(naive),]
memory <- fcs.data[as.logical(memory),]
save.gate()
# naive effector / treg
cat('Naive\n')
res <- flowClust::flowClust(naive,varNames=c('CD25','FOXP3'),B=B,level=level,u.cutoff=u.cutoff,trans=0,lambda=1,control=list(B.lambda=0),usePrior='yes', prior=prior[['Naive Eff / Treg']], K=prior[['Naive Eff / Treg']]$K)
classification <- MAP(naive, res)
print(cd25.mfi <- tapply(naive[,'CD25'], classification, mean))
print(foxp3.mfi <- tapply(naive[,'FOXP3'], classification, mean))
# naive eff
naive.eff <- as.numeric(classification==which.min(cd25.mfi))
CLR[as.logical(CLR[,'Naive']),'Naive Eff'] <- naive.eff
# naive treg
naive.tregs <- as.numeric(classification==which.max(foxp3.mfi))
CLR[as.logical(CLR[,'Naive']),'Naive Treg'] <- naive.tregs
smoothPlot( naive[,c('CD25','FOXP3')], classification=naive.eff+naive.tregs*2, ellipses=TRUE, chulls=TRUE, outliers=TRUE )
title(nextElem(figure.labels), adj=0)
save.gate()
# memory effector / treg
cat('Memory\n')
res <- flowClust::flowClust(memory,varNames=c('CD25','FOXP3'),K=prior[['Memory Eff / Treg']]$K,B=B,level=level,u.cutoff=u.cutoff,trans=0,lambda=1,control=list(B.lambda=0),usePrior='yes',prior=prior[['Memory Eff / Treg']])
classification <- MAP(memory, res)
print(cd25.mfi <- tapply(memory[,'CD25'], classification, mean))
print(foxp3.mfi <- tapply(memory[,'FOXP3'], classification, mean))
# memory eff
memory.eff <- as.numeric(classification==which.min(cd25.mfi))
CLR[as.logical(CLR[,'Memory']),'Memory Eff'] <- memory.eff
# memory treg
memory.tregs <- as.numeric(classification==which.max(foxp3.mfi))
CLR[as.logical(CLR[,'Memory']),'Memory Treg'] <- memory.tregs
smoothPlot( memory[,c('CD25','FOXP3')], classification=memory.eff+memory.tregs*2, ellipses=TRUE, chulls=TRUE, outliers=TRUE )
title(nextElem(figure.labels), adj=0)
dev.off()



