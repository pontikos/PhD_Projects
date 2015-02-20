#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tree))
suppressPackageStartupMessages(library(rpart))
#suppressPackageStartupMessages(library(mvpart))
suppressPackageStartupMessages(library(partykit))
suppressPackageStartupMessages(library(diptest))
source('~nikolas/bin/FCS/fcs.R',chdir=T)

option_list <- list( 
make_option(c("--in.file"), default=NULL, help = ".RData file")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

if (!is.null(opt$channels)) {
    channels <- unlist(strsplit(opt$channels, ","))
} else  {
    channels <- NULL
}

file.name <- opt$in.file


# Follows the manual gating using CART.

load(file.name)
load('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/transform-w1.RData')
fcs.data <- applyTransforms(fcs.data,transforms)
fcs.data <- baseline.relative.pstat5(fcs.data) 
clr <- 'magnetic-manual-gates2'
gate.dir <- sprintf('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/%s/CLR/',clr)
load(file.path(gate.dir,basename(file.name)))

CLR <- matrix(0,nrow(fcs.data),9)
colnames(CLR) <- c("Lymphocytes", "Single cells", "CD4", "Memory", "Naive", "Naive Eff", "Naive Treg", "Memory Eff", "Memory Treg")

### gating plot
pdf(file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3','cart-manual-gates','Plots',sprintf('%s.pdf',gsub('.RData','',basename(file.name)))))
par(mfrow=c(3,2), cex.lab=2, cex.main=2, las=1, mar=c(6,6,2,1))
figure.labels <- iter(paste(letters,')',sep=''))
#lymphocytes
t <- tree::tree( diff.PSTAT5.4 ~ FSCA + SSCA, data=data.frame(fcs.data) )
t <- prune.tree( t, best=3 )
w <- as.factor(t$where)
levels(w) <- 1:length(unique(w))
smoothPlot( fcs.data[,c('SSCA','FSCA')], classification=w, ellipses=TRUE, chulls=FALSE)
title(nextElem(figure.labels), adj=0) 
lymphocytes <- points.in.ellipse(fcs.data[,c('SSCA','FSCA')], w==2) 
CLR[,'Lymphocytes'] <- lymphocytes
fcs.data <- fcs.data[as.logical(lymphocytes),]
#single cells
t <- tree::tree( diff.PSTAT5.4 ~ SSCH + SSCW, data=data.frame(fcs.data) )
t <- prune.tree( t, best=2 )
w <- as.factor(t$where)
levels(w) <- 1:length(unique(w))
p1 <- points.in.ellipse(fcs.data[,c('SSCH','SSCW')], w==1)
p2 <- points.in.ellipse(fcs.data[,c('SSCH','SSCW')], w==2)
single.cells <- as.logical(p1) | as.logical(p2) 
CLR[as.logical(lymphocytes),'Single cells'] <- as.numeric(single.cells)
smoothPlot( fcs.data[,c('SSCH','SSCW')], classification=single.cells, ellipses=TRUE, chulls=FALSE)
title(nextElem(figure.labels), adj=0) 
fcs.data <- fcs.data[as.logical(single.cells),]
# CD4+ lymphocytes
t <- tree::tree( diff.PSTAT5.3 ~ CD4 + SSCA, data=data.frame(fcs.data), minsize = 100)
if (!'CD4' %in% t$frame$var)
t <- tree::tree( diff.PSTAT5.4 ~ CD4 + SSCA, data=data.frame(fcs.data), minsize = 100)
print(t)
cat('CD4\n')
cd4.positive <- unique(as.numeric(gsub('>|<','',t$frame[which(t$frame$var=='CD4'),'splits'])))
cd4 <- as.numeric(fcs.data[,'CD4'] > cd4.positive)
CLR[as.logical(CLR[,'Single cells']),'CD4'] <- cd4
smoothPlot( fcs.data[,c('CD4','SSCA')], classification=cd4, ellipses=TRUE, chulls=FALSE)
title(nextElem(figure.labels), adj=0) 
fcs.data <- fcs.data[as.logical(cd4),]
# memory / naive
t <- tree::tree( diff.PSTAT5.4 ~ CD45RA + SSCA, data=data.frame(fcs.data) )
print(t)
cat('CD45RA\n')
print(cd45ra.positive <- (min(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='CD45RA'),'splits'])))))
memory <- as.numeric(fcs.data[,'CD45RA']<cd45ra.positive)
CLR[as.logical(CLR[,'CD4']),'Memory'] <- memory
print(cd45ra.positive <- (max(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='CD45RA'),'splits'])))))
naive <- as.numeric(fcs.data[,'CD45RA']>cd45ra.positive)
CLR[as.logical(CLR[,'CD4']),'Naive'] <- naive
smoothPlot( fcs.data[,c('CD45RA','SSCA')], classification=memory+naive*2, ellipses=TRUE, chulls=FALSE)
title(nextElem(figure.labels), adj=0)
naive <- fcs.data[as.logical(naive),]
memory <- fcs.data[as.logical(memory),]
# naive effector / treg
t <- tree::tree( diff.PSTAT5.2 ~ CD25 + FOXP3, data=data.frame(naive), minsize = 100 )
if (!'CD25' %in% t$frame$var)
t <- tree::tree( diff.PSTAT5.3 ~ CD25 + FOXP3, data=data.frame(naive), minsize = 100 ) 
if (!'CD25' %in% t$frame$var)
t <- tree::tree( diff.PSTAT5.4 ~ CD25 + FOXP3, data=data.frame(naive), minsize = 100 )  
cat('Naive\n')
print(t)
# naive eff
print(foxp3.positive <- min(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='FOXP3'),'splits'])))))
print(cd25.positive <- min(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='CD25'),'splits'])))))
naive.eff <- as.numeric(naive[,'FOXP3'] < foxp3.positive & naive[,'CD25'] < cd25.positive)
CLR[as.logical(CLR[,'Naive']),'Naive Eff'] <- naive.eff
# naive treg
print(foxp3.positive <- max(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='FOXP3'),'splits'])))))
print(cd25.positive <- max(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='CD25'),'splits'])))))
naive.tregs <- as.numeric(naive[,'FOXP3'] > foxp3.positive & naive[,'CD25'] > cd25.positive)
CLR[as.logical(CLR[,'Naive']),'Naive Treg'] <- naive.tregs
smoothPlot( naive[,c('CD25','FOXP3')], classification=naive.eff+naive.tregs*2, ellipses=FALSE, chulls=TRUE, outliers=TRUE )
title(nextElem(figure.labels), adj=0)
# memory effector / treg
t <- tree::tree( diff.PSTAT5.2 ~ CD25 + FOXP3, data=data.frame(memory) )
cat('Memory\n')
print(t)
# memory eff
print(foxp3.positive <- min(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='FOXP3'),'splits'])))))
print(cd25.positive <- min(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='CD25'),'splits'])))))
memory.eff <- as.numeric(memory[,'FOXP3'] < foxp3.positive & memory[,'CD25'] < cd25.positive)
CLR[as.logical(CLR[,'Memory']),'Memory Eff'] <- memory.eff
# memory treg
print(foxp3.positive <- max(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='FOXP3'),'splits'])))))
print(cd25.positive <- max(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='CD25'),'splits'])))))
memory.tregs <- as.numeric(memory[,'FOXP3'] > foxp3.positive & memory[,'CD25'] > cd25.positive)
CLR[as.logical(CLR[,'Memory']),'Memory Treg'] <- memory.tregs
smoothPlot( memory[,c('CD25','FOXP3')], classification=memory.eff+memory.tregs*2, ellipses=FALSE, chulls=TRUE, outliers=TRUE )
title(nextElem(figure.labels), adj=0)
dev.off()

#CLR file
save(CLR, file=file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cart-manual-gates/CLR/',basename(file.name)))
#save(CLR, file=file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR/',basename(file.name)))

#
load(file.name)
load('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/transform-w1.RData')
fcs.data <- applyTransforms(fcs.data,transforms)
fcs.data <- baseline.relative.pstat5(fcs.data) 
clr <- 'magnetic-manual-gates2'
gate.dir <- sprintf('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/%s/CLR/',clr)
load(file.path(gate.dir,basename(file.name)))

CLR <- matrix(0,nrow(fcs.data),9)
colnames(CLR) <- c("Lymphocytes", "Single cells", "CD4", "Memory", "Naive", "Naive Eff", "Naive Treg", "Memory Eff", "Memory Treg")

### tree plot
pdf(file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3','cart-manual-gates','Plots',sprintf('tree-%s.pdf',gsub('.RData','',basename(file.name)))))
par(mfrow=c(3,2), cex.lab=2, cex.main=2, las=1, mar=c(6,6,2,1))
figure.labels <- iter(paste(letters,')',sep=''))
#lymphocytes
t <- tree::tree( diff.PSTAT5.4 ~ FSCA + SSCA, data=data.frame(fcs.data) )
t <- prune.tree( t, best=3 )
w <- as.factor(t$where)
levels(w) <- 1:length(unique(w))
plot(t)
text(t)
title(nextElem(figure.labels), adj=0) 
lymphocytes <- points.in.ellipse(fcs.data[,c('SSCA','FSCA')], w==2) 
CLR[,'Lymphocytes'] <- lymphocytes
fcs.data <- fcs.data[as.logical(lymphocytes),]
#single cells
t <- tree::tree( diff.PSTAT5.4 ~ SSCH + SSCW, data=data.frame(fcs.data) )
t <- prune.tree( t, best=2 )
w <- as.factor(t$where)
levels(w) <- 1:length(unique(w))
p1 <- points.in.ellipse(fcs.data[,c('SSCH','SSCW')], w==1)
p2 <- points.in.ellipse(fcs.data[,c('SSCH','SSCW')], w==2)
single.cells <- as.logical(p1) | as.logical(p2) 
CLR[as.logical(lymphocytes),'Single cells'] <- as.numeric(single.cells)
plot(t)
text(t)
title(nextElem(figure.labels), adj=0) 
fcs.data <- fcs.data[as.logical(single.cells),]
# CD4+ lymphocytes
t <- tree::tree( diff.PSTAT5.3 ~ CD4 + SSCA, data=data.frame(fcs.data) )
if (!'CD4' %in% t$frame$var)
t <- tree::tree( diff.PSTAT5.4 ~ CD4 + SSCA, data=data.frame(fcs.data) )
print(t)
cat('CD4\n')
cd4.positive <- unique(as.numeric(gsub('>|<','',t$frame[which(t$frame$var=='CD4'),'splits'])))
cd4 <- as.numeric(fcs.data[,'CD4'] > cd4.positive)
CLR[as.logical(CLR[,'Single cells']),'CD4'] <- cd4
plot(t)
text(t)
title(nextElem(figure.labels), adj=0) 
fcs.data <- fcs.data[as.logical(cd4),]
# memory / naive
t <- tree::tree( diff.PSTAT5.4 ~ CD45RA + SSCA, data=data.frame(fcs.data) )
print(t)
cat('CD45RA\n')
print(cd45ra.positive <- (min(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='CD45RA'),'splits'])))))
memory <- as.numeric(fcs.data[,'CD45RA']<cd45ra.positive)
CLR[as.logical(CLR[,'CD4']),'Memory'] <- memory
print(cd45ra.positive <- (max(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='CD45RA'),'splits'])))))
naive <- as.numeric(fcs.data[,'CD45RA']>cd45ra.positive)
CLR[as.logical(CLR[,'CD4']),'Naive'] <- naive
plot(t)
text(t)
title(nextElem(figure.labels), adj=0)
naive <- fcs.data[as.logical(naive),]
memory <- fcs.data[as.logical(memory),]
# naive effector / treg
t <- tree::tree( diff.PSTAT5.2 ~ CD25 + FOXP3, data=data.frame(naive) )
if (!'CD25' %in% t$frame$var)
t <- tree::tree( diff.PSTAT5.3 ~ CD25 + FOXP3, data=data.frame(naive) )
if (!'CD25' %in% t$frame$var)
t <- tree::tree( diff.PSTAT5.4 ~ CD25 + FOXP3, data=data.frame(naive) )
cat('Naive\n')
print(t)
# naive eff
print(foxp3.positive <- min(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='FOXP3'),'splits'])))))
print(cd25.positive <- min(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='CD25'),'splits'])))))
naive.eff <- as.numeric(naive[,'FOXP3'] < foxp3.positive & naive[,'CD25'] < cd25.positive)
CLR[as.logical(CLR[,'Naive']),'Naive Eff'] <- naive.eff
# naive treg
print(foxp3.positive <- max(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='FOXP3'),'splits'])))))
print(cd25.positive <- max(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='CD25'),'splits'])))))
naive.tregs <- as.numeric(naive[,'FOXP3'] > foxp3.positive & naive[,'CD25'] > cd25.positive)
CLR[as.logical(CLR[,'Naive']),'Naive Treg'] <- naive.tregs
plot(t)
text(t)
title(nextElem(figure.labels), adj=0)
# memory effector / treg
t <- tree::tree( diff.PSTAT5.2 ~ CD25 + FOXP3, data=data.frame(memory) )
cat('Memory\n')
print(t)
# memory eff
print(foxp3.positive <- min(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='FOXP3'),'splits'])))))
print(cd25.positive <- min(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='CD25'),'splits'])))))
memory.eff <- as.numeric(memory[,'FOXP3'] < foxp3.positive & memory[,'CD25'] < cd25.positive)
CLR[as.logical(CLR[,'Memory']),'Memory Eff'] <- memory.eff
# memory treg
print(foxp3.positive <- max(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='FOXP3'),'splits'])))))
print(cd25.positive <- max(unique(as.numeric(gsub('>|<','', t$frame[which(t$frame$var=='CD25'),'splits'])))))
memory.tregs <- as.numeric(memory[,'FOXP3'] > foxp3.positive & memory[,'CD25'] > cd25.positive)
CLR[as.logical(CLR[,'Memory']),'Memory Treg'] <- memory.tregs
plot(t)
text(t)
title(nextElem(figure.labels), adj=0)
dev.off()





