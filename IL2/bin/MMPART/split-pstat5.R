#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(mixtools))
suppressPackageStartupMessages(library(RANN))
source('~nikolas/bin/FCS/fcs.R',chdir=T)
source('~nikolas/Projects/IL2/bin/MMPART/mmpart-functions.R',chdir=T)

# only works on joined files
# divide cells as low responders and high responders on pSTAT5 response (i.e baseline subtracted) at 1000U 
# within responders further divide on low/high on pSTAT5 response at 10U
# repeat on pSTAT5 response at 0.1U
# cluster ones which are consistently high, these are the most sensitive cell populations 

#for x in RData/pstat5-join/*.RData; do x=`basename $x`; Rscript ~nikolas/Projects/IL2/bin/split-pstat5.R --in.file  RData/pstat5-join/$x --plot.file first-response-cells/Plots/${x%.RData}.pdf --clr.file CLR/$x  --out.file first-response-cells/$x ; done

option_list <- list( 
make_option(c("--in.file"), default=NULL, help = ".RData file"),
make_option(c("--out.file"), default=NULL, help = ".RData file"),
make_option(c("--plot.file"), default=NULL, help = ".RData file"),
#make_option(c("--cell.subset"), default=NULL, help = "any of the subsets defined in the CLR file"),
make_option(c("--clr.file"),  help = "exclude Lymphocytes identified in the CLR file")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

#file.name <- 'CB01477E_2012-09-26.RData'
file.name <- opt$in.file
print(basename(file.name))
print(clr <- opt$clr.file)

print(plot.file <- opt$plot.file)
print(out.file <- opt$out.file)

cell.subset <- 'Lymphocytes'
#cell.subset <- opt$cell.subset

#file.name <- '~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/RData/pstat5-join/CB00406Q_2012-06-12.RData'
#clr <- '~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/CLR/CB00406Q_2012-06-12.RData'

print(load(file.name))
print(load(clr)) 
#remove cell subset
#fcs.data <- fcs.data[!as.logical(as.logical(CLR[,cell.subset])),] 
load('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/transform-w1.RData')
fcs.data <- applyTransforms(fcs.data, transforms)
fcs.data <- baseline.relative.pstat5(fcs.data) 
#if (!is.null(cell.subset)) fcs.data <- fcs.data[which(as.logical(CLR[,cell.subset])),] 

#X <- split.response(fcs.data[which(!as.logical(CLR[,'CD4'])),], 4) 
#X <- split.response(fcs.data[which(as.logical(CLR[,'CD4'])),], 4)
X <- split.response(fcs.data, 4)

hi.x <- X$d2$d2$d2$d
nn <- nn2(fcs.data,query=hi.x,k=1)

# these are ungated
hi.x <- hi.x[which(CLR[nn$nn.idx,cell.subset]==0),]

#smoothPlot(fcs.data[,c('FSCA','SSCA')])
#points(hi.x[,c('FSCA','SSCA')], col='purple', pch=20)

#hc <- hclust(dist(hi.x[,c('FSCA','SSCA','CD4','CD45RA','FOXP3','CD25')]))
#pairs(hi.x[,c('FSCA','SSCA','CD4','CD45RA','FOXP3','CD25')], col=cutree(hc,k=8), pch=20)
#points(hi.x[,c('SSCH','SSCW')], col=cutree(hc,k=6))
#plotClusters(fcs.data[,c('FSCA','SSCA','CD4','CD45RA','FOXP3','CD25')], outliers=TRUE, plot.points=hi.x, classification=cutree(hc,k=8), upper=contourPlot, ellipses=FALSE, chulls=FALSE)

#save( hi.x, file=file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/first-response-cells2/',basename(file.name)) )
save( hi.x, file=out.file )

#pdf(file.path('~/Plots-tree-pstat5/',gsub('.RData','.pdf',basename(file.name))))
#plotClusters(fcs.data[,c('diff.PSTAT5.4', as.character(unique(f$var[which(f$splits[,'cutleft']!='')])))], classification=t$where, ellipses=FALSE, chull.lwd=2)
#dev.off()

#plotClusters(d[,c('diff.PSTAT5.4','CD45RA','CD25')],classification=t$where, ellipses=FALSE, chull.lwd=2) 
#plotClusters(d[,c('diff.PSTAT5.4','CD45RA','FOXP3','CD25')],classification=prune.tree(t,best=3)$where, ellipses=FALSE, chull.lwd=2)


## the recursive partitioning of pSTAT5 starts at the highest dose then moves
## down to lower doses to find the subset which are the first responders
zoom <- 2
#pdf('~nikolas/dunwich/Projects/IL2/Paper/IL2stim/article/figures/pstat5-rpart.pdf',width=5*zoom,height=3*zoom)
pdf(plot.file,width=5*zoom,height=3*zoom)
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


