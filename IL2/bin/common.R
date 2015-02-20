library(tree)
library(iterators)
library(mclust)
library(mixtools)
library(RANN)
source('~nikolas/bin/FCS/fcs.R')
#source('~nikolas/bin/FCS/transforms.R')
#source('~nikolas/bin/FCS/density-estimation.R')
#source('~nikolas/Projects/IL2/bin/gating/manual-gates.R')
library(spade)
library(iterators)
library(flowCore)


REPEATS <- structure(list(individual = c("CB00165D", "CB00366X", "CB00396E", 
"CB00406Q", "CB01484M", "CB01494Y", "CB01495Z", "CB01498C", "CB01503H", 
"CB01504J"), pch = c("a", "b", "c", "d", "e", "f", "g", "h", 
"i", "j"), day1 = c("2012-11-29", "2012-11-07", "2012-09-25", 
"2012-10-16", "2012-09-25", "2012-10-09", "2012-10-09", "2012-10-16", 
"2012-11-07", "2012-11-07"), day2 = c("2013-03-07", "2013-03-27", 
"2013-03-11", "2013-01-22", "2013-03-11", "2013-01-29", "2013-01-29", 
"2013-01-22", "2013-03-07", "2013-03-27"), col = c("#0066FFFF", 
"#FF9900FF", "#00FFFFFF", "#FF0099FF", "#33FF00FF", "#CCFF00FF", 
"#CC00FFFF", "#3300FFFF", "#00FF66FF", "#FF0000FF"), day.diff = structure(c(98, 
140, 167, 98, 167, 112, 112, 98, 120, 140), units = "days", class = "difftime"), 
t1d = c('control', 'case', 'control', 'control', 'case', 'case', 'control', 'case', 'case', 'control')),
.Names = c("individual", 
"pch", "day1", "day2", "col", "day.diff", "t1d"), row.names = c("CB00165D", 
"CB00366X", "CB00396E", "CB00406Q", "CB01484M", "CB01494Y", "CB01495Z", 
"CB01498C", "CB01503H", "CB01504J"), class = "data.frame")

#witthout the decimal dot in 01U
DOSES_ <- c( '0U', '01U', '10U', '1000U')
# with the decimal dot in 0.1U
DOSES <- c( '0U', '0.1U', '10U', '1000U')

# Eff should be TEffs; Treg should be TRegs
CLR.CELL.TYPES <- c("Lymphocytes", "Single cells", "CD4", "Memory", "Memory Eff", "Memory Treg", "Naive", "Naive Eff", "Naive Treg")
CELL.TYPES <- c('Memory Eff', 'Memory Treg', 'Naive Eff', 'Naive Treg')
CELL.TYPES.COL <- c('black','red','darkgreen','darkblue') 
NEW.CELL.TYPES.COL <- c('purple','pink','lightblue','gray','orange')
COLS <- c(CELL.TYPES.COL,NEW.CELL.TYPES.COL)

#BASE.DIR <- '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/Lymphocytes/RData'

blues4 <- blues9[5:9]

# this is the individual I use as an example in the IL2 Chapter
individual <- 'CB00086S'
date <- '2012-09-18' 
base.dir <- '~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/' 
 
 
CORE.FMARKERS <- c('CD25','CD3','CD4','CD45RA','CD56','CD8','FOXP3')
CORE.MARKERS <- c('SSCA','FSCA','CD25','CD3','CD4','CD45RA','CD56','CD8','FOXP3')

FILES <- sprintf('%s.RData', paste(individual, DOSES_, date, sep='_'))

load("~/dunwich/Projects/IL2/transforms.RData")
#scale scatter to be on a similar scale as fluorescence
transforms[['SSCA']] <- function(x) 10*(x-min(x))/(max(x)-min(x))
transforms[['FSCA']] <- function(x) 10*(x-min(x))/(max(x)-min(x))



#################### univariate clusters
plotLymphocytesClusters <- function(file.name=NULL, CLR=CLR, cells=CELL.TYPES, fcs.data=fcs.data, cols=COLS[1:length(cells)]) {
    if(!is.null(file.name)) pdf(file.name)
    r <- c(1,1,2,2,3,3)
    layout(rbind(r,3+r,6+r,rep(10,length(r))),heights=c(4,4,4,1)) 
    par(mai=rep(0.5, 4)) 
    for (marker in CORE.MARKERS)
    smoothPlot1D(fcs.data[,marker], posteriors=CLR[as.logical(CLR[,'Single cells']),cells],outliers=TRUE, clusters.col=cols, main=marker, col='white')
    par(mai=c(0,0,0,0))
    plot.new()
    legend('center',ncol=3,legend=cells,fill=cols,cex=1.5,bty='n')
    if (!is.null(file.name)) dev.off()
}


#################### dose response
plotLymphocytesDoseResponse <- function(file.name=NULL, CLR=CLR, cells=CELL.TYPES, fcs.data=fcs.data, cols=COLS[1:length(cells)]) {
    #cells <- c(CELL.TYPES,cells)
    if (!is.null(file.name)) pdf(file.name)
    ylim <- range( sapply(cells, function(cell.type) colMedians(fcs.data[as.logical(CLR[,cell.type]),paste('diff.PSTAT5',1:4,sep='.')])) )
    plot(NULL, xlim=c(0,3), ylim=ylim, xaxt='n', xlab='dose', ylab='pSTAT5 MFI', main='dose response')
    axis(1, at=0:3, labels=DOSES)
    i <- 1
    for (cell.type in cells) {
        mfi <- fcs.data[which(as.logical(CLR[,cell.type])),paste('diff.PSTAT5',1:4,sep='.')]
        lines(0:3, colMedians( mfi ), col=cols[[i]], lwd=3)
        i <- i+1
    }
    legend('topleft', legend=cells, fill=cols, cex=1.5, bty='n')
    if (!is.null(file.name)) dev.off()
}



#################### univariate clusters
plotNonLymphocytesClusters <- function(file.name=NULL, CLR=CLR, cells=c('Lymphocytes'), fcs.data=fcs.data, cols=c('black',NEW.CELL.TYPES.COL)[1:length(cells)]) {
    #cells <- c('Lymphocytes',cells)
    if(!is.null(file.name)) pdf(file.name)
    r <- c(1,1,2,2,3,3)
    layout(rbind(r,3+r,6+r,rep(10,length(r))),heights=c(4,4,4,1)) 
    par(mai=rep(0.5, 4)) 
    for (marker in CORE.MARKERS)
    smoothPlot1D(fcs.data[,marker], posteriors=CLR[,cells],outliers=TRUE, clusters.col=cols, main=marker, col='white')
    par(mai=c(0,0,0,0))
    plot.new()
    legend('center',ncol=3,legend=cells,fill=cols,cex=1.5,bty='n')
    if (!is.null(file.name)) dev.off()
}


#################### dose response
plotNonLymphocytesDoseResponse <- function(file.name=NULL, CLR=CLR, cells=c('Lymphocytes'), fcs.data=fcs.data, cols=c('black',NEW.CELL.TYPES.COL)[1:length(cells)] ) {
#cells <- c('Lymphocytes',cells)
    if (!is.null(file.name)) pdf(file.name)
    ylim <- range( sapply(cells, function(cell.type) colMedians(fcs.data[as.logical(CLR[,cell.type]),paste('diff.PSTAT5',1:4,sep='.')])) )
    plot(NULL, xlim=c(0,3), ylim=ylim, xaxt='n', xlab='dose', ylab='pSTAT5 MFI', main='dose response')
    axis(1, at=0:3, labels=DOSES)
    i <- 1
    for (cell.type in cells) {
        mfi <- fcs.data[which(as.logical(CLR[,cell.type])),paste('diff.PSTAT5',1:4,sep='.')]
        lines(0:3, colMedians( mfi ), col=cols[[i]], lwd=3)
        i <- i+1
    }
    legend('topleft', legend=cells, fill=cols, cex=1.5, bty='n')
    if (!is.null(file.name)) dev.off()
}


