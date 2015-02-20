#!/usr/bin/env Rscript
source('~nikolas/bin/FCS/fcs.R')
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("flowCore"))

option_list <- list( 
make_option(c("-f","--in.file"), help = ""),
make_option(c("--out.dir"), default='~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All-pstat5-normalised', help = "")
#make_option(c("-g","--gate"), default=NULL, help = "")

)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

out.dir <- opt$out.dir
#print( out.dir <- sprintf('%s-%s', opt$out.dir, ifelse(!is.null(opt$gate), opt$gate, 'ungated')) ) 

dir.create(out.dir, showWarnings=FALSE)

sliding.window.peaks <- function(d, span=40) {
  y <- d$y
  x <- d$x
  ## sliding a window of size span and returns locations where the middle point in the window is maximum.  
  ## returns the indexes of the peaks
  ind <- c()
  for( i in 1:(length(y)-span)) {
    mid <- i+span%/%2
    if ( y[mid]==max(y[i:(i+span)]) & y[mid]!=y[i] & y[mid]!=y[i+span] ) ind <- c(ind, mid)
  }
  return(x[ind])
}


# returns the top K sliding window peaks
top.sliding.window.peaks <- function(d, K, span=40) {
  y <- d$y
  x <- d$x
  ## sliding a window of size span and returns locations where the middle point in the window is maximum.  
  ## returns the indexes of the peaks
  ind <- c()
  for( i in 1:(length(y)-span)) {
    mid <- i+span%/%2
    if ( y[mid]==max(y[i:(i+span)]) & y[mid]!=y[i] & y[mid]!=y[i+span] ) ind <- c(ind, mid)
  }
  peaks <- cbind(x=x[ind],y=y[ind])
  top.peaks <- peaks[order(peaks[,'y'],decreasing=TRUE)[1:K],]
  top.peaks <- top.peaks[order(top.peaks[,'x']),]
  return(top.peaks)
}


f <- opt$in.file

# ~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/

print(f)
load(f)

chan <- paste('PSTAT5',1:4,sep='.')

# define w in whole sample
#find the peaks in the sample stimulated at 10U
#we expect 2 positive peaks 
#x <- fcs.data[,'PSTAT5.4']
#w <- seq(.1,.9,.1)
#score <- sapply(w, function(w) {
    #(p <- top.sliding.window.peaks(density(logicleTransform(w=w)(x)),2))
    #if ( is.na(p[2,'x']) || is.na(p[1,'x']) ) return(0)
    ##neither of the peaks should be in the negatives
    #if ( (p[2,'x'] < 0) || (p[1,'x'] < 0) ) return(0)
    #return( ( p[2,'x'] - p[1,'x'] ) / abs(p[2,'y'] - p[1,'y']) )
#})
#names(score) <- w
#print(w.best <- as.numeric(names(which.max(score))))
w.best <- .6
w.best <- 1
lgcl <- logicleTransform(w=w.best)
invlgcl <- inverseLogicleTransform(trans=lgcl)


pstat5.normalise <- function(gated.fcs.data,trans.chan='PSTAT5.3',from=-1,to=2,n=512) {
  x <- gated.fcs.data[,trans.chan]
  # x location of the peaks in the gated data
  # depending on the subset, pSTAT5 may not be
  # bimodal
  peaks.x <- (top.sliding.window.peaks(normalised.density(x,from=from,to=to,n=n),2)[,'x'])
  pstat5.norm.transform <- function(y) (cbind(1,y)%*%coefficients(lm(c(0,1) ~ (peaks.x))))
  pstat5.norm <- apply(fcs.data[,chan],2, pstat5.norm.transform)
  return( list(pstat5.norm=pstat5.norm, x=x, peaks.x=(peaks.x), trans.chan=trans.chan, pstat5.norm.transform=pstat5.norm.transform) )
}

fcs.data <- apply(fcs.data, 2, lgcl)
fcs.data <- baseline.relative.pstat5(fcs.data,REPLACE=TRUE)


pstat5.normalisation <- list()
#pstat5.normalisation[['ungated']] <- pstat5.normalise(fcs.data,trans.chan='PSTAT5.4')
load(file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR',basename(f)))
pstat5.normalisation[['Lymphocytes']] <- pstat5.normalise(fcs.data[which(as.logical(CLR[,'Single cells'])),],trans.chan='PSTAT5.4')
#pstat5.normalisation[['Lymphocytes.10U']] <- pstat5.normalise(fcs.data[which(as.logical(CLR[,'Single cells'])),],trans.chan='PSTAT5.3')
#pstat5.normalisation[['Lymphocytes.1000U']] <- pstat5.normalise(fcs.data[which(as.logical(CLR[,'Single cells'])),],trans.chan='PSTAT5.4')
pstat5.normalisation[['CD4']] <- pstat5.normalise(fcs.data[which(as.logical(CLR[,'CD4'])),],trans.chan='PSTAT5.3')

pstat5 <- fcs.data[,paste('PSTAT5',1:4,sep='.')]

print(out.file <- file.path(out.dir,basename(f)))
save(pstat5, lgcl, invlgcl, pstat5.normalisation, file=out.file)


