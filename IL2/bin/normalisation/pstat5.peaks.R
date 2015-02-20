#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

option_list <- list( 
make_option(c("-f","--fcs"), help = "fcsfile to parse"),
make_option(c("-p","--parameter"), default='pSTAT5', help=" [default %default]"),
make_option(c("-k", "--peaks"), default=2, help=" [default %default]"),
make_option(c("--plot"), default=NULL, help=" [default %default]")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

if (is.null(opt$fcs)) stop("Not FCS file specified on command line!")
if (is.null(opt$parameter)) stop("Specify at least on parameter!")

p <- opt$parameter
k <- tolower(opt$peaks)

#suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowCore, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(source('~nikolas/bin/FCS/fcs.R'))

fcs.data <- suppressWarnings(read.FCS(opt$fcs, channels=p))
x <- sort(getChannels(fcs.data, p))

suppressPackageStartupMessages(library("cluster"))
#peak locations
peaks.loc <- sort(pam(rep(pam(sample(x,round(length(x)/1000)),k=k)$medoids,5000),k=k)$medoids)
#peak heights using the ecdf 
e <- ecdf(x)
peaks.ecdf.height <- e(peaks.loc)
#peak heights using the density 
dens.fun <- function(x,...) {
    d <- density(x,...)
    return(splinefun(d$x,d$y))
}
d <- dens.fun(x)
peaks.dens.height <- d(peaks.loc)

if (!is.null(opt$plot)) {
png(opt$plot)
par(mfrow=c(2,1))
plot(density(x))
abline(v=peaks.loc)
abline(h=peaks.dens.height)
plot(ecdf(x))
abline(v=peaks.loc)
abline(h=peaks.ecdf.height)
dev.off()
}


cat(c('>fcsFile', paste(p, 'loc', 1:k, sep='.'), paste(p, 'ecdf.height', 1:k, sep='.'), paste(p, 'dens.height', 1:k, sep='.')), sep=',')
cat('\n')
cat('>', opt$fcs, peaks.loc, peaks.ecdf.height, peaks.dens.height, sep=',')
cat('\n')

