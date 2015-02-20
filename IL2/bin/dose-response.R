#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

option_list <- list( 
make_option(c('-d', "--dir"), help = "dir to find fcsfiles to parse"),
make_option(c('-p', "--param"), default='pstat5', help="parameter [default %default]")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowCore, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))


flowset <- read.flowSet(path=opt$dir, pattern='*.fcs')


style <- data.frame(conc=as.character(paste('_', c('0', '01', '10', '1000'), 'U.fcs', sep='')),
                    col=c('black', 'blue', 'green', 'red'),
                    lwd=seq(.5, by=.75, length=4),
                    stringsAsFactors=FALSE)

trans <- logicleTransform()

d <- list()
for (conc in style$conc) {
    f <- flowset[[grep(conc, flowset@phenoData@data$name)]]
    d[[conc]] <- trans( f@exprs[,grep(opt$param, f@parameters@data$desc, ignore.case=T)] )
}

baseline.conc <- d[['_0U.fcs']]
baseline.conc.ecdf <- ecdf(baseline.conc)
high.conc <- d[['_1000U.fcs']]
high.conc.ecdf <- ecdf(high.conc) 

x <- seq(min(baseline.conc), max(baseline.conc), .1)
ecdf.diff <- function(x) return(high.conc.ecdf(x)-baseline.conc.ecdf(x))
area <- sum(abs(ecdf.diff(x)))
cat('>>', opt$dir, area, '\n', sep=',')



plot.dose.response <- function() {
    pdf(file.path(opt$dir,'plot.pdf'))
    plot(ecdf(d[['_0U.fcs']]), main=basename(opt$dir), xlab=opt$param, ylab='', lwd=0, col='white')
    sapply(1:length(d), function(i) lines(
                                          ecdf(d[[style$conc[i]]]),
                                          lwd=style$lwd[i],
                                          col=style$col[i]
                                          #lty=style$lty[i],
                                          #pch=style$pch[i]
                                          )) 
    lines(x, abs(ecdf.diff(x)), lty=1, lwd=.25) 
    #legend(0,1, style$conc, cex=.8, col=style$col, pch=style$pch, lty=style$lty)
    legend(x="bottomright", legend=style$conc, col=style$col, bty='n')
    dev.off()
}


