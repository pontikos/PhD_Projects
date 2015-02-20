#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

option_list <- list( 
make_option(c("--dir"), help = "dir to find fcsfiles to parse"),
make_option(c("--param", default='pstat5'), help="parameter [default %default]"),
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowCore, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))


flowset <- read.flowSet(path=opt$dir)


conc.col <- data.frame(conc=paste('_', c('0', '01', '10', '1000'), 'U.fcs', sep=''), col=c('black', 'blue', 'green', 'red'))
trans <- logicleTransform()

d <- list()
for (conc in conc.col$conc) {
    f <- flowset[[grep(conc, flowset@phenoData@data$name)]]
    d[[conc]] <- trans( f@exprs[,grep('stat5', f@parameters@data$desc, ignore.case=T)] )
}

pdf('plot.pdf')
plot(ecdf(d[[1]]), main=opt$dir)
sapply(1:length(d), function(i) lines(ecdf(d[[i]]), col=conc.col$col[i]))
dev.off()

