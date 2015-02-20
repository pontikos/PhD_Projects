#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

option_list <- list( 
make_option(c('-f', "--fcs"), help = "fcsfile to parse")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowCore, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))

trans <- logicleTransform()

suppressWarnings(read.FCS(opt$fcs, alter.names=TRUE)) -> f

cd25 <- trans( f@exprs[,grep('cd25', f@parameters@data$desc, ignore.case=T)] )
pstat5 <- trans( f@exprs[,grep('pstat5', f@parameters@data$desc, ignore.case=T)] )


pdf.file <- gsub('.fcs', '.pdf', opt$fcs)
print(pdf.file)
pdf(pdf.file)
plot(ecdf(cd25), main=basename(opt$fcs))
lines(ecdf(pstat5), col='blue')
dev.off()


