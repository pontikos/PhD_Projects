#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(mixtools))
suppressPackageStartupMessages(library(cluster))
source('~nikolas/bin/FCS/fcs.R')

option_list <- list( 
make_option(c("--fcsFile"), help = ""),
make_option(c("--channel"), help = ""),
make_option(c('--RData'), help='')
) 
OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

X <- read.FCS(opt$fcsFile, channel=opt$channel, TRANS=log10)[,1]

res.pam <- pam(X,2)
res <- list( mu=res.pam$medoids[,1], lambda=as.numeric(prop.table(table(res.pam$clustering))), sigsqrd=as.numeric(by(X,res.pam$clustering,var)) )
m <- mixtools::normalmixEM2comp( X, mu=res$mu, sigsqrd=res$sigsqrd, lambda=res$lambda ) 

write(m, file=opt$RData) 


