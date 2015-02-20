#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))

option_list <- list( 
make_option(c("--in.file"), default=NULL, help = ".RData CLR file")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

file.name <- opt$in.file
print(basename(file.name))
load(file.name)

i <- which((CLR[,'Memory']==1) & (CLR[,'Naive']==1))
CLR[i,'Naive'] <- 0
CLR[i,'Memory'] <- 0

i <- which((CLR[,'Naive Eff']==1) & (CLR[,'Naive Treg']==1))
CLR[i,'Naive Eff'] <- 0
CLR[i,'Naive Treg'] <- 0

i <- which((CLR[,'Memory Eff']==1) & (CLR[,'Memory Treg']==1))
CLR[i,'Memory Eff'] <- 0
CLR[i,'Memory Treg'] <- 0

save(CLR, file=file.name)

