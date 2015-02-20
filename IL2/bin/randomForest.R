#!/usr/bin/env Rscript
source('~nikolas/bin/FCS/fcs.R')
library(randomForest) 
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("flowCore"))

option_list <- list( 
make_option(c("--in.file"), help = "List of files to join to file. First file is taken as reference file.")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

load( opt$in.file )

print(dim(fcs.data))
print(clr.file <- file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR/',basename(opt$in.file)))
print(load(clr.file))
print(dim(fcs.data <- fcs.data[which(as.logical(CLR[,'CD4'])),]))

#print(load('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/transforms.RData'))
print(load('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/transform-w1.RData'))
fcs.data <- applyTransforms(fcs.data, transforms)
fcs.data <- baseline.relative.pstat5(fcs.data)

print(head(fcs.data))

#cd4.lymph <- as.data.frame(fcs.data[which(as.logical(CLR[,'CD4'])),])
print(form <- formula(paste('diff.PSTAT5.4','~',paste(grep('PSTAT5',colnames(fcs.data),value=TRUE,invert=TRUE),collapse='+'))))
#would like to use proximity but too large
RF <- randomForest(form, data=fcs.data, proximity=FALSE)

head(getTree(RF,1))

save(RF, file=file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/RF/',basename(opt$in.file)))

print(importance(RF))

