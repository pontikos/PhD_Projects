library(iterators)
source('~nikolas/bin/FCS/fcs.R')



setwd('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/RData')

MARKERS <- c('CD25','CD3','CD4','CD45RA','CD56','CD8','FOXP3')

load('~/dunwich/Projects/IL2/transforms.RData')

for (f in list.files(pattern='.*.RData')) {
    print(load(f))
    fcs.data <- applyTransforms(fcs.data,transforms)
    break
}



