source('~nikolas/Projects/IL2/bin/common.R')
library(pls)

BASE.DIR <- base.dir

# partial least squares regression 
individual.date <- do.call('rbind',strsplit(gsub('.RData','',list.files(file.path(base.dir,'RData','pstat5-join'))),'_')) 

######### Lymphocytes
source('~nikolas/Projects/IL2/bin/PLSR/thesis-lymphocytes.R')

######### Non-lymphocytes
source('~nikolas/Projects/IL2/bin/PLSR/thesis-nonlymphocytes.R')


