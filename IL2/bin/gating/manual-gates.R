#source('~nikolas/bin/FCS/transforms.R')
source('~nikolas/bin/FCS/fcs.R')

# the manual gates
# these are in fact ellipses derived from a CLR file exported from FlowJo
# I prefer dealing with CLR file (essentially a bitmap) rather than parsing
# the gate coordinates out of the worspace file

#manual gating

#fcs.data <- read.FCS('/chiswick/data/store/facs/Tony-FCS/PSTAT5/CD25/CD45RA/CD4/FOXP3/CB00165D_0U_2012-11-29.fcs', channels=c('SSCH','SSCW',"FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3", 'PSTAT5'), TRANS=logicleTransform(w=1))
#fcs.data <- read.FCS('/chiswick/data/store/facs/Tony-FCS/PSTAT5/CD25/CD45RA/CD4/FOXP3/CB00366X_0U_2012-11-07.fcs', channels=c('SSCH','SSCW',"FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3", 'PSTAT5'), TRANS=logicleTransform(w=1))


#CLR <- read.csv('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CLR/CB00165D_0U_2012-11-29.clr')
#load( '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/CB00165D_2012-11-29.RData' )

CLR.CELL.TYPES <- c("Lymphocytes", "Single Cells", "CD4", "Memory", "Memory Eff", "Memory Treg", "Naive", "Naive Eff", "Naive Treg")
CLR <- read.csv('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CLR/CB00366X_2012-11-07.clr')
colnames(CLR) <- CLR.CELL.TYPES


#load('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/transforms.RData')
load('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/transform-w1.RData')

#par(mfrow=c(3,7))
#for (n in colnames(fcs.data)) plot(density(fcs.data[,n]), main=n, xlim=range(fcs.data[,n][percentile.filter(fcs.data[,n])]))

load( '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/RData/pstat5-join/CB00366X_2012-11-07.RData' )

print(colnames(fcs.data))
print(names(transforms))
fcs.data <- applyTransforms(fcs.data, transforms)

fcs.data <- baseline.relative.pstat5(fcs.data)

e.lymphocytes <- classification.to.ellipse(fcs.data[,c('FSCA','SSCA')],CLR[,1])
e.singlecells <- classification.to.ellipse(fcs.data[,c('SSCH','SSCW')],CLR[,2])
e.cd4 <- classification.to.ellipse(fcs.data[,c('CD4','SSCA')],CLR[,3])
e.memory <- classification.to.ellipse(fcs.data[,c('CD45RA','SSCA')],CLR[,4])
e.memory.conv <- classification.to.ellipse(fcs.data[,c('CD25','FOXP3')],CLR[,5])
e.memory.tregs <- classification.to.ellipse(fcs.data[,c('CD25','FOXP3')],CLR[,6])
e.naive <- classification.to.ellipse(fcs.data[,c('CD45RA','SSCA')],CLR[,7])
e.naive.conv <- classification.to.ellipse(fcs.data[,c('CD25','FOXP3')],CLR[,8])
e.naive.tregs <- classification.to.ellipse(fcs.data[,c('CD25','FOXP3')],CLR[,9])

#evil hack
e.memory.tregs$Sigma<-e.naive.tregs$Sigma


# intersection of elliptical gates
g.cd4 <- list("Lymphocytes"=e.lymphocytes, "Single cells"=e.singlecells, "CD4"=e.cd4)
g.lymphocytes.singlecells.cd4.memory <- g.cd4
g.lymphocytes.singlecells.cd4.naive <- g.cd4
#
g.lymphocytes.singlecells.cd4.memory[["Memory"]] <- e.memory
g.memory.conv <- g.lymphocytes.singlecells.cd4.memory
g.memory.conv[['Memory Eff']] <- e.memory.conv
g.memory.tregs <- g.lymphocytes.singlecells.cd4.memory
g.memory.tregs[['Memory TReg']] <- e.memory.tregs
#
g.lymphocytes.singlecells.cd4.naive[["Naive"]] <- e.naive
g.naive.conv <- g.lymphocytes.singlecells.cd4.naive
g.naive.conv[['Naive Eff']] <- e.naive.conv
g.naive.tregs <- g.lymphocytes.singlecells.cd4.naive
g.naive.tregs[['Naive TReg']] <- e.naive.tregs
#
gates <- list('Memory Eff'=g.memory.conv, 'Memory Treg'=g.memory.tregs, 'Naive Eff'=g.naive.conv, 'Naive Treg'=g.naive.tregs)




