# CD45RA/FOXP3 spillover

load("~/dunwich/Projects/IL2/transforms.RData")
#scale scatter to be on a similar scale as fluorescence
transforms[['SSCA']] <- function(x) 10*(x-min(x))/(max(x)-min(x))
transforms[['FSCA']] <- function(x) 10*(x-min(x))/(max(x)-min(x))
 
old.comp <- read.csv('~/Projects/IL2/oldcompmatrix.csv',check.names=FALSE)
comp <- read.csv('~/Projects/IL2/compmatrix.csv',check.names=FALSE)

#3 files from the same batch
#I021196N_CB00086S_0U.fcs
#I021228Y_CB01482K_0U.fcs
#I021212F_CB01481J_0U.fcs


#f1 <- flowCore::read.FCS('/chiswick/data/store/facs/Tony/Tony/180912_IL2_sens_T1D_pSTAT5_TC2_NKCD/I021196N_CB00086S_0U.fcs')
d1 <- read.FCS('/chiswick/data/store/facs/Tony/Tony/180912_IL2_sens_T1D_pSTAT5_TC2_NKCD/I021196N_CB00086S_0U.fcs', channels=CORE.MARKERS, spillover=old.comp)
d2 <- read.FCS('/chiswick/data/store/facs/Tony/Tony/180912_IL2_sens_T1D_pSTAT5_TC2_NKCD/I021196N_CB00086S_0U.fcs', channels=CORE.MARKERS, spillover=comp)

load('~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/CLR/CB00086S_0U_2012-09-18.RData')

figure.labels <- iter(paste(letters,')',sep=''))
par(mfrow=c(2,2))
smoothPlot(d1[,c('CD45RA','FOXP3')],posteriors=CLR[,c('Memory Eff','Memory Treg','Naive Eff','Naive Treg')],ellipses=FALSE,ellipses.lwd=2)
title(nextElem(figure.labels), adj=0)
smoothPlot(d2[,c('CD45RA','FOXP3')],posteriors=CLR[,c('Memory Eff','Memory Treg','Naive Eff','Naive Treg')])
title(nextElem(figure.labels), adj=0)
smoothPlot(d1[,c('CD8','FOXP3')])
title(nextElem(figure.labels), adj=0)
smoothPlot(d2[,c('CD8','FOXP3')])
title(nextElem(figure.labels), adj=0)



#f2 <- flowCore::read.FCS('/chiswick/data/store/facs/Tony/Tony/180912_IL2_sens_T1D_pSTAT5_TC2_Treg/I021196N_CB00086S_0U.fcs')
d2 <- read.FCS('/chiswick/data/store/facs/Tony/Tony/180912_IL2_sens_T1D_pSTAT5_TC2_Treg/I021196N_CB00086S_0U.fcs',channels=grep('CD8|CD3',CORE.MARKERS,invert=T,value=T))

d3 <- applyTransforms(read.FCS('/chiswick/data/store/facs/Tony/Tony/22_August_2102_DGAP_NKCD8/KA866533H_0_U.fcs', channels=CORE.MARKERS, TRANS=NULL),transforms)

figure.labels <- iter(paste(letters,')',sep=''))
par(mfrow=c(2,2))
smoothPlot(d1[,c('CD45RA','FOXP3')])
title(nextElem(figure.labels), adj=0)
smoothPlot(d2[,c('CD45RA','FOXP3')])
title(nextElem(figure.labels), adj=0)
smoothPlot(d3[,c('CD45RA','FOXP3')])
title(nextElem(figure.labels), adj=0)
 
smoothPlot(d3[,c('SSCA','FSCA')])

