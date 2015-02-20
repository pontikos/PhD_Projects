# these are the individuals I will use as examples in the final section to see if
# I can identify the same cells in other samples.

# According to Tony's spreadsheet, panel 2 is the one containing:
# pSTAT5 - AF 4888
# CD3 - PCPCy5.5
# CD56 - bv421
# CD4 - AF700
# FOXP3 - PE
# CD45RA - PE Cy7
# CD25 - APC
# CD8 - Bv 605

panel <-list("FSC-A"="FSCA", "FSC-H"="FSCH", "FSC-W"="FSCW", "SSC-A"="SSCA", "SSC-H"="SSCH", "SSC-W"="SSCW", "Alexa Fluor 488-A"="PSTAT5", "PerCP-Cy5-5-A"="CD3", "APC-A"="CD25", "Alexa Fluor 700-A"="CD4", "APC-Cy7-A"="APC-Cy7-A", "Pacific Blue-A"="CD56", "Qdot 605-A"="CD8", "PE YG-A"="FOXP3", "PE-Cy7 YG-A"="CD45RA", "Time"="Time") 
f <- '/dunwich/scratch/nikolas/FCS.Tony/pSTAT5_DGAP_KM00782Z/pSTAT5_DGAP_KM00782Z_0U.fcs'
f1 <- flowCore::read.FCS(f)
d <- f1@exprs
colnames(d) <- panel[colnames(f1)]
d <- d[,CORE.MARKERS]
d <- applyTransforms(d,transforms)
 

figure.labels <- iter(paste(letters,')',sep=''))
par(mfrow=c(2,2))
#a) Lymphocytes
channels <- c('FSCA','SSCA')
smoothPlot(d[,channels],outliers=FALSE)
print(nrow(d))
title(paste(nextElem(figure.labels),"Lymphocytes"), adj=0)
G <- structure(list(FSCA=c(101222.172487118, 76692.169242046, 66434.167885016, 64650.1676490108, 86504.1705400747, 114602.174257157, 119062.17484717,  107912.173372137), SSCA=c(39527.3275550547, 35859.2021732815,  28522.9514097351, 13850.4498826422, 4221.62075548753, 16143.0282462505,  31274.045446065, 37234.7491914464)), .Names = c("FSCA", "SSCA"))
in.poly <- point.in.polygon(d[,'FSCA'],d[,'SSCA'],G$FSCA,G$SSCA)
d <- d[as.logical(in.poly)),]  
print(nrow(d))
100*prop.table(table(as.logical(in.poly)))
ch <- chull(d[,channels]))
#d <- d[f3,][c(ch,ch[1]),]
lines(d[c(ch,ch[1]),channels],col='red',lwd=4)
#points(d[f3,channels])
#b) CD4+
channels <- c('CD4','SSCA')
smoothPlot(d[,channels],outliers=FALSE)
title(paste(nextElem(figure.labels),"CD4+"), adj=0)
G <- structure(list(CD4= c(2.61208110909558, 2.42158681605733, 2.37197892724529,  2.31244946067083, 2.25887294075382, 2.29062198959353, 2.33626124730061,  2.46127312710697, 2.65176742014522, 2.76090477553172, 2.76090477553172,  2.75495182887427, 2.71328120227216, 2.62993994906792, 2.58628500691332 ), SSCA= c(28332.3360029655, 27254.8870739755, 25682.3940424767,  23178.0532886081, 18810.0170900002, 16014.4739228911, 13772.215340939,  13073.3295491618, 13131.5700318099, 14791.4237872809, 18664.4158833799,  22479.1674968308, 25362.0713879121, 27080.1656260312, 27458.7287632439 )), .Names = c("CD4", "SSCA"))
d <- d[as.logical(in.poly <- point.in.polygon(d[,channels[1]],d[,channels[2]],G$CD4,G$SSCA)),]
print(nrow(d))
100*prop.table(table(as.logical(in.poly)))
ch<-chull(d[,channels])
lines(d[c(ch,ch[1]),channels],col='red',lwd=4)
#c) non-Tregs
channels <- c('CD25','CD127')
smoothPlot(d[,channels],outliers=FALSE)
title(paste(nextElem(figure.labels),"non-Tregs"), adj=0)
#G <- locator(type='l')
G <- structure(list(CD25 = c(-0.0615496613501853, 0.464992919811841,  1.21249533414008, 1.68732391179512, 2.15274994335798, 2.29378813474067,  2.33139831910938, 2.30319068083284, 1.786050645763, 1.24070297241661,  0.859899855683361, 0.389772551074409, 0.0512808917559632, -0.103861118764991,  -0.0991598457189014), CD127 = c(0.839362777435828, 0.95498778496685,  1.01922390026186, 1.04491834637987, 1.15626094622456, 1.49028874575862,  1.79862209917468, 1.96563599894171, 1.96563599894171, 1.96135359125538,  1.94850636819638, 1.95278877588271, 1.95278877588271, 1.96563599894171,  0.916446115789843)), .Names = c("CD25", "CD127"))
d <- d[as.logical(in.poly <- point.in.polygon(d[,channels[1]],d[,channels[2]],G$CD25,G$CD127)),]
100*prop.table(table(as.logical(in.poly)))
print(nrow(d))
ch<-chull(d[,channels])
lines(d[c(ch,ch[1]),channels],col='red',lwd=4)
#d) Memory / Naive
channels <- c('CD25','CD45RA')
smoothPlot(d[,channels],outliers=FALSE)
title(paste(nextElem(figure.labels),"Memory / Naive"), adj=0)
abline(h=c(.75,1),col='red',lwd=4)
cat('prop memory', 100*prop.table(table(d[,'CD45RA']<.75)), '\n')
abline(v=1,col='red',lwd=4)
cat('prop naive cd25pos', 100*prop.table(table(d[,'CD45RA']>1 & d[,'CD25']>1)), '\n')
#G <- locator(type='l')



X <- d[,CORE.FMARKERS]
X<-5*X[,'CD56']
X<-2*X[,'CD45RA']
res <- kmeans(X , 6)
plotClusters(d[,CORE.FMARKERS],classification=res$cluster,chulls=FALSE)

#, classification=d[,'CD56']>1.5)
#G <- locator(type='l')
#lines(x=g$FSCA,y=g$SSCA,lwd=2,col="red")
#d<-d[which(d[,'CD56']>2.6 & d[,'CD3']<2.5),] 
#plotClusters(d[,c('CD3','CD4','CD25','CD56')])
 

# non-lymphocytes
smoothPlot(d[,c('FSCA','SSCA')])
#G <- locator(type='l')
G <- structure(list(FSCA= c(189273.430370576, 208261.662097808, 219839.852175389,  212429.810525737, 185105.281942647, 165653.922612311, 165653.922612311,  184642.154339543, 204093.513669879), SSCA= c(61338.4455984501,  61338.4455984501, 48892.6820326653, 34891.1980211574, 24001.1549010956,  27112.5957925418, 45781.2411412191, 60301.2986346348, 59782.725152727 )), .Names = c("FSCA", "SSCA"))
#lines(x=g$FSCA,y=g$SSCA,lwd=2,col="red")
#100*prop.table(table(as.logical(in.poly)))
d <- d[as.logical(in.poly <- point.in.polygon(d[,'FSCA'],d[,'SSCA'],G$FSCA,G$SSCA)),]  
ch <- chull(d[,c('FSCA','SSCA')])
#d <- d[f3,][c(ch,ch[1]),]
lines(d[c(ch,ch[1]),names(G)],col='red',lwd=4)
plotClusters(d[,c('CD3','CD4','CD25','CD56')])
#, classification=d[,'CD56']>1.5)

 


#
# pSTAT5 - AF 4888
# CD3 - PCPCy5.5
# CD56 - bv421
# CD4 - AF700
# FOXP3 - PE
# CD45RA - PE Cy7
# CD25 - APC
# CD8 - Bv 605


# Individuals with panel 2 are:
# KM00782Z  23.07.2012
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KM00782Z/23_July_2102_pSTAT5_D-GAP-KM00782Z_1000U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KM00782Z/23_July_2102_pSTAT5_D-GAP-KM00782Z_10U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KM00782Z/23_July_2102_pSTAT5_D-GAP-KM00782Z_0_1U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KM00782Z/23_July_2102_pSTAT5_D-GAP-KM00782Z_0_U.fcs
read.FCS('/dunwich/scratch/nikolas/FCS.Tony/pSTAT5_DGAP_KM00782Z/pSTAT5_DGAP_KM00782Z_01U.fcs')
read.FCS('/dunwich/scratch/nikolas/FCS.Tony/pSTAT5_DGAP_KM00782Z/pSTAT5_DGAP_KM00782Z_10U.fcs')
read.FCS('/dunwich/scratch/nikolas/FCS.Tony/pSTAT5_DGAP_KM00782Z/pSTAT5_DGAP_KM00782Z_1000U.fcs')
# KM00783A	23.07.2012
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KM00783A/23_July_2102_pSTAT5_D-GAP-KM00783A_0_U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KM00783A/23_July_2102_pSTAT5_D-GAP-KM00783A_0_1U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KM00783A/23_July_2102_pSTAT5_D-GAP-KM00783A_10U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KM00783A/23_July_2102_pSTAT5_D-GAP-KM00783A_1000U.fcs
# KM00784B	23.07.2012
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KM00784B/23_July_2102_pSTAT5_D-GAP-KM00784B_0_1U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KM00784B/23_July_2102_pSTAT5_D-GAP-KM00784B_10U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KM00784B/23_July_2102_pSTAT5_D-GAP-KM00784B_1000U.fcs
# KA866533H KM00853B  22.08.2012
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KA866533H/22_August_2102_DGAP_Treg-KA866533H_10U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KA866533H/22_August_2102_DGAP_Treg-KA866533H_0_1U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KA866533H/22_August_2102_DGAP_NKCD8-KA866533H_10U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KA866533H/22_August_2102_DGAP_NKCD8-KA866533H_0_1U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KA866533H/22_August_2102_DGAP_NKCD8-KA866533H_1000U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KA866533H/22_August_2102_DGAP_NKCD8-KA866533H_0_U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KA866533H/22_August_2102_DGAP_Treg-KA866533H_0.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KA866533H/22_August_2102_DGAP_Treg-KA866533H_1000U.fcs
/dunwich/scratch/nikolas/FCS.Tony/DGAP_NKCD8_KA866533H/DGAP_NKCD8_KA866533H_1000U.fcs
/dunwich/scratch/nikolas/FCS.Tony/DGAP_NKCD8_KA866533H/DGAP_NKCD8_KA866533H_0U.fcs
/dunwich/scratch/nikolas/FCS.Tony/DGAP_NKCD8_KA866533H/DGAP_NKCD8_KA866533H_01U.fcs
/dunwich/scratch/nikolas/FCS.Tony/DGAP_NKCD8_KA866533H/DGAP_NKCD8_KA866533H_10U.fcs
/dunwich/scratch/nikolas/FCS.Tony/DGAP_Treg_KA866533H/DGAP_Treg_KA866533H_1000U.fcs
/dunwich/scratch/nikolas/FCS.Tony/DGAP_Treg_KA866533H/DGAP_Treg_KA866533H_01U.fcs
/dunwich/scratch/nikolas/FCS.Tony/DGAP_Treg_KA866533H/DGAP_Treg_KA866533H_10U.fcs
# KA118577M KM00854C	22.08.2012
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KA118577M/22_August_2102_DGAP_NKCD8-KA118577M_1000U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KA118577M/22_August_2102_DGAP_Treg-KA118577M_0.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KA118577M/22_August_2102_DGAP_Treg-KA118577M_1000U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KA118577M/22_August_2102_DGAP_NKCD8-KA118577M_0_1U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KA118577M/22_August_2102_DGAP_NKCD8-KA118577M_10U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KA118577M/22_August_2102_DGAP_NKCD8-KA118577M_0_U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KA118577M/22_August_2102_DGAP_Treg-KA118577M_0_1U.fcs
/dunwich/scratch/nikolas/Projects/IL2/clean.Tony/KA118577M/22_August_2102_DGAP_Treg-KA118577M_10U.fcs
/dunwich/scratch/nikolas/FCS.Tony/DGAP_Treg_KA118577M/DGAP_Treg_KA118577M_1000U.fcs
/dunwich/scratch/nikolas/FCS.Tony/DGAP_Treg_KA118577M/DGAP_Treg_KA118577M_01U.fcs
/dunwich/scratch/nikolas/FCS.Tony/DGAP_Treg_KA118577M/DGAP_Treg_KA118577M_0U.fcs
/dunwich/scratch/nikolas/FCS.Tony/DGAP_Treg_KA118577M/DGAP_Treg_KA118577M_10U.fcs
/dunwich/scratch/nikolas/FCS.Tony/DGAP_NKCD8_KA118577M/DGAP_NKCD8_KA118577M_0U.fcs
/dunwich/scratch/nikolas/FCS.Tony/DGAP_NKCD8_KA118577M/DGAP_NKCD8_KA118577M_10U.fcs
/dunwich/scratch/nikolas/FCS.Tony/DGAP_NKCD8_KA118577M/DGAP_NKCD8_KA118577M_01U.fcs
/dunwich/scratch/nikolas/FCS.Tony/DGAP_NKCD8_KA118577M/DGAP_NKCD8_KA118577M_1000U.fcs





 # this is the individual I use as an example in the IL2 Chapter
individual <- 'CB00086S'
date <- '2012-09-18' 
base.dir <- '~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/' 
 
#ok I021445J CB01486P 26.09.2012 T1D
individual <- 'CB00618W'
date <- '2012-05-28' 
  
#witthout the decimal dot in 01U
DOSES_ <- c( '0U', '01U', '10U', '1000U')
# with the decimal dot in 0.1U
DOSES <- c( '0U', '0.1U', '10U', '1000U')


CORE.FMARKERS <- c('CD25','CD3','CD4','CD45RA','CD56','CD8','FOXP3')
CORE.MARKERS <- c('SSCA','FSCA','CD25','CD3','CD4','CD45RA','CD56','CD8','FOXP3')

FILES <- sprintf('%s.RData', paste(individual, DOSES_, date, sep='_'))

load("~/dunwich/Projects/IL2/transforms.RData")
#scale scatter to be on a similar scale as fluorescence
transforms[['SSCA']] <- function(x) 10*(x-min(x))/(max(x)-min(x))
transforms[['FSCA']] <- function(x) 10*(x-min(x))/(max(x)-min(x))

load(FILES[[1]])

#CB01090J 28.05.2012
#CB01110F 28.05.2012
#CB01325P 28.05.2012
#CB01388H 28.05.2012
s <- '/chiswick/data/store/facs/Tony-FCS/chiswick-facs-Tony/28_May_2102_pSTAT5_DGAP/CB00618W_0.fcs'

s <- '/chiswick/data/store/facs/Tony/Tony/28_May_2102_pSTAT5_DGAP/CB01090J_1000.fcs'

f <- read.FCS(s)
f <- read.FCS(s,channels=CORE.MARKERS,TRANS=transforms)

smoothPlot()


