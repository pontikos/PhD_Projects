source('~/bin/FCS/fcs.R')
channels <- c('SSCA','FSCA','CD4','CD45RA')
load("~/dunwich/Projects/IL2/transforms.RData")

#d <- read.FCS('/chiswick/data/store/facs/Marcin/Marcin_CD25_CBR200/CAD100_2008aug05_Treg/CAD100_2008aug05_Treg_I009114F_009.fcs',channels=channels,TRANS=NULL)
d2 <- read.FCS('/chiswick/data/store/calliope-ARCHIVE/FOR_VP/Calli_200/CAD200_Treg/CAD100_050808_Treg/CAD100_2008aug05_Treg_I009114F_009.fcs',channels=channels,TRANS=NULL) 
d3 <- read.FCS('/chiswick/data/store/calliope-ARCHIVE/FOR_VP/Calli_200/CAD200_Treg/CAD70_060508_Treg/CAD70_2008may06_Treg_I007853K_015.fcs',channels=channels,TRANS=NULL)  
d4 <- read.FCS('/chiswick/data/store/calliope-ARCHIVE/FOR_VP/Calli_200/CAD200_Treg/CAD67_160408_Treg/CAD67_2008apr16_Treg_I007693L_009.fcs',channels=channels,TRANS=NULL)


par(mfrow=c(2,2))
smoothPlot(d[,c('FSCA','SSCA')],outliers=FALSE)
smoothPlot(d2[,c('FSCA','SSCA')],outliers=FALSE)
smoothPlot(d3[,c('FSCA','SSCA')],outliers=FALSE)
smoothPlot(d4[,c('FSCA','SSCA')],outliers=FALSE)
smoothPlot(d5[,c('FSCA','SSCA')],outliers=FALSE)


transforms <- list()
for (chan in c('CD4','CD127','CD25','CD45RA')) {
  transforms[[chan]] <- log10
}
for (chan in c('SSCA','FSCA')) {
  transforms[[chan]] <- identity
}

channels <- c('SSCA','FSCA','CD4','CD127','CD25','CD45RA')
#d <- read.FCS('/chiswick//data/store/calliope-ARCHIVE/FOR_VP/Calli_200/CAD200_Treg/CAD116_091008_Treg/CAD116_2008oct09_Treg_I009546A_013.fcs',channels=channels,TRANS=transforms)
dim(d <- read.FCS('/chiswick//data/store/calliope-ARCHIVE/FOR_VP/Calli_200/CAD200_Treg/CAD100_050808_Treg/CAD100_2008aug05_Treg_I009114F_009.fcs',channels=channels,TRANS=transforms))
dim(d <- d[-which(d==0,arr.ind=T)[,1],])

figure.labels <- iter(paste(letters,')',sep=''))
par(mfrow=c(2,2))
 #a) Lymphocytes
channels <- c('FSCA','SSCA')
smoothPlot(d[,channels],outliers=FALSE)
print(nrow(d))
title(paste(nextElem(figure.labels),"Lymphocytes"), adj=0)
#G <- locator(type='l')
#G <- structure(list(FSCA = c(928.414544519474, 911.84663262999, 887.547476125493,  728.064889525607, 602.640966626337, 451.789841302946, 384.198114977502,  263.228164733921, 225.871373397841, 244.935541602152, 282.886943462597,  295.912455198282, 602.640966626337, 818.466630350764, 856.15287549311,  895.574381451552, 903.673881439695, 928.414544519474, 928.414544519474 ), SSCA = c(87, 70, 54, 37, 26, 17, 15, 19, 68, 104, 140, 152,  167, 176, 172, 155, 150, 131, 87)), .Names = c("FSCA", "SSCA" ), row.names = c(NA, -19L), class = "data.frame")
G <- structure(list(FSCA= c(641.285658710046, 515.844912141632, 390.404165573217,  411.310956667953, 822.477848197755, 1031.54575914511, 1393.93013812053,  1582.09125797315, 1533.30874541877, 1373.0233470258, 954.887525131081,  759.757474913548, 571.596355060927, 446.155608492513), SSCA= c(19.9784114636972,  66.1633018735767, 140.572291978382, 220.112936573175, 248.337036268101,  256.034518003081, 256.034518003081, 140.572291978382, 55.8999928936035,  40.5050294236437, 17.412584218704, 17.412584218704, 35.3733749336571,  86.6899198335231)), .Names = c("FSCA","SSCA"))
#lines(d[f3,channels][c(ch,ch[1]),],col='red',lwd=4)
#lines(d[f3,channels][c(ch,ch[1]),],col='red',lwd=4)
d <- d[as.logical(in.poly <- point.in.polygon(d[,'FSCA'],d[,'SSCA'],G$FSCA,G$SSCA)),]  
print(nrow(d))
100*prop.table(table(as.logical(in.poly)))
ch <- chull(d[,c('FSCA','SSCA')])
#d <- d[f3,][c(ch,ch[1]),]
lines(d[c(ch,ch[1]),channels],col='red',lwd=4)
#points(d[f3,channels])
#b) CD4+
channels <- c('CD4','SSCA')
smoothPlot(d[,channels],outliers=FALSE)
title(paste(nextElem(figure.labels),"CD4+"), adj=0)
#G <- locator(type='l')
#G <- structure(list(x = c(2.11952316242512, 2.27959888532937, 2.35747356133685,  2.26661977266146, 2.05462759908556, 1.89022550529201, 1.86426727995618,  2.12384953331443, 2.32718896511172), y = c(213.215637015057,  174.76622183342, 56.8546819430669, -9.79097103843732, -35.423914492862,  20.9685611068723, 128.626923615456, 210.652342669615, 144.006689688111 )), .Names = c("x", "y"))
#G <- structure(list(x = c(1.54090714954455, 1.0153243179039, 0.924358058581476,  1.04059272327124, 1.52574610632415, 1.93509427327504, 2.07154366225868,  2.0513289379648, 1.98057740293626, 1.65714181423431), y = c(385.267445156431,  383.109414739122, 289.594763322387, 242.118094141583, 239.240720251837,  245.714811503765, 295.349511101878, 334.194058613445, 371.599919180139,  384.188429947776)), .Names = c("x", "y"))
G <- structure(list(CD4 = c(1.92888659641913, 2.04974755372363, 2.20641916504426,  2.35861444461289, 2.4123304256371, 2.47499907016536, 2.46157007490931,  2.39890143038105, 2.2646114778205, 1.97812624569133, 1.79459664385859,  1.74088066283437, 1.8483126248828, 1.8706942836429), SSCA = c(125.865951539838,  155.652850257628, 167.071161432781, 147.709677266217, 131.326882971433,  113.951192052722, 83.6678450229689, 56.859636176958, 46.4342216257315,  43.4555317539525, 56.859636176958, 75.2282237195951, 106.008019061311,  121.39791673217)), .Names = c("CD4", "SSCA"))
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

smoothPlot(d[,c('CD4','FSCA')])
#G <- locator(type='l')
G <- structure(list(x = c(2.18646095657097, 1.96304926329938, 1.89860358254796,  2.12201527581955, 2.3153523180738, 2.34972334780789, 2.25949939475591,  2.1606826842704), y = c(935.76050813135, 797.686550245687, 457.104120794386,  208.570996600193, 328.235093434434, 622.792870257181, 871.325994451374,  954.170369182771)), .Names = c("x", "y"))
f <- point.in.polygon(d[,'CD4'],d[,'FSCA'],G$x,G$y)
smoothPlot(d[,c('CD4','SSCA')])
#G <- locator(type='l')
G <- structure(list(x = c(2.13490441196983, 2.006013050467, 1.91149271869825,  1.9501601271491, 2.10912613966926, 2.28527766705648, 2.32394507550733,  2.33683421165761, 2.2380175011721, 2.13490441196983, 2.06616235250165 ), y = c(179.059700653282, 143.138119223376, 73.8607836085566,  17.412584218704, 4.58344799373746, 4.58344799373746, 58.4658201385968,  132.874810243403, 166.230564428315, 179.059700653282, 173.928046163295 )), .Names = c("x", "y"))
f2 <- point.in.polygon(d[,'CD4'],d[,'SSCA'],G$x,G$y)
f3 <- as.logical(f) & as.logical(f2)
      

