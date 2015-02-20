library(pls)
library(iterators)
library(flowCore)
source('~nikolas/bin/FCS/fcs.R')
source('~nikolas/Projects/IL2/bin/common.R')

# partial least squares regression

# this individual looks ok
individual <- 'CB00086S'
day <- '2012-09-18'
BASE.DIR <- '~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5'
individual.date <- do.call('rbind',strsplit(gsub('.RData','',list.files(file.path(BASE.DIR,'RData','pstat5-join'))),'_')) 
setwd(BASE.DIR) 
#load('~/dunwich/Projects/IL2/transforms.RData')

#### All lymphocytes
print(f <- file.path(BASE.DIR,'RData','pstat5-join', sprintf('%s_%s.RData',individual,day)))
print(load(f))
fcs.data <- baseline.relative.pstat5(fcs.data)
#fcs.data <- applyTransforms(fcs.data,transforms)
print(load(file.path(BASE.DIR,'CLR','CB00086S_0U_2012-09-18.RData')))
print(dim(fcs.data <- fcs.data[as.logical(CLR[,'Single cells']),]))
MARKERS <- c('CD25','CD3','CD4','CD45RA','CD56','CD8','FOXP3')

#
fm <- as.formula(paste(sprintf('diff.PSTAT5.%d',4), paste(MARKERS, collapse='+'), sep='~'))
plsr(fm, data=as.data.frame(fcs.data),scale=FALSE)->p
Y <- p$scores[,1:2]
smoothPlot(Y, outliers=TRUE)

#g<-locator(type='l')
#g <- structure(list(x = c(-0.0103174435677124, -0.1599450850475, -0.339498254823246, -0.504088660451012, -0.548976952894949, -0.53401418874697, -0.474163132155055, -0.369423783119203, -0.264684434083352, -0.1599450850475, -0.115056792603564, -0.100094028455585), y = c(-1.62767441406014, -1.59984835807974, -1.73897863798175, -1.98941314180538, -2.2120215896486, -2.35115186955062, -2.35115186955062, -2.32332581357022, -2.2120215896486, -1.98941314180538, -1.85028286190336, -1.59984835807974)), .Names = c("x", "y"))
g <- structure(list(x = c(0.549386221184538, 0.299526136195612, 0.103689853366454, 0.0834309275565415, 0.21173745768599, 0.468350517944886, 0.542633245914567), y = c(-1.4566307230717, -1.58756893468617, -1.92800828488379, -2.16369706578983, -2.22916617159706, -1.92800828488379, -1.4566307230717)), .Names = c("x", "y"))
lines(g$x,g$y,col='lightblue')
j <- as.logical(point.in.polygon(Y[,1],Y[,2],g$x,g$y))

#g<-locator(type='l')
g <- structure(list(x = c(-0.444233868802116, -0.92602314405666, -1.55234920188757, -1.27933527924333, -0.861784574022721, -0.492412796327571), y = c(-0.163647369673267, 0.787426504230537, -0.269322244551468, -1.07949628528434, -1.1499462018698, -0.128422411380534)), .Names = c("x", "y"))
lines(g$x,g$y,col='pink')
k <- as.logical(point.in.polygon(Y[,1],Y[,2],g$x,g$y))

#g<-locator(type='l')
g <- structure(list(x = c(-0.701188148937873, -0.379995298768177, 0.439046469164547, 0.053615048960912, -0.701188148937873), y = c(0.857876420816004, 0.0829273383758669, 1.03400121227967, 1.70327541984161, 0.998776253986938)), .Names = c("x", "y"))
#lines(g$x,g$y,col='lightblue')
lines(g$x,g$y,col='purple')
l <- as.logical(point.in.polygon(Y[,1],Y[,2],g$x,g$y))

#identification of pink and blue high density clusters
print(dim(clr <- CLR[as.logical(CLR[,'Single cells']),]))
clr <- cbind(clr,lightblue=as.numeric(j))
clr <- cbind(clr,pink=as.numeric(k))
clr <- cbind(clr,purple=as.numeric(l))

pdf('~/Thesis/figures/plsr-lymphocytes.pdf',width=10,height=10)
figure.labels <- iter(paste(letters,')',sep=''))
par(mfrow=c(2,2),mai=c(.75,.75,.5,.1))
for (i in 2:4) {
    fm <- as.formula(paste(sprintf('diff.PSTAT5.%d',i), paste(MARKERS, collapse='+'), sep='~'))
    plsr(fm, data=as.data.frame(fcs.data),scale=FALSE)->p
    # TODO: get the loadings and the variance explained
    print(p)
    Y <- p$scores[,1:2]
    xlab <- 'Comp 1'
    #xlab <- ''
    ylab <- 'Comp 2'
    #ylab <- ''
    smoothPlot(Y,posteriors=clr[,c('Memory Eff','Memory Treg','Naive Eff','Naive Treg','purple','pink','lightblue')],chulls=FALSE,clusters.col=c('black','red','darkgreen','blue','purple','pink','lightblue'), outliers=TRUE,ellipse.lwd=3,ylab=ylab,xlab=xlab)
    xlab <- paste(round(p$coefficients[,,1],2),'*',rownames(p$coefficients),collapse=' + ',sep='')
    ylab <- paste(round(p$coefficients[,,2],2),'*',rownames(p$coefficients),collapse=' + ',sep='')
    print(p$loadings)
    print(p$coefficients[,,1])
    print(p$coefficients[,,2])
    #plot(NULL,xlim=range(p$loadings[,1]),ylim=range(p$loadings[,2]))
    #for (m in MARKERS) {
        #segments(0,0,p$loadings[m,1],p$loadings[m,2])
        #text(p$loadings[m,1],p$loadings[m,2], label=m)
    #}
    #segments(0,0,p$loadings[m,1],p$loadings[m,2])
    #segments(0,0,p$Yloadings[,1],p$Yloadings[,2])
    #text(p$Yloadings[,1],p$Yloadings[,2],label=m)
    title(paste(nextElem(figure.labels),DOSES[[i]], sep='\t'), adj=0)
}
cell.types <- c('Memory Eff','Memory Treg','Naive Eff','Naive Treg')
ylim <- range(sapply(c(cell.types,'purple','pink','lightblue'), function(cell.type) colMedians(fcs.data[as.logical(clr[,cell.type]),paste('diff.PSTAT5',1:4,sep='.')]) ))
plot(NULL, xlim=c(0,3), ylim=ylim, xaxt='n', xlab='dose', ylab='pSTAT5 MFI')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
cols <- c('black','red','darkgreen','blue','purple','pink','lightblue')
for (cell.type in c(cell.types,'purple','pink','lightblue')) {
    mfi <- fcs.data[which(as.logical(clr[,cell.type])),paste('diff.PSTAT5',1:4,sep='.')]
    lines(0:3, colMedians( mfi ), col=cols[[i]], lwd=3)
    #polygon(
            #c(0:3,3:0),
            #c( colQuantile( mfi, prob=.75 ),
            #rev(colQuantile( mfi, prob=.25)) ),
            #col=do.call('rgb', c(as.list(t(col2rgb(i)/255)),alpha=.1)),
            #border=i,
            #lwd=.5)
    i <- i+1
} 
#a <- xx[,paste('diff.PSTAT5',1:4,sep='.')]
#lines(0:3, colMedians( a ), col='purple', lwd=2)
#polygon(
        #c(0:3,3:0),
        #c( colQuantile( a, prob=.75 ),
        #rev(colQuantile( a, prob=.25)) ),
        #col=do.call('rgb', c(as.list(t(col2rgb('purple')/255)),alpha=.1)),
        #border='purple',
        #lwd=.5)
dev.off()

#multivariate view
plotClusters(fcs.data[,MARKERS],posteriors=clr[,c('Memory Eff','Memory Treg','Naive Eff','Naive Treg','pink','lightblue')],chulls=FALSE,clusters.col=c('black','red','darkgreen','blue','pink','lightblue'), outliers=TRUE,ellipse.lwd=3)

#pink and purple do not respond so do not include them
pdf('~/Thesis/figures/plsr-lymphocytes-clusters.pdf',width=10,height=10)
#univariate view
par(mfrow=c(3,3))
for (marker in c(MARKERS,'diff.PSTAT5.3','diff.PSTAT5.4'))
smoothPlot1D(fcs.data[,marker],posteriors=clr[,c('Memory Eff','Memory Treg','Naive Eff','Naive Treg','lightblue')],chulls=FALSE,clusters.col=c('black','red','darkgreen','blue','lightblue'), outliers=TRUE,ellipse.lwd=3, main=marker)
dev.off()


#### Everything else
print(f <- file.path(BASE.DIR,'RData','pstat5-join', sprintf('%s_%s.RData',individual,day)))
print(load(f))
fcs.data <- baseline.relative.pstat5(fcs.data)
fcs.data <- transform.scatter(applyTransforms(fcs.data,transforms))
print(load(file.path(BASE.DIR,'CLR','CB00086S_0U_2012-09-18.RData')))
MARKERS <- c('FSCA','SSCA','CD25','CD3','CD4','CD45RA','CD56','CD8','FOXP3')

par(mfrow=c(1,1))
fm <- as.formula(paste(sprintf('diff.PSTAT5.%d',2), paste(MARKERS, collapse='+'), sep='~'))
plsr(fm, data=as.data.frame(fcs.data),scale=FALSE)->p
Y <- p$scores[,1:2]
smoothPlot(Y,outliers=T)

#g<-locator(type='l')
#dput(g)
g <- structure(list(x = c(-3.42781386640277, -2.83331075199373, -2.23880763758468, -2.29825794902559, -3.04138684203689, -4.02231698081181, -4.73572071810266, -4.94379680814583, -4.46819431661859, -3.54671448928458, -3.2197377763596), y = c(0.876407526174889, 0.876407526174889, 0.614897660608943, 0.291856061968656, 0.0303461964027103, 0.153409662551391, 0.322621928505827, 0.584131794071773, 0.799492859831964, 0.922556325980644, 0.861024592906304)), .Names = c("x", "y"))
lines(g$x,g$y,col='purple')
j <- as.logical(point.in.polygon(Y[,1],Y[,2],g$x,g$y))

#g<-locator(type='l')
#dput(g)
g <- structure(list(x = c(9.59180433915526, 8.46224842177809, 7.68939437304633, 7.89747046308949, 8.61087420038035, 9.35400309339165, 9.79988042919843, 10.3349332321666, 10.097131986403, 9.59180433915526), y = c(-0.000419670134459622, 0.445685394654507, 0.799492859831964, 1.07638565866649, 0.953322192517815, 0.707195260220453, 0.430302461385922, 0.153409662551391, -0.0773343364773852, -0.0158026034030448)), .Names = c("x", "y"))
lines(g$x,g$y,col='pink')
k <- as.logical(point.in.polygon(Y[,1],Y[,2],g$x,g$y))

#g<-locator(type='l')
#dput(g)
g <- structure(list(x = c(2.99281976921489, 4.62770333383976, 5.57890831689423, 6.35176236562598, 7.18406672579864, 7.30296734868045, 6.64901392283051, 5.81670956265785, 5.25193160396925, 4.5385278666784, 3.854849285108, 3.08199523637625, 2.33886634336495, 1.684912917515, 1.0012343359446, 1.0012343359446, 1.80381354039681, 2.30914118764449, 2.96309461349444), y = c(0.799492859831964, 0.99947099232357, 1.01485392559215, 0.968705125786399, 0.830258726369134, 0.368770728311582, 0.0611120629398807, -0.185014869357481, -0.384993001849086, -0.538822334534937, -0.508056467997767, -0.569588201072107, -0.554205267803522, -0.461907668192012, -0.384993001849086, -0.0465684699402148, 0.384153661580167, 0.676429393683283, 0.799492859831964)), .Names = c("x", "y"))
lines(g$x,g$y,col='lightblue')
l <- as.logical(point.in.polygon(Y[,1],Y[,2],g$x,g$y))

#g<-locator(type='l')
#dput(g)
g <- structure(list(x = c(-1.37677812169157, -0.782275007282528, -0.128321581432581, -0.068871269991677, -0.722824695841624, -1.43622843313247, -1.88210576893926, -1.94155608038016, -1.91183092465971, -1.46595358885293, -1.28760265453021), y = c(-0.323461268774746, -0.308078335506161, -0.523439401266352, -1.20028846508409, -1.38488366430712, -1.36950073103853, -1.01569326586107, -0.815715133369468, -0.492673534729182, -0.323461268774746, -0.338844202043331)), .Names = c("x", "y"))
lines(g$x,g$y,col='orange')
m <- as.logical(point.in.polygon(Y[,1],Y[,2],g$x,g$y))

par(mfrow=c(1,1))
fm <- as.formula(paste(sprintf('diff.PSTAT5.%d',4), paste(MARKERS, collapse='+'), sep='~'))
plsr(fm, data=as.data.frame(fcs.data),scale=FALSE)->p
Y <- p$scores[,1:2]
smoothPlot(Y,outliers=T)

#g<-locator(type='l')
g <- structure(list(x = c(3.90086813263438, 3.51311752200631, 3.34985410700502, 3.39066996075534, 3.63556508325728, 4.10494740138599, 4.65596142701535, 4.737593134516, 4.61514557326503, 4.43147423138858, 4.2069870357618, 3.92127605950954, 3.8396443520089, 3.6763809370076, 3.6763809370076), y = c(6.32942929401107, 5.61260526180101, 4.6316881650925, 3.80168139095453, 3.57531590709872, 4.06577445545298, 4.93350881023358, 5.72578800372891, 6.36715687465371, 6.63124993915215, 6.59352235850952, 6.32942929401107, 6.17851897144054, 5.87669832629946, 5.87669832629946)), .Names = c("x", "y"))
lines(g$x,g$y,col='yellow')
n <- as.logical(point.in.polygon(Y[,1],Y[,2],g$x,g$y))

#identification of pink and blue high density clusters
clr <- CLR
clr <- cbind(clr,purple=as.numeric(j))
clr <- cbind(clr,pink=as.numeric(k))
clr <- cbind(clr,lightblue=as.numeric(l))
clr <- cbind(clr,orange=as.numeric(m))
clr <- cbind(clr,yellow=as.numeric(n))

pdf('~/Thesis/figures/plsr-nonlymphocytes.pdf',width=10,height=10)
figure.labels <- iter(paste(letters,')',sep=''))
par(mfrow=c(2,2),mai=c(.75,.75,.5,.1))
for (i in 2:4) {
    fm <- as.formula(paste(sprintf('diff.PSTAT5.%d',i), paste(MARKERS, collapse='+'), sep='~'))
    plsr(fm, data=as.data.frame(fcs.data),scale=FALSE)->p
    Y <- p$scores[,1:2]
    smoothPlot(Y,posteriors=clr[,c('Lymphocytes','purple','pink','lightblue','orange','yellow')],clusters.col=c('black','purple','pink','lightblue','orange','yellow'),chulls=FALSE,outliers=TRUE,ellipse.lwd=3)
    #smoothPlot(Y,classification=clr[,'Lymphocytes'],chulls=TRUE,outliers=TRUE,ellipse.lwd=3)
    title(paste(nextElem(figure.labels),DOSES[[i]], sep='\t'), adj=0)
}
#dose response
ylim <- range(sapply(c('Lymphocytes','purple','pink','lightblue', 'orange', 'yellow'), function(cell.type) colMedians(fcs.data[as.logical(clr[,cell.type]),paste('diff.PSTAT5',1:4,sep='.')]) ))
plot(NULL, xlim=c(0,3), ylim=ylim, xaxt='n', xlab='dose', ylab='pSTAT5 MFI')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
cols <- c('black','purple','pink','lightblue','orange','yellow')
for (cell.type in c('Lymphocytes','purple','pink','lightblue','orange','yellow')) {
    mfi <- fcs.data[which(as.logical(clr[,cell.type])),paste('diff.PSTAT5',1:4,sep='.')]
    lines(0:3, colMedians( mfi ), col=cols[[i]], lwd=3)
    i <- i+1
} 
dev.off()


#MARKERS <- c('FSCA','SSCA','CD25','CD3','CD4','CD45RA','CD56','CD8','FOXP3')

#univariate view
#yellow and pink are interesting
pdf('~/Thesis/figures/plsr-nonlymphocytes-clusters.pdf',width=10,height=10)
#univariate view
par(mfrow=c(3,3)) 
for (marker in MARKERS)
#smoothPlot1D(fcs.data[,marker],posteriors=clr[,c('Lymphocytes','purple','pink','lightblue','orange','yellow')],chulls=FALSE,clusters.col=c('black','purple','pink','lightblue','orange','yellow'), outliers=TRUE,ellipse.lwd=3, main=marker)
smoothPlot1D(fcs.data[,marker],posteriors=clr[,c('Lymphocytes','pink','yellow')],chulls=FALSE,clusters.col=c('black','pink','yellow'), outliers=TRUE,ellipse.lwd=3, main=marker)
dev.off()

plotClusters(fcs.data[,c('SSCA','FSCA')],posteriors=clr[,c('Lymphocytes','purple','pink','lightblue','orange')],chulls=FALSE,clusters.col=c('black','purple','pink','lightblue','orange'), outliers=TRUE,ellipse.lwd=3)


