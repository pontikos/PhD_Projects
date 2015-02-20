### Look at evolution of cell frequencies over different people (with a few recalls over a year).
### Only a handful are recalled.


library(ggplot2)
library(reshape)

setwd('~nikolas/Projects/IL2RA/')

d<-read.csv('Calli_CD25bright_CBR200.csv')

x<-d[,c( 'date', 'Lymphocytes.CD4.CD127hi.CD45RAneg.freqPar', 'Lymphocytes.CD4.CD127hi.CD45RApos.freqPar', 'Lymphocytes.CD4.CD127hi.freqPar', 'Lymphocytes.CD4.CD127hiCD25neg.CD45RAneg.freqPar', 'Lymphocytes.CD4.CD127hiCD25neg.CD45RApos.freqPar', 'Lymphocytes.CD4.CD127intlowCD25hi.CD45RAneg.freqPar', 'Lymphocytes.CD4.CD127intlowCD25hi.CD45RApos.freqPar', 'Lymphocytes.CD4.CD127intlowCD25hi.freqPar', 'Lymphocytes.CD4.CD127intlowCD25hi.HLADRpos.freqPar', 'Lymphocytes.CD4.CD127intlowCD25pos.CD45RAneg.freqPar', 'Lymphocytes.CD4.CD127intlowCD25pos.CD45RApos.freqPar', 'Lymphocytes.CD4.CD127intlowCD25pos.freqPar', 'Lymphocytes.CD4.CD127intlowCD25pos.HLADRneg.freqPar', 'Lymphocytes.CD4.CD127intlowCD25pos.HLADRpos.freqPar', 'Lymphocytes.CD4.CD127negCD25neg.freqPar', 'Lymphocytes.CD4.CD45RAneg.freqPar', 'Lymphocytes.CD4.CD45RApos.freqPar', 'Lymphocytes.CD4.freqPar')]

head(x<-melt(x, 'date'))
x <- x[x$value < 100,]
x <- na.omit(x)
x$variable <- gsub('Lymphocytes.(.*).freqPar','\\1',x$variable)
ggplot(x,aes(x=as.Date(date),y=as.numeric(value),colour=variable)) + geom_point() + geom_smooth(method='loess') + theme(legend.position="none") + ylim(0,100) + facet_wrap(~variable) + ylab('freqpar') + xlab('')


ct <- d[,c('date', paste('Lymphocytes',unique(x$variable),'ct',sep='.'))]


# MeanPE_norm
mean.pe <- d[,c('date',unlist(sapply(paste('Lymphocytes',unique(x$variable),'MeanPE_norm',sep='.'), function(x) grep(x, colnames(d),value=T))))]
mean.pe <- melt(mean.pe, 'date')
mean.pe <- na.omit(mean.pe)
mean.pe$variable <- gsub('Lymphocytes.(.*)_norm','\\1',mean.pe$variable)
ggplot(mean.pe,aes(x=as.Date(date),y=as.numeric(value),colour=variable)) + geom_point() + geom_smooth(method='loess') + theme(legend.position="none") + facet_wrap(~variable) + ylab('MeanPE_norm') + xlab('')

# MeanAlexa488_norm
mean.alexa <- d[,c('date',unlist(sapply(paste('Lymphocytes',unique(x$variable),'MeanAlexa488_norm',sep='.'), function(x) grep(x, colnames(d),value=T))))]
mean.alexa <- melt(mean.alexa, 'date')
mean.alexa <- na.omit(mean.alexa)
mean.alexa$variable <- gsub('Lymphocytes.(.*)_norm','\\1',mean.alexa$variable)
ggplot(mean.alexa,aes(x=as.Date(date),y=as.numeric(value),colour=variable)) + geom_point() + geom_smooth(method='loess') + theme(legend.position="none") + facet_wrap(~variable) + ylab('MeanAlexa488_norm') + xlab('')


# MeanAPC_norm
mean.apc <- d[,c('date',unlist(sapply(paste('Lymphocytes',unique(x$variable),'MeanAPC_norm',sep='.'), function(x) grep(x, colnames(d),value=T))))]
mean.apc <- melt(mean.apc, 'date')
mean.apc <- na.omit(mean.apc)
mean.apc$variable <- gsub('Lymphocytes.(.*).MeanAPC_norm','\\1',mean.apc$variable)
ggplot(mean.apc,aes(x=as.Date(date),y=as.numeric(value),colour=variable)) + geom_point() + geom_smooth(method='loess') + theme(legend.position="none") + facet_wrap(~variable) + ylab('MeanAPC_norm') + xlab('')


# MeanAPC
mean.apc <- d[,c('date',unlist(sapply(paste('Lymphocytes',unique(x$variable),'MeanAPC$',sep='.'), function(x) grep(x, colnames(d),value=T))))]
mean.apc <- melt(mean.apc, 'date')
mean.apc <- na.omit(mean.apc)
mean.apc$variable <- gsub('Lymphocytes.(.*).MeanAPC','\\1',mean.apc$variable)
ggplot(mean.apc,aes(x=as.Date(date),y=as.numeric(value),colour=variable)) + geom_point() + geom_smooth(method='loess') + theme(legend.position="none") + facet_wrap(~variable) + ylab('MeanAPC') + xlab('')


