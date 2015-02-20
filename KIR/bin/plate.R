library(plyr)
library(reshape)
library(foreign)
library(cluster) 
library(lattice)
library(latticeExtra)
library(mitools)

qc.experiment <- function(d,ref.max.diff.threshold=2) {
    #drop samples which have less than 4 background wells
    drop.samples <- as.character(subset(d, ref.nwells!=4, 'sample'))
    cat('Dropping',length(subset(d, ref.nwells!=4)[,1]),'samples with less than 4 background wells:', drop.samples, '\n')
    d <- subset(d, ref.nwells==4)
    #drop samples which have background well diff greater than threshold
    cat('Dropping',length(subset(d, ref.max.diff >= ref.max.diff.threshold)[,1]),'samples with background well diff greater than', ref.max.diff.threshold, '\n')
    d <- subset(d, ref.max.diff < ref.max.diff.threshold)
    #drop samples which have background median outside of [22; 30]
    cat('Dropping',length(subset(d, !(22 < ref.median & ref.median < 30))[,1]),'samples with background well median outside of [22; 30] \n')
    d <- subset(d, (22 < ref.median & ref.median < 30))
    return(d)
}


qc.experiment <- function(d, ref.max.diff.threshold=2) {
    #drop samples which have less than 4 background wells
    drop.samples <- as.character(subset(d, ref.nwells!=4, 'sample'))
    cat('Dropping',length(subset(d, ref.nwells!=4)[,1]),'samples with less than 4 background wells:', drop.samples, '\n')
    d <- subset(d, ref.nwells==4)
    #drop samples from plate 22
    cat('Dropping',length(subset(d, box384==22)[,1]),'samples from plate 22\n')
    d <- subset(d, box384 != 22)
    #drop samples which have background well diff greater than threshold
    cat('Dropping',length(subset(d, ref.max.diff >= ref.max.diff.threshold)[,1]),'samples with background well diff greater than', ref.max.diff.threshold, '\n')
    d <- subset(d, ref.max.diff < ref.max.diff.threshold)
    #drop samples which have background median outside of [22; 30]
    cat('Dropping',length(subset(d, !(22 < ref.median & ref.median < 30))[,1]),'samples with background well median outside of [22; 30] \n')
    d <- subset(d, (22 < ref.median & ref.median < 30))
    return(d)
}




### using melt/cast
load.experiment <- function(experiment.file, QC=NULL, trans=function(x)x) {
    options(stringsAsFactors=FALSE)
    d <- read.csv(experiment.file)
    #nwells
    ref <- data.matrix(d[,paste('CrossingPoint_ref_well384',1:4,sep='.')])
    d$ref.nwells <- rowSums(!is.na(ref))
    DL1 <- data.matrix(d[,paste('CrossingPoint_DL1_well384',1:4,sep='.')])
    d$DL1.nwells <- rowSums(!is.na(DL1))
    DS1 <- data.matrix(d[,paste('CrossingPoint_DS1_well384',1:4,sep='.')])
    d$DS1.nwells <- rowSums(!is.na(DS1))
    #well median
    d$ref.median <- apply(ref, 1, median, na.rm=T)
    #
    d$DL1.median <- apply(trans(DL1-ref), 1, median, na.rm=T)
    d$DL1.ct.median <- apply(DL1, 1, median, na.rm=T)
    #
    d$DS1.median <- apply(trans(DS1-ref), 1, median, na.rm=T)
    d$DS1.ct.median <- apply(DS1, 1, median, na.rm=T)
    #well var
    d$ref.var <- apply(ref, 1, var, na.rm=T)
    d$DL1.var <- apply(DL1-ref, 1, var, na.rm=T)
    d$DS1.var <- apply(DS1-ref, 1, var, na.rm=T)
    #well max diff
    d$ref.max.diff <- apply(ref, 1, function(x) ifelse(all(is.na(sort(x))),NA,diff(range(sort(x)))) )
    d$ref.ct.range <- apply(ref, 1, function(x) ifelse(all(is.na(sort(x))),NA,diff(range(sort(x)))) )
    #
    d$DL1.max.diff <- apply(DL1-ref, 1, function(x) ifelse(all(is.na(sort(x))),NA,diff(range(sort(x)))) )
    d$DL1.ct.range <- apply(DL1, 1, function(x) ifelse(all(is.na(sort(x))),NA,diff(range(sort(x)))) )
    #
    d$DS1.max.diff <- apply(DS1-ref, 1, function(x) ifelse(all(is.na(sort(x))),NA,diff(range(sort(x)))) )
    d$DS1.ct.range <- apply(DS1, 1, function(x) ifelse(all(is.na(sort(x))),NA,diff(range(sort(x)))) )
    #well mean
    d$ref.mean <- rowMeans(ref, na.rm=T)
    d$DL1.mean <- rowMeans(DL1-ref, na.rm=T)
    d$DS1.mean <- rowMeans(DS1-ref, na.rm=T)
    #t1d
    d$t1d <- d$t1d-1
    cat(sum(table(d[!duplicated(d$sample),'t1d'])), '=', table(d[!duplicated(d$sample),'t1d']), '\n')
    #QC
    if (!is.null(QC)) d <- QC(d)
    cat(sum(table(d[!duplicated(d$sample),'t1d'])), '=', table(d[!duplicated(d$sample),'t1d']), '\n')
    #normalise peak align
    d <- normalise.peak.align(d)
    #the leftovers plate contains no calibrators
    d <- normalise.well.align(d)
    #if all wells are NA and background is not NA then copy number state is 0
    d[is.na(d$DL1.median), 'DL1.cn.state'] <- 0
    #if there is only 1 or less foreground well then copy number state is 0
    d[which( d[,'DL1.nwells'] <= 1 ), 'DL1.cn.state']  <- 0
    #if there is only 1 or less foreground well then copy number state is 0
    d[is.na(d$DS1.median), 'DS1.cn.state'] <- 0
    #if there is only 1 or less foreground well then copy number state is 0
    d[which( d[,'DS1.nwells'] <= 1 ), 'DS1.cn.state']  <- 0
    #ca:co
    table(d$t1d)
    pch <- seq(1,length(unique(d$box384))) 
    sym <- data.frame(Rows(trellis.par.get('superpose.symbol'), pch))
    sym$pch <- pch
    sym$box384 <- unique(d$box384)
    d <- merge(d, sym)
    cat(sum(table(d[!duplicated(d$sample),'t1d'])), '=', table(d[!duplicated(d$sample),'t1d']), '\n')
    return(d)
}

### transform which maps quantiles to position 1 to 5
normalise.quantile <- function(d) {
    quantile.fun <- function(x) {
            co <- coefficients(lm(1:5 ~ quantile(x, na.rm=T)))
            return(co[[1]]+co[[2]]*x)
    }
    ddply(d, .(box384), function(x) {
                   x$DL1.quantilenorm.median <- quantile.fun(x$DL1.median)
                   x$DS1.quantilenorm.median <- quantile.fun(x$DS1.median)
                   return(x)
            })
} 

 
### transform which maps peaks (these are assumed to correspond to cn 1 and 2) to position 1 and 2
normalise.peak.align <- function(d) {
    peak.align.fun <- function(x) {
            x.sort <- sort(x)
            #identify peaks
            #peaks <- sort(as.vector(pam(x.sort,k=3)$medoids))[2:3]
            #co <- coefficients( lm( c(1, 2) ~ peaks ) )
            #return(list(x=co[[1]]+co[[2]]*x,peak1=peaks[[1]],peak2=peaks[[2]], peaknorm.b=co[[2]], peaknorm.a=co[[1]]))
            peaks <- sort(as.vector(pam(x.sort[-2<x.sort & x.sort<2],k=2)$medoids))
            co <- coefficients( lm( c(1, 2) ~ peaks ) )
            return(list(x=co[[1]]+co[[2]]*x,peak1=peaks[[1]],peak2=peaks[[2]], peaknorm.b=co[[2]], peaknorm.a=co[[1]]))
    }
    ddply(d, .(box384), function(x) {
                   p <- peak.align.fun(x$DL1.median)
                   x$DL1.peaknorm.median <- p$x
                   x$DL1.peaknorm.1 <- p$peak1
                   x$DL1.peaknorm.2 <- p$peak2
                   x$DL1.peaknorm.b <- p$peaknorm.b
                   x$DL1.peaknorm.a <- p$peaknorm.a
                   p <- peak.align.fun(x$DS1.median)
                   x$DS1.peaknorm.median <- p$x
                   x$DS1.peaknorm.1 <- p$peak1
                   x$DS1.peaknorm.2 <- p$peak2
                   x$DS1.peaknorm.b <- p$peaknorm.b
                   x$DS1.peaknorm.a <- p$peaknorm.a
                   return(x)
            })
} 

### transform which maps wells with known copy numbers to position 1 and 2
normalise.well.align <- function(d) {
    well.align.fun <- function(x, gene) {
            calibrators <- x[which((x$sample_type == 'calibrator')),] 
            cn <- calibrators[,sprintf('%s.cn.state',gene)]
            medians <- calibrators[,sprintf('%s.median',gene)]
            co <- coefficients( lm(cn ~ medians) )
            return(co[[1]]+co[[2]]*x[,sprintf('%s.median',gene)])
    }
    ddply(d, .(box384), function(x) {
                   #for leftovers do no normalisation because we have no calibrators on that plate
                   if (unique(x$box384) == 'leftovers') {
                    x$DL1.wellnorm.median <- x$DL1.median
                    x$DS1.wellnorm.median <- x$DS1.median
                   } else {
                    x$DL1.wellnorm.median <- well.align.fun(x,'DL1')
                    x$DS1.wellnorm.median <- well.align.fun(x,'DS1')
                   }
                   return(x)
            })
} 

### 
# Assign to zero
zero.cn.trans <- function(d, sd.jitter=.02) {
    d <- subset(d, sample_type %in% c('calibrator', 'sample'))
        for (x in c('DS1', 'DL1')) {
            nwells <- sprintf('%s.nwells',x)
            cn.state <- sprintf('%s.cn.state',x)
            for (normtrans in c('peaknorm.', 'wellnorm.')) {
                field <- sprintf('%s.%smedian',x,normtrans)
                # if cn state is zero then set to zero
                d[which(d[,cn.state] == 0), field] <- 0
                # anything which is zero or less is centered around zero with some jitter
                i <- which(d[,field] <=0)
                d[i, field] <- rnorm(length(i),0,sd.jitter)
                # if cn state is zero then set to max Ct
                d[which(d[,cn.state] == 0), paste(x,'median',sep='.')] <- NA
            }
    }
    return(d)
}

###
process.params <- function(params, formal.params) {
    params <- c(params, formal.params)
    params <- params[names(params)]
    return(params)
}

###
panel.extra.peaknorm <- function(d, trans, field, subscripts) {
                   f <- unlist(strsplit(field, '\\.'))[[1]]
                   p1 <- d[subscripts, paste(f,'peaknorm','1',sep='.')]
                   p2 <- d[subscripts, paste(f,'peaknorm','2',sep='.')]
                   panel.abline(h=trans(p1), lwd=.5, lty=1)
                   panel.abline(h=trans(p2), lwd=.5, lty=1)
}


###
panel.extra.abline <- function(d, trans, field, subscripts, h=1:2) {
                   f <- unlist(strsplit(field, '\\.'))[[1]]
                   p1 <- d[subscripts, paste(f,'peaknorm','1',sep='.')]
                   p2 <- d[subscripts, paste(f,'peaknorm','2',sep='.')]
                   a <- d[subscripts, paste(f,'peaknorm','a',sep='.')]
                   b <- d[subscripts, paste(f,'peaknorm','b',sep='.')]
                   panel.abline(h=trans(p1*b+a), lwd=.5, lty=1)
                   panel.abline(h=trans(p2*b+a), lwd=.5, lty=1)
}


###
plot.experiment <- function(
                                    field,
                                    experiment,
                                    plot.fun=stripplot,
                                    relation="same",
                                    trans=identity,
                                    panel.extra=NULL,
                                    ...) {
    d <- experiment
    x <- trans(experiment[,field])
    params <- process.params(formal.params=as.list(formals()), params=as.list(sys.calls()[[1]])[-1])
    panel.call <- match.fun(sprintf('panel.%s', as.character(params$plot.fun)))
    panel.function <- function(x,y,subscripts, cex, pch, ...) {
               print(subscripts)
               #panel.abline(h=1, lwd=.5, lty=2)
               #panel.abline(h=2, lwd=.5, lty=2)
               #
               if (!is.null(panel.extra)) panel.extra(d, trans, field, subscripts)
               #panel.call(...)
               #pch <- d[subscripts, 'pch']
               fill <- c('blue', 'red')[d[subscripts, 't1d']+1]
               col <- c('blue', 'red')[d[subscripts, 't1d']+1]
               #alpha <- d[subscripts, 'alpha']
               #panel.violin(...)
               #panel.call(fill=fill, ...)
               #browser()
               print(table(d[subscripts,c('t1d','box384')]))
               #panel.stripplot(fill=d[subscripts,'fill'], col=d[subscripts,'col'], subscripts=subscripts, ...)
               set.seed(1234)
               if (length(calib <- which(d[subscripts,'sample_type']=='calibrator'))>1)
               {
                   panel.points(jitter(as.numeric(x[-calib])),y[-calib],col=col[-calib],fill=fill[-calib],pch=pch,cex=cex,alpha=.5)
                   panel.points(as.numeric(x[calib]),y[calib],col='black',fill='black',pch=pch,cex=cex)
               } else {
                   panel.points(jitter(as.numeric(x)),y,col=col,fill=fill,pch=pch,cex=cex,alpha=.5)
                   #panel.stripplot(x,y,subscripts=subscripts, ...)
               }
              }
    par(mgp=c(1.5, 1, 0))
    print(plot.fun(
              x ~ factor(t1d) | factor(box384),
              data=experiment,
              #layout=c(1,18),
              layout=c(18,1),
              #strip=F, strip.left=T,
              #xlab=paste(params[names(params)][-1], sep=" ", collapse=" "),
              #auto.key=list(columns=2),
              groups=as.vector(experiment$t1d),
              #scales=list(y=list(relation=relation, draw=F)),
              panel=panel.function,
              #col=c("black", "red"),
              pch=21,
              cex=.5,
              jitter.data=TRUE,
              ...
              ))
}
###
well.mapping.96.to.384 <- function( ) {
    g <- as.data.frame(expand.grid(x=toupper(letters[1:8]), y=1:12))
    g.96 <- paste(g$x, sprintf("%.2d",g$y), sep="")
    g.384 <- t( apply( expand.grid(1:8, 1:12) , 1, function(x) {
                      d <- as.data.frame(expand.grid(x=toupper(letters[(2*x[1]-1):(2*x[1])]), y=(2*x[2]-1):(2*x[2])))
                      sort(paste(d$x, sprintf('%.1d',d$y), sep=''))
                  }) )
    g<-cbind(g.96, g.384)
    colnames(g) <- c('well96', paste('well384', 1:4, sep='.'))
    return(as.data.frame(g, stringsAsFactors=F))
}
###
experiment.to.dta <- function(
                              file.name='../lc_repeats/lc_repeats_2011-11-21.tab',
                              pheno.file.name='../lc_repeats/phenotype_repeats.tab',
                              out.file.name='../lc_repeats/lc_repeats.dta'
                              ) { 
    d <- read.csv(file.name, sep='\t')

    g <- well.mapping.96.to.384()

    data <- ddply(d, c('Experiment.Name'), function(d) {
                    experimentname <- as.character(d[1,'Experiment.Name'])
                    Experiment.Name <- gsub('CAC[O0]_?(.*)_.*', '\\1', experimentname, ignore.case=FALSE)
                    ldply(apply(g, 1, function(x) {
                        d2 <- data.frame(experimentname=experimentname, Experiment.Name=Experiment.Name, plate_well=x[1], stringsAsFactors=FALSE)
                        d <- as.data.frame(t( d )[,d$Position %in% x[-1]], stringsAsFactors=FALSE)
                        d2[,paste('ct_ref',1:4,sep='')] <- as.numeric(d['CrossingPoint', d['Analysis.Name',] == 'aq2dm_dfo'])
                        d2[,paste('ct_tar1',1:4,sep='')] <- as.numeric(d['CrossingPoint', d['Analysis.Name',] == 'aq2dm_cy5'])
                        d2[,paste('ct_tar2',1:4,sep='')] <- as.numeric(d['CrossingPoint', d['Analysis.Name',] == 'aq2dm_fam'])
                        d2[,paste('lc_position',1:4,sep='_')] <-  x[-1]
                        return(d2)
                        }))
            })
    pheno <- read.csv(pheno.file.name, sep='\t')
    pheno$Experiment.Name <- gsub('CAC[O0]_?(.*)', '\\1', pheno$box, ignore.case=FALSE)
    data <- merge(data, pheno, by=c('Experiment.Name', 'plate_well'))
    write.dta(data, file=out.file.name)
}
###
###
load.dta.experiment <- function(
                               file.name='../lc_scale_up/lc_scale_up.dta',
                               replace.na=T,
                               #background.correct.fun=background.substract,
                               #normalise.plate.fun=normalise.peak.align,
                               #filter.expression = expression(sample_type %in% c('calibrator', 'sample'))
                               filter.expression = T
                               ) {
    d<-read.dta(file.name)
    #t1d: 0/1
    d$t1d <- d$t1d-1
    #standardise plate names
    d$experimentname <- gsub('CAC[O0]_?(.*)_.*', '\\1', d$experimentname, ignore.case=FALSE)
    #light.cycler$Experiment.Name <- as.numeric(gsub('leftovers', '0', light.cycler$Experiment.Name))
    #exclude plates
    d$rownames <- as.numeric(rownames(d))
    colnames(d) <- gsub('tar1', 'DL1', colnames(d))
    colnames(d) <- gsub('tar2', 'DS1', colnames(d))
    #filter out some data
    d <- subset(d, eval(filter.expression))
    d <- ddply(d, .(rownames),
               function(x) {
                  ref <- as.numeric(x[paste('ct_ref',1:4, sep='')])
                  ok.wells <- which(!is.na(ref))
                  ref <- ref[ok.wells]
                  x$ref.nwells <- length(ref)
                  #if all background is NA then ignore this sample
                  if ( x$ref.nwells == 0 ) return(NULL)
                  DL1 <- as.numeric(x[,paste('ct_DL1',1:4, sep='')])[ok.wells]
                  DS1 <- as.numeric(x[,paste('ct_DS1',1:4, sep='')])[ok.wells]
                  #number of non NA wells
                  x$DL1.nwells <- length(na.omit(DL1))
                  x$DS1.nwells <- length(na.omit(DS1))
                  #well median
                  x$ref.median <- median(ref, na.rm=T)
                  x$DL1.median <- median(DL1-ref, na.rm=T)
                  x$DS1.median <- median(DS1-ref, na.rm=T)
                  #well mean
                  x$ref.mean <- mean(ref, na.rm=T)
                  x$DL1.mean <- mean(DL1-ref, na.rm=T)
                  x$DS1.mean <- mean(DS1-ref, na.rm=T)
                  #if all wells are NA and background is not NA then copy number state is 0
                  if (is.na(x$DL1.median)) x$DL1.cn.state <- 0
                  #if all wells are NA and background is not NA then copy number state is 0
                  if (is.na(x$DS1.median)) x$DS1.cn.state <- 0
                  return(x)
            })
    #with(calibrators, sample_type == 'calibrator' & plate_well c('H08', 'H09', 'H10', 'H11')
    d[with(d, sample_type=='calibrator' & plate_well %in% c('H08', 'H09', 'H10', 'H11')), 'DL1.cn.state'] <- c(2,2,1,1)
    d[with(d, sample_type=='calibrator' & plate_well %in% c('H08', 'H09', 'H10', 'H11')), 'DS1.cn.state'] <- c(1,1,2,2)
    #if (is.function(normalise.plate.fun)) d<-normalise.plate.fun(d)
    d <- normalise.peak.align(d)
    d <- normalise.well.align(d)
    ten_dnas_lookup <- read.csv('../lc_scale_up/ten_dnas_lookup.csv', header=T, colClasses=c('character', 'character', 'numeric', 'numeric'))
    d <- merge(d, ten_dnas_lookup, by=c('plate_box', 'src_well'), all.x=T)
    d$DL1.cn.state <- apply(cbind( d$DL1.cn.state.x, d$DL1.cn.state.y ), 1, function(x) ifelse(all(is.na(x)), NA, x[!is.na(x)]))
    d$DS1.cn.state <- apply(cbind( d$DS1.cn.state.x, d$DS1.cn.state.y ), 1, function(x) ifelse(all(is.na(x)), NA, x[!is.na(x)]))
    return(d)
}

