source('~nikolas/bin/FCS/fcs.R')
library(reshape)
library(iterators)
CELL.TYPES <- c('Memory Eff','Memory Treg','Naive Eff','Naive Treg')
DOSES <- c('0U','01U','10U','1000U')

# not joined
base.dir <- '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/RData/'
f <- list.files(base.dir, pattern='.*_0U.*.RData')
individuals <- sapply(strsplit(f,'_'),'[[',1)
dates <- gsub('.RData','',sapply(strsplit(f,'_'),'[[',3))
individuals.dates <- unique(cbind(individuals,dates))
pstat5.mfi <- data.frame()
for (i in 1:nrow(individuals.dates)) {
    individual <- individuals.dates[i,1]
    day <- individuals.dates[i,2]
    cat(individual, day, '\n')
    for (dose in DOSES) {
      print(dose)
      filename <- sprintf(base.dir,'%s_%s_%s.RData',individual,dose,day)
      if (!file.exists(filename)) break
      print((load(filename)))
      print(dim(fcs.data))
      print(load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3','magnetic-manual-gates2/CLR', sprintf('%s.RData',paste(individual, dose, day, sep='_')))))
      for (gate in CELL.TYPES) {
          print(length(g <- which(as.logical(CLR[,gate]))))
          pstat5mfi <- median(fcs.data[g,'PSTAT5'])
          pstat5.mfi <- rbind( pstat5.mfi,
                        data.frame(individual=as.character(individual),
                                   date=day,
                                   cell.type=gate,
                                   dose=dose,
                                   pstat5.mfi=pstat5mfi) )
      }
    }
}
pstat5.mfi <- cast(pstat5.mfi, individual + date + cell.type ~ dose)
pstat5.mfi <- na.omit(pstat5.mfi)
save(pstat5.mfi, file=file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/','pstat5.mfi.RData'))


# joined
fcs.files <- list.files(path='~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/', pattern='.*.RData', full.names=TRUE)
pstat5.mfi <- data.frame()
for (i in  1:length(fcs.files)) {
#for (i in  1:3) {
    print(f <- fcs.files[[i]])
    individual <- unlist(strsplit(gsub('.RData','',basename(f)), '_'))[[1]]
    date <- unlist(strsplit(gsub('.RData','',basename(f)), '_'))[[2]]
    print(load(file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/',basename(f))))
    print(load(file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR',basename(f))))
    for (gate in CELL.TYPES) {
        print(gate)
        d <- fcs.data[which(as.logical(CLR[,gate])),]
        mfi <- sapply(paste('PSTAT5', 1:4, sep='.'), function(chan) {
            x <- d[,chan]
            return( median(x) )
            })
        names(mfi) <- DOSES
        base.mfi <- sapply(paste('PSTAT5', 1:4, sep='.'), function(chan) {
            x <- d[,chan] - d[,'PSTAT5.1']
            return( median(x) )
            })
        names(base.mfi) <- paste('base', DOSES, sep='.')
        pstat5.mfi <- rbind( pstat5.mfi,
                            data.frame(individual=as.character(individual),
                                       date=date,
                                       cell.type=gate,
                                       t(mfi),
                                       t(base.mfi),
                                       auc.pstat5=sum((mfi)),
                                       base.auc.pstat5=sum((mfi))
                            ))
    }
}
save(pstat5.mfi, file=file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/','nn.pstat5.mfi.RData'))



# joined
### peak normalised pSTAT5 MFI per cell type
### ungated normalisation
fcs.files <- list.files(path='~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/', pattern='.*.RData', full.names=TRUE)
pstat5.mfi <- data.frame()
for (i in  1:length(fcs.files)) {
    print(f <- fcs.files[[i]])
    individual <- unlist(strsplit(gsub('.RData','',basename(f)), '_'))[[1]]
    date <- unlist(strsplit(gsub('.RData','',basename(f)), '_'))[[2]]
    print(load(file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All-pstat5-normalised/',basename(f))))
    print(load(file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR',basename(f))))
    for (gate in CELL.TYPES) {
        print(gate)
        d <- pstat5[which(as.logical(CLR[,gate])),]
        mfi <- sapply(paste('PSTAT5', 1:4, sep='.'), function(chan) {
            x <- invlgcl(d[,chan])
            return( median(x) )
            })
        names(mfi) <- DOSES
        #d.norm <- pstat5.normalisation$Lymphocytes.10U$pstat5.norm[which(as.logical(CLR[,gate])),]
        d.norm <- pstat5.normalisation$Lymphocytes$pstat5.norm[which(as.logical(CLR[,gate])),]
        norm.mfi <- sapply(paste('PSTAT5', 1:4, sep='.'), function(chan) {
             x <- invlgcl(d.norm[,chan])
             return( median(x) )
            })
        names(norm.mfi) <- paste('norm',DOSES,sep='.')
        #d.norm <- pstat5.normalisation$Lymphocytes.1000U$pstat5.norm[which(as.logical(CLR[,gate])),]
        d.norm <- pstat5.normalisation$CD4$pstat5.norm[which(as.logical(CLR[,gate])),]
        norm2.mfi <- sapply(paste('PSTAT5', 1:4, sep='.'), function(chan) {
             x <- invlgcl(d.norm[,chan])
             return( median(x) )
            })
        names(norm2.mfi) <- paste('norm2',DOSES,sep='.')
        pstat5.mfi <- rbind( pstat5.mfi,
                            data.frame(individual=as.character(individual),
                                       date=date,
                                       cell.type=gate,
                                       t(mfi),
                                       t(norm.mfi),
                                       t(norm2.mfi),
                                       auc.pstat5=sum((mfi)),
                                       auc.norm.pstat5=sum((norm.mfi)),
                                       auc.norm2.pstat5=sum((norm2.mfi)))
                            )
    }
}
save(pstat5.mfi, file=file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/','nn.peaknorm.pstat5.mfi.RData'))



pdf('~nikolas/Thesis/figures/pstat5-mfi-auc-celltypes.pdf')
par(mfrow=c(2,2))
figure.labels <- iter(paste(letters,')',sep=''))
#a) pSTAT5 MFI
load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/','pstat5.mfi.RData'))
plot(NULL, xlim=c(0,3), ylim=range(pstat5.mfi[,DOSES]), xaxt='n', xlab='dose', ylab='pSTAT5 MFI')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
colMedians <- function(x) apply(x, 2, median)
for (cell.type in CELL.TYPES) {
    mfi <- pstat5.mfi[grep(cell.type,pstat5.mfi$cell.type),DOSES]
    lines(0:3, colMedians( mfi ), col=i, lwd=3)
    polygon(
            c(0:3,3:0),
            c( colQuantile( mfi, prob=.75 ),
            rev(colQuantile( mfi, prob=.25)) ),
            col=do.call('rgb', c(as.list(t(col2rgb(i)/255)),alpha=.1)),
            border=i,
            lwd=.5)
    i <- i+1
}
legend('topleft', CELL.TYPES, text.col=1:4, bty='n')
#b) NN pSTAT5 MFI
load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/','nn.pstat5.mfi.RData'))
plot(NULL, xlim=c(0,3), ylim=range(pstat5.mfi[,paste('X',DOSES,sep='')]), xaxt='n', xlab='dose', ylab='NN pSTAT5 MFI')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
colMedians <- function(x) apply(x, 2, median)
for (cell.type in CELL.TYPES) {
    mfi <- pstat5.mfi[grep(cell.type,pstat5.mfi$cell.type),paste('X',DOSES,sep='')]
    lines(0:3, colMedians( mfi ), col=i, lwd=3)
    polygon(
            c(0:3,3:0),
            c( colQuantile( mfi, prob=.75 ),
            rev(colQuantile( mfi, prob=.25)) ),
            col=do.call('rgb', c(as.list(t(col2rgb(i)/255)),alpha=.1)),
            border=i,
            lwd=.5)
    i <- i+1
}
legend('topleft', CELL.TYPES, text.col=1:4, bty='n')
#c) base pSTAT5 MFI
load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/','pstat5.mfi.RData'))
plot(NULL, xlim=c(0,3), ylim=range(pstat5.mfi[,DOSES]), xaxt='n', xlab='dose', ylab='base pSTAT5 MFI')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
colMedians <- function(x) apply(x, 2, median)
for (cell.type in CELL.TYPES) {
    mfi <- pstat5.mfi[grep(cell.type,pstat5.mfi$cell.type),DOSES] - pstat5.mfi[grep(cell.type,pstat5.mfi$cell.type),'0U']
    lines(0:3, colMedians( mfi ), col=i, lwd=3)
    polygon(
            c(0:3,3:0),
            c( colQuantile( mfi, prob=.75 ),
            rev(colQuantile( mfi, prob=.25)) ),
            col=do.call('rgb', c(as.list(t(col2rgb(i)/255)),alpha=.1)),
            border=i,
            lwd=.5)
    i <- i+1
}
legend('topleft', CELL.TYPES, text.col=1:4, bty='n')
#d) NN base pSTAT5 MFI
load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/','nn.pstat5.mfi.RData'))
plot(NULL, xlim=c(0,3), ylim=range(pstat5.mfi[,paste('base',DOSES,sep='.')]), xaxt='n', xlab='dose', ylab='NN base pSTAT5 MFI')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
colMedians <- function(x) apply(x, 2, median)
for (cell.type in CELL.TYPES) {
    mfi <- pstat5.mfi[grep(cell.type,pstat5.mfi$cell.type),paste('base',DOSES,sep='.')]
    lines(0:3, colMedians( mfi ), col=i, lwd=3)
    polygon(
            c(0:3,3:0),
            c( colQuantile( mfi, prob=.75 ),
            rev(colQuantile( mfi, prob=.25)) ),
            col=do.call('rgb', c(as.list(t(col2rgb(i)/255)),alpha=.1)),
            border=i,
            lwd=.5)
    i <- i+1
}
legend('topleft', CELL.TYPES, text.col=1:4, bty='n')
dev.off()



pdf('~nikolas/Thesis/figures/pstat5-mfi-peaknorm-auc-celltypes.pdf')
#
load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/','nn.peaknorm.pstat5.mfi.RData'))
figure.labels <- iter(paste(letters,')',sep=''))
par(mfrow=c(1,2))
plot(NULL, xlim=c(0,3), ylim=quantile(as.matrix(pstat5.mfi[,paste('norm',DOSES,sep='.')]),probs=c(.1,.95)), xaxt='n', xlab='dose', ylab='Lymph norm pSTAT5 MFI')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
colMedians <- function(x) apply(x, 2, median)
for (cell.type in CELL.TYPES) {
    mfi <- pstat5.mfi[grep(cell.type,pstat5.mfi$cell.type),paste('norm',DOSES,sep='.')]
    lines(0:3, colMedians( mfi ), col=i, lwd=3)
    polygon(
            c(0:3,3:0),
            c( colQuantile( mfi, prob=.75 ),
            rev(colQuantile( mfi, prob=.25)) ),
            col=do.call('rgb', c(as.list(t(col2rgb(i)/255)),alpha=.1)),
            border=i,
            lwd=.5)
    i <- i+1
}
legend('topleft', CELL.TYPES, text.col=1:4, bty='n')
#
plot(NULL, xlim=c(0,3), ylim=quantile(as.matrix(pstat5.mfi[,paste('norm2',DOSES,sep='.')]),probs=c(.1,.95)), xaxt='n', xlab='dose', ylab='CD4+ norm pSTAT5 MFI')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
colMedians <- function(x) apply(x, 2, median)
for (cell.type in CELL.TYPES) {
    mfi <- pstat5.mfi[grep(cell.type,pstat5.mfi$cell.type),paste('norm2',DOSES,sep='.')]
    lines(0:3, colMedians( mfi ), col=i, lwd=3)
    polygon(
            c(0:3,3:0),
            c( colQuantile( mfi, prob=.75 ),
            rev(colQuantile( mfi, prob=.25)) ),
            col=do.call('rgb', c(as.list(t(col2rgb(i)/255)),alpha=.1)),
            border=i,
            lwd=.5)
    i <- i+1
}
legend('topleft', CELL.TYPES, text.col=1:4, bty='n')
dev.off()


# boxplots of pSTAT5 MFI
print(load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/','nn.pstat5.mfi.RData')))
nn.pstat5.mfi <- pstat5.mfi
print(load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/','nn.peaknorm.pstat5.mfi.RData')))
nn.peaknorm.mfi <- pstat5.mfi
print(load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/','pstat5.mfi.RData')))
boxplot(X0U ~ cell.type, data=data.frame(pstat5.mfi))
boxplot(X01U ~ cell.type, data=data.frame(pstat5.mfi))
boxplot(X10U ~ cell.type, data=data.frame(pstat5.mfi))
boxplot(X1000U ~ cell.type, data=data.frame(pstat5.mfi))



