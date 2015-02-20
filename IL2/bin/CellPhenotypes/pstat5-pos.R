source('~nikolas/bin/FCS/fcs.R')
library(reshape)
library(iterators)
CELL.TYPES <- c('Memory Eff','Memory Treg','Naive Eff','Naive Treg')
DOSES <- c('0U','01U','10U','1000U')

# not joined
base.dir <- '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/RData'
f <- list.files(base.dir, pattern='.*_0U.*.RData')
individuals <- sapply(strsplit(f,'_'),'[[',1)
dates <- gsub('.RData','',sapply(strsplit(f,'_'),'[[',3))
individuals.dates <- unique(cbind(individuals,dates))
pstat5.pos <- data.frame()
for (i in 1:nrow(individuals.dates)) {
    individual <- individuals.dates[i,1]
    day <- individuals.dates[i,2]
    print(filename <- file.path(base.dir,sprintf('%s_0U_%s.RData',individual,day)))
    print(file.exists(filename))
    print((load(filename)))
    print(load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3','magnetic-manual-gates2/CLR', sprintf('%s.RData',paste(individual, '0U', day, sep='_')))))
    baseline.pstat5 <- fcs.data[,'PSTAT5']
    baseline.clr <- CLR
    for (dose in DOSES) {
      print(dose)
      filename <- file.path(base.dir,sprintf('%s_%s_%s.RData',individual,dose,day))
      if (!file.exists(filename)) break
      print((load(filename)))
      print(dim(fcs.data))
      print(load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3','magnetic-manual-gates2/CLR', sprintf('%s.RData',paste(individual, dose, day, sep='_')))))
      for (gate in CELL.TYPES) {
          print(length(g <- which(as.logical(CLR[,gate]))))
          pstat5pos <- quantile(baseline.pstat5[which(as.logical(baseline.clr[,gate]))],probs=.99)
          pct.pstat5pos <- 100*length(which(fcs.data[g,'PSTAT5'] >= pstat5pos))/length(g)
          pstat5.pos <- rbind( pstat5.pos,
                        data.frame(individual=as.character(individual),
                                   date=day,
                                   cell.type=gate,
                                   dose=dose,
                                   pstat5.pos=pct.pstat5pos) )
      }
    }
}
#pstat5.pos[,'0U'] <- 1
pstat5.pos <- cast(pstat5.pos, individual + date + cell.type ~ dose)
pstat5.pos <- na.omit(pstat5.pos)
save(pstat5.pos, file='~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/pstat5.pos.RData')


# joined
base.dir <- '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/RData/pstat5-join/'
setwd(base.dir)
filenames <- list.files(base.dir, pattern='.*.RData')
pstat5.pos <- data.frame()
for (filename in filenames) {
    print(filename)
    print(individual <- unlist(strsplit(gsub('.RData','',basename(filename)), '_'))[[1]])
    print(day <- unlist(strsplit(gsub('.RData','',basename(filename)), '_'))[[2]])
    print((load(filename)))
    print(load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3','magnetic-manual-gates2/CLR', basename(filename))))
    baseline.pstat5 <- fcs.data[,'PSTAT5.1']
    for (gate in CELL.TYPES) {
          print(length(g <- which(as.logical(CLR[,gate]))))
          pstat5pos <- quantile(baseline.pstat5[g],probs=.99)
          pct.pstat5pos <- sapply(paste('PSTAT5',1:4,sep='.'), function(pstat5) 100*length(which(fcs.data[g,pstat5] >= pstat5pos))/length(g))
          pstat5.pos <- rbind( pstat5.pos,
                        data.frame(individual=as.character(individual),
                                   date=day,
                                   cell.type=gate,
                                   pstat5.pos=t(pct.pstat5pos)) )
      }
}
colnames(pstat5.pos) <- gsub('pstat5.pos.PSTAT5.1','0U',colnames(pstat5.pos))
colnames(pstat5.pos) <- gsub('pstat5.pos.PSTAT5.2','01U',colnames(pstat5.pos))
colnames(pstat5.pos) <- gsub('pstat5.pos.PSTAT5.3','10U',colnames(pstat5.pos))
colnames(pstat5.pos) <- gsub('pstat5.pos.PSTAT5.4','1000U',colnames(pstat5.pos))
save(pstat5.pos, file='~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/nn.pstat5.pos.RData')

# joined base
base.dir <- '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/RData/pstat5-join'
setwd(base.dir)
filenames <- list.files(base.dir, pattern='.*.RData')
pstat5.pos <- data.frame()
for (filename in filenames) {
    print(filename)
    print(individual <- unlist(strsplit(gsub('.RData','',basename(filename)), '_'))[[1]])
    print(day <- unlist(strsplit(gsub('.RData','',basename(filename)), '_'))[[2]])
    print((load(filename)))
    fcs.data <- baseline.relative.pstat5(fcs.data)
    print(load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3','magnetic-manual-gates2/CLR', basename(filename))))
    for (gate in CELL.TYPES) {
          print(length(g <- which(as.logical(CLR[,gate]))))
          pct.pstat5pos <- sapply(paste('diff.PSTAT5',1:4,sep='.'), function(pstat5) 100*length(which(fcs.data[g,pstat5] > 0))/length(g))
          pstat5.pos <- rbind( pstat5.pos,
                        data.frame(individual=as.character(individual),
                                   date=day,
                                   cell.type=gate,
                                   pstat5.pos=t(pct.pstat5pos)) )
      }
}
colnames(pstat5.pos) <- gsub('pstat5.pos.diff.PSTAT5.1','0U',colnames(pstat5.pos))
colnames(pstat5.pos) <- gsub('pstat5.pos.diff.PSTAT5.2','01U',colnames(pstat5.pos))
colnames(pstat5.pos) <- gsub('pstat5.pos.diff.PSTAT5.3','10U',colnames(pstat5.pos))
colnames(pstat5.pos) <- gsub('pstat5.pos.diff.PSTAT5.4','1000U',colnames(pstat5.pos))
save(pstat5.pos, file='~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/nn.base.pstat5.pos.RData')



zoom <- 5
pdf('~nikolas/Thesis/figures/pstat5-pos-auc-celltypes.pdf', width=2*zoom,height=1*zoom)
par(mfrow=c(1,2))
figure.labels <- iter(paste(letters,')',sep=''))
#a) % pSTAT5+
load('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/pstat5.pos.RData')
plot(NULL, xlim=c(0,3), ylim=range(pstat5.pos[,DOSES]), xaxt='n', xlab='dose', ylab='% pSTAT5+')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
colMedians <- function(x) apply(x, 2, median)
for (cell.type in CELL.TYPES) {
    pct.pos <- pstat5.pos[grep(cell.type,pstat5.pos$cell.type),DOSES]
    lines(0:3, colMedians( pct.pos ), col=i, lwd=3)
    polygon(
            c(0:3,3:0),
            c( colQuantile( pct.pos, prob=.75 ),
            rev(colQuantile( pct.pos, prob=.25)) ),
            col=do.call('rgb', c(as.list(t(col2rgb(i)/255)),alpha=.1)),
            border=i,
            lwd=.5)
    i <- i+1
}
legend('topleft', CELL.TYPES, text.col=1:4, bty='n') 
#b) NN % pSTAT5+
load('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/nn.pstat5.pos.RData')
plot(NULL, xlim=c(0,3), ylim=range(pstat5.pos[,DOSES]), xaxt='n', xlab='dose', ylab='NN % pSTAT5+')
title(nextElem(figure.labels), adj=0)
axis(1, at=0:3, labels=DOSES)
i <- 1
colMedians <- function(x) apply(x, 2, median)
for (cell.type in CELL.TYPES) {
    pct.pos <- pstat5.pos[grep(cell.type,pstat5.pos$cell.type),DOSES]
    lines(0:3, colMedians( pct.pos ), col=i, lwd=3)
    polygon(
            c(0:3,3:0),
            c( colQuantile( pct.pos, prob=.75 ),
            rev(colQuantile( pct.pos, prob=.25)) ),
            col=do.call('rgb', c(as.list(t(col2rgb(i)/255)),alpha=.1)),
            border=i,
            lwd=.5)
    i <- i+1
}
legend('topleft', CELL.TYPES, text.col=1:4, bty='n')
dev.off()


