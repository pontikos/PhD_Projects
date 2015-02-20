source('~nikolas/bin/FCS/fcs.R')

beads <- read.csv('~nikolas/Projects/flowBeads/Beads.Stats/beads.fcs2.kmedoids.stats') 
beads$date <- as.Date(beads$date)
#beads.mef <- colMeans(beads[,paste('mfi',2:6,sep='.')])
beads.mef <- dakomef[,'APC'][2:6]
beads <- do.call('rbind',by(beads,beads$date,function(x) colMeans(x[,paste('mfi',1:6,sep='.')])))


memory.name <- 'lymphocytes-cd4-cd127hi-127hi-cd45raneg.fcs'
fcsFiles <- unique(unlist(lapply(strsplit(list.files(pattern='cad.*.fcs', path='~nikolas/dunwich/Projects/IL2RA/FCS.Gated/Calli/all.gated'), '-'), function(x) x[[1]])))
d <- data.frame()
for (fcsFile in fcsFiles) { 
  memory<-read.FCS(memory.fcsFile<-paste(file.path('~nikolas/dunwich/Projects/IL2RA/FCS.Gated/Calli/all.gated',fcsFile),memory.name, sep='-'), channel=c('CD25'),TRANS=log10)
  print(memory.fcsFile)
  memory.date <- getDate(memory.fcsFile)
  beads.mfi <- as.numeric(beads[memory.date,paste('mfi',1:6,sep='.')])
  MEF <- function(x) {
        return(cbind(1,x)%*%coefficients(lm(log10(beads.mef) ~ log10(beads.mfi[2:6])))) 
  }
  memory.cd25.mfi <- mean(memory[,'CD25'])
  memory.cd25.mef <- MEF(mean(memory[,'CD25']))
  #
  d <- rbind(d,data.frame(fcsFile=fcsFile,
                            date=as.Date(memory.date),
                            mfi.1=beads.mfi[1],
                            mef.1=MEF(log10(beads.mfi[1])),
                            mfi.2=beads.mfi[2],
                            mef.2=MEF(log10(beads.mfi[2])),
                            memory.cd25.mfi=memory.cd25.mfi,
                            memory.cd25.mef=memory.cd25.mef
                            ))
}

x <- read.csv('~/Projects/IL2RA/Calli_CD25bright_CBR200.csv')
x$fcsFile <- tolower(x$fcsFile)
d$fcsFile <- as.character(d$fcsFile)
print(dim(x <- merge(x[,c('individual','fcsFile')], d, by='fcsFile')))
x$fcsFile <- gsub('.fcs','',x$fcsFile)

#recalled individuals
r <- read.csv('~/Projects/IL2RA/CellPhenotypes/recalled.individuals.pch') 
#date mismatch always fucks things up on merge
print(dim(r <- merge(r[,c('individual','fcsFile','pch')], x, all.x=TRUE, all.y=FALSE)))
r <- r[order(r$individual,as.Date(r$date)),]
r$pch <- as.character(r$pch)
r$date <- as.Date(r$date)
#get pch for individuals
print(dim(x<-merge(x, r[c(TRUE,FALSE),c('individual','pch')], all.x=TRUE)))
print(table(x$pch <- ifelse(is.na(x$pch),'x',x$pch)))
r <- na.omit(r)


library(scales)
library(ggplot2)


# no mef correction
pdf('~nikolas/GoogleDrive/Phd/Thesis/IL2RA/figures/CD25-MFI-time-effect-repeatability.pdf')
g <- ggplot(d, aes(date, log10(mfi.1))) + scale_x_date(labels=date_format( "%b")) + xlab("") + ylab("MFI")# + theme(panel.background = element_rect(fill='white',col='black'))
g <- g + geom_point(col='blue') + geom_smooth(col='blue',method='loess')
g <- g + geom_point(data=d, aes(date, log10(mfi.2)), col='blue')+geom_smooth(data=d, aes(date, log10(mfi.2)), col='blue',method='loess')
g <- g+geom_point(data=d, aes(date, memory.cd25.mfi), col='black')+geom_smooth(data=d, aes(date, (memory.cd25.mfi)), col='black', method='loess')
g <- g + geom_point(data=r, aes(x=date, y=(memory.cd25.mfi), group=individual, col=individual))+geom_line(data=r, aes(x=date, y=(memory.cd25.mfi), group=individual, col=individual),size=.25)+guides(colour=FALSE)
print(g)
dev.off()
#
pdf('~nikolas/GoogleDrive/Phd/Thesis/IL2RA/figures/CD25-MFI-repeatability.pdf')
rsquared <- round(cor((r$memory.cd25.mfi[c(TRUE,FALSE)]),y=(r$memory.cd25.mfi[c(FALSE,TRUE)]))**2,digits=3)
ggplot(r, aes(x=(memory.cd25.mfi[c(TRUE,FALSE)]),y=(memory.cd25.mfi[c(FALSE,TRUE)]),col=individual[c(TRUE,FALSE)])) + geom_point() + geom_abline(b=1,a=0) + xlab('MFI Day 1') + ylab('MFI Day2') + guides(colour=FALSE) + xlim(range((r$memory.cd25.mfi))) + ylim(range((r$memory.cd25.mfi))) + ggtitle(bquote(r^2 == .(rsquared)))
dev.off() 

# mef correction
pdf('~nikolas/GoogleDrive/Phd/Thesis/IL2RA/figures/CD25-MFI-time-effect-beads-normalised.pdf')
# mef 
g <- ggplot(d, aes(date, (mef.1))) + scale_x_date(labels=date_format( "%b")) + xlab("") + ylab("MEF")# + theme(panel.background = element_rect(fill='white',col='black'))
g <- g + geom_point(col='blue') + geom_smooth(col='blue',method='loess')
g <- g + geom_point(data=d, aes(date, (mef.2)), col='blue')+geom_smooth(data=d, aes(date, (mef.2)), col='blue',method='loess')
g <- g+geom_point(data=d, aes(date, (memory.cd25.mef)), col='black')+geom_smooth(data=d, aes(date, (memory.cd25.mef)), col='black', method='loess')
#g <- g+geom_hline(yintercept=log10(y)[1:2],col='blue',linetype='dashed') 
g <- g + geom_point(data=r, aes(x=date, y=(memory.cd25.mef), group=individual, col=individual))+geom_line(data=r, aes(x=date, y=(memory.cd25.mef), group=individual, col=individual),size=.25)+guides(colour=FALSE)
print(g)
dev.off()
#
pdf('~nikolas/GoogleDrive/Phd/Thesis/IL2RA/figures/CD25-MFI-beads-normalised.pdf')
rsquared <- round(cor((r$memory.cd25.mef[c(TRUE,FALSE)]),y=(r$memory.cd25.mef[c(FALSE,TRUE)]))**2,digits=3)
ggplot(r, aes(x=(memory.cd25.mef[c(TRUE,FALSE)]),y=(memory.cd25.mef[c(FALSE,TRUE)]),col=individual[c(TRUE,FALSE)])) + geom_point() + geom_abline(b=1,a=0) + xlab('MEF Day 1') + ylab('MEF Day2') + guides(colour=FALSE) + xlim(range((r$memory.cd25.mef))) + ylim(range((r$memory.cd25.mef))) + ggtitle(bquote(r^2 == .(rsquared)))
dev.off()


