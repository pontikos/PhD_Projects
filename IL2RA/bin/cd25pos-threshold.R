library(flowBeads)
source('~nikolas/bin/FCS/fcs.R')

PHENO <- '% CD25+ naive'
par(cex.lab=1.5, cex.main=2)

naive.name <- 'lymphocytes-cd4-cd127hi-127hi-cd45rapos.fcs'
naive.cd25pos.name <- 'lymphocytes-cd4-cd127hicd25pos-cd45rapos.fcs'

b <- read.csv('~nikolas/Projects/flowBeads/Beads.Stats/beads.fcs2.kmedoids.stats')

blank.quantiles <- data.frame()
for (bf in b$file.name) {
    bf <- BeadFlowFrame(file.path('~/dunwich/Projects/IL2RA/FCS.Beads/FCS2.0/FCS',sprintf('%s.FCS',toupper(bf))))
    gb <- gateBeads(bf, 6)
    apc.1 <- as.numeric(gb@exprs[which(gb@labels==1),'APC'])
    q.apc <- quantile(apc.1,probs=seq(0,1,.01))
    mean.apc <- mean(apc.1)
    blank.quantiles <- rbind(blank.quantiles, data.frame(date=flowBeads::getDate(bf), mean.apc=mean.apc, t(q.apc)))
}
colnames(blank.quantiles) 


blank.quantiles <- blank.quantiles[-which(duplicated(blank.quantiles$date)),]
#blank.quantiles[which(duplicated(blank.quantiles$date,fromLast=TRUE)),]


fcsFiles <- unique(unlist(lapply(strsplit(list.files(pattern='cad.*.fcs', path='~nikolas/dunwich/Projects/IL2RA/FCS.Gated/Calli/all.gated'), '-'), function(x) x[[1]])))
d <- data.frame()
for (fcsFile in fcsFiles) { 
    naive.num <- nrow(naive<-read.FCS(naive.fcsFile<-paste(file.path('~nikolas/dunwich/Projects/IL2RA/FCS.Gated/Calli/all.gated',fcsFile), naive.name, sep='-'),channel='CD25',TRANS=log10))
    #cd25.gate <- log10(b[which(as.character(b$date) == getDate(naive.fcsFile)),'mfi.1'])
    cd25.gates <- log10(blank.quantiles[which(as.character(blank.quantiles$date) == getDate(naive.fcsFile)),grep('X',colnames(blank.quantiles))])
    naive.cd25pos.num <- nrow(naive.cd25pos<-read.FCS(naive.cd25pos.fcsFile<-paste(file.path('~nikolas/dunwich/Projects/IL2RA/FCS.Gated/Calli/all.gated',fcsFile),naive.cd25pos.name, sep='-'),channel='CD25',TRANS=log10))
    manual.naive.cd25pos.ratio <- 100*naive.cd25pos.num/naive.num
    bead.naive.cd25pos.ratios <- sapply(cd25.gates, function(cd25.gate) 100*(length(which(naive[,'CD25']>cd25.gate)))/naive.num)
    d <- rbind(d,data.frame(fcsFile=fcsFile,
                            date=getDate(naive.fcsFile),
                            manual.naive.cd25pos.ratio=round(manual.naive.cd25pos.ratio),
                            t(round(bead.naive.cd25pos.ratios)),
                            naive.min.cd25pos=min(naive.cd25pos[,'CD25']),
                            mean.apc=mean.apc,
                            (cd25.gates)))
}
colnames(d) <- c('fcsFile', 'date', 'manual.naive.cd25pos.ratio', paste('q',seq(0,100,1),sep='.'), 'naive.min.cd25pos', 'mean.apc', paste('g',seq(0,100,1),sep='.'))

#best gate position agreement with manual
pct <- 40:99
ms.diff <- sapply(paste('g',pct,sep='.'), function(i) mean((d[,'naive.min.cd25pos']-d[,i])**2))
pdf( '~nikolas/GoogleDrive/PhD/Thesis/IL2RA/figures/cd25pos-gate-agreement.pdf' )
plot(pct, ms.diff,
     ylab='MSD with position of manual gate',
     xlab='APC percentile of blank bead population')
i <- as.numeric(which.min(ms.diff))
print(min.diff.gate <- paste('g',pct,sep='.')[[i]])
points(pct[[i]],ms.diff[[i]], pch=20)
text(pct[[i]], ms.diff[[i]], as.character(pct[[i]]),pos=3)
dev.off()

d$bead.naive.cd25pos.ratio <- d[,gsub('g','q',min.diff.gate)]

#we need this merge to get the individual names
x <- read.csv('~/Projects/IL2RA/Calli_CD25bright_CBR200.csv')
x$fcsFile <- tolower(x$fcsFile)
d$fcsFile <- as.character(d$fcsFile)
print(dim(x <- merge(x, d, by='fcsFile')))
x$fcsFile <- gsub('.fcs','',x$fcsFile)
x$date <- x$date.y
#recalled individuals
r <- read.csv('~/Projects/IL2RA/CellPhenotypes/recalled.individuals.pch') 
#date mismatch always fucks things up on merge
print(dim(r <- merge(r[,c('individual','fcsFile','pch')], x, all.x=TRUE, all.y=FALSE)))
r <- r[order(r$individual,as.Date(r$date)),]
r$pch <- as.character(r$pch)
#get pch for individuals
print(dim(x<-merge(x, r[c(TRUE,FALSE),c('individual','pch')], all.x=TRUE)))
print(table(x$pch <- ifelse(is.na(x$pch),'x',x$pch)))
R <- na.omit(cbind(r[c(TRUE,FALSE),c('individual','pch','manual.naive.cd25pos.ratio','bead.naive.cd25pos.ratio')],r[c(FALSE,TRUE),c('manual.naive.cd25pos.ratio','bead.naive.cd25pos.ratio')]))
colnames(R) <- c("individual", "pch", "manual.naive.cd25pos.ratio.day1", "bead.naive.cd25pos.ratio.day1", "manual.naive.cd25pos.ratio.day2", "bead.naive.cd25pos.ratio.day2")
d <- x


# % naive cd25pos phenotype agreement
pdf( '~nikolas/GoogleDrive/PhD/Thesis/IL2RA/figures/naive-cd25pos-beads-manual-agreement.pdf' )
xlim <- range(d[,c('manual.naive.cd25pos.ratio','bead.naive.cd25pos.ratio')])
main <- round(cor(d$manual.naive.cd25pos.ratio,d$bead.naive.cd25pos.ratio)**2,digits=3)
plot( d$bead.naive.cd25pos.ratio, d$manual.naive.cd25pos.ratio,
     xlab=paste('beads.thresh',PHENO), ylab=paste('manual',PHENO),
     xlim=xlim, ylim=xlim, main=bquote(r^2 == .(main)), pch=d$pch )
abline(b=1,a=0)
dev.off()


#CD25+ gate over time
pdf( '~nikolas/GoogleDrive/PhD/Thesis/IL2RA/figures/cd25pos-gates.pdf' )
d <- d[order(as.Date(d$date)),]
plot(as.Date(d$date), d$naive.min.cd25pos, pch='-', cex=2, xlab='', ylab='CD25+ threshold')
manual.cd25pos.gate.range <- sapply(unique(as.Date(d$date)), function(x) range(d[which(as.Date(d$date)==x),'naive.min.cd25pos']))
segments(x0=unique(as.Date(d$date)), y0=manual.cd25pos.gate.range[1,], y1=manual.cd25pos.gate.range[2,])
lines(as.Date(d$date), d[,min.diff.gate], col='blue')
dev.off()

#repeatability
pdf('~/GoogleDrive/PhD/Thesis/IL2RA/figures/repeatability-cd25pos-naive.pdf')
xlim <- range(R[,grep('cd25',colnames(R))])
plot(R[,'manual.naive.cd25pos.ratio.day1'],R[,'manual.naive.cd25pos.ratio.day2'], pch=as.character(R$pch), xlim=xlim, ylim=xlim,
     xlab=paste('day 1:',PHENO), ylab=paste('day 2:',PHENO), cex=1.5, cex.lab=1.5)
abline(b=1,a=0)
points(R[,'bead.naive.cd25pos.ratio.day1'],R[,'bead.naive.cd25pos.ratio.day2'], pch=as.character(R$pch), col='red', cex=1.5)
rp <- vector('expression',2)
rp[1] <- substitute(expression(r^2 == r2),list(r2=( round(cor(R[,'manual.naive.cd25pos.ratio.day1'],R[,'manual.naive.cd25pos.ratio.day2'])**2,3) )))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=( round(cor(R[,'bead.naive.cd25pos.ratio.day1'],R[,'bead.naive.cd25pos.ratio.day2'])**2,3) )))[2]
legend('bottomright', legend=rp, text.col=c('black','red'))
dev.off()


fcsFile.noisy <- 'cad64_2008mar26_treg_i007576j_017.fcs'
fcsFile.good <- 'cad116_2008oct09_treg_i009546a_013.fcs'
for ( fcsFile in c(fcsFile.noisy, fcsFile.good) ) {
    naive.num <- nrow(naive<-read.FCS(naive.fcsFile<-paste(file.path('~nikolas/dunwich/Projects/IL2RA/FCS.Gated/Calli/all.gated',fcsFile), naive.name, sep='-'),channel='CD25',TRANS=log10))
    #naive
    plot(density(naive))
}

library(ggplot2) 
ggplot(data=x, aes(x=as.Date(x$date), y=x$naive.min.cd25pos)) + geom_point() + geom_smooth(col='black') + geom_point(data=x, aes(x=as.Date(x$date), y=x$g.86), col='blue') + geom_smooth(data=x, aes(x=as.Date(x$date), y=x$g.86), col='blue')

ggplot(data=x, aes(x=as.Date(x$date), y=x$naive.min.cd25pos)) + geom_point() + geom_smooth(col='black') + geom_point(data=x, aes(x=as.Date(x$date), y=x$g.50), col='blue') + geom_smooth(data=x, aes(x=as.Date(x$date), y=x$g.50), col='blue')


#association tests
library(nlme)
library(beeswarm)
rm(box)

x$rs12722495 <- as.numeric(x$DIL9620.CD25.rs12722495)-1
x$rs2104286 <- as.numeric(x$DIL8103.CD25.rs2104286)-1
x$rs11594656 <- as.numeric(x$DIL10847.CD25.rs11594656)-1

f <- function(snp,x) {
    pheno1 <- 'manual.naive.cd25pos.ratio'
    pheno2 <- 'bead.naive.cd25pos.ratio'
    print(dim(x1 <- na.omit(x[,c(pheno1,snp,'individual','pch')])))
    x1$col <- 'black'
    print(dim(x2 <- na.omit(x[,c(pheno2,snp,'individual','pch')])))
    x2$col <- 'red'
    colnames(x1) <- colnames(x2) <- c('pheno',snp,'individual','pch','col')
    print(head(x <- rbind(x1,x2)))
    #x[,pheno] <- log10(x[,pheno])
    #fm <- formula(paste('pheno ~',snp,sep=''))
    fm <- formula(sprintf('pheno ~ %s +  age + sex',snp))
    m1<-lme(fm, random=~ 1|as.factor(individual), data=x1)
    print(p.value.1 <- summary(m1)$tTable[,'p-value'][[2]])
    m2<-lme(fm, random=~ 1|as.factor(individual), data=x2)
    print(p.value.2 <- summary(m2)$tTable[,'p-value'][[2]])
    #m2<-lm(fm, data=x1)
    #print(summary(m2)$coefficients[,"Pr(>|t|)"])
    beeswarm(fm, data=x, pwcol=x$col, pwpch=x$pch, ylab=PHENO)
    #, main=sprintf('p-value=%f',p.value))
    #abline(intervals(m1)$fixed[,1], col=unique(x1$col), lty=2)
    abline(intervals(m1)$fixed[,2], col=unique(x1$col))
    #abline(intervals(m1)$fixed[,3], col=unique(x1$col), lty=2)
    #abline(intervals(m2)$fixed[,1], col=unique(x2$col), lty=2)
    abline(intervals(m2)$fixed[,2], col=unique(x2$col))
    #abline(intervals(m2)$fixed[,3], col=unique(x2$col), lty=2)
    return( list( manual=m1,
                  auto=m2,
                  manual.lower=round(intervals(m1)$fixed[,1][2],digits=3),
                  manual.fixed=round(intervals(m1)$fixed[,2][2],digits=3),
                  manual.upper=round(intervals(m1)$fixed[,3][2],digits=3),
                  manual.pvalue=format.pval(p.value.1),
                  auto.lower=round(intervals(m2)$fixed[,1][2],digits=3),
                  auto.fixed=round(intervals(m2)$fixed[,2][2],digits=3),
                  auto.upper=round(intervals(m2)$fixed[,3][2],digits=3),
                  auto.pvalue=format.pval(p.value.2)
                  ) )
} 


summary(lm(manual.naive.cd25pos.ratio ~ rs12722495 + rs2104286 + rs11594656 + age + sex, data=x))
summary(lm(bead.naive.cd25pos.ratio ~  rs12722495 + rs2104286 + rs11594656 + age + sex, data=x))



m.rs12722495 <- f('rs12722495',x)
m.rs2104286 <- f('rs2104286',x)
m.rs11594656<- f('rs11594656',x)
m.sex <- f('sex',x)
m.age <- f('age',x)

pval.col <- function(pval) if (as.numeric(pval)<.05) sprintf('\\textcolor{red}{%s}',pval) else pval

params <- list(
               #rs12722495
               manual_rs12722495_lower = m.rs12722495$manual.lower,
               manual_rs12722495_effect = m.rs12722495$manual.fixed,
               manual_rs12722495_upper = m.rs12722495$manual.upper,
               manual_rs12722495_pvalue = pval.col(m.rs12722495$manual.pvalue),
               #
               auto_rs12722495_lower   = m.rs12722495$auto.lower,
               auto_rs12722495_effect   = m.rs12722495$auto.fixed,
               auto_rs12722495_upper   = m.rs12722495$auto.upper,
               auto_rs12722495_pvalue   = pval.col(m.rs12722495$auto.pvalue),
               #rs2104286
               manual_rs2104286_lower  = m.rs2104286$manual.lower,
               manual_rs2104286_effect  = m.rs2104286$manual.fixed,
               manual_rs2104286_upper  = m.rs2104286$manual.upper,
               manual_rs2104286_pvalue  = pval.col(m.rs2104286$manual.pvalue),
               #
               auto_rs2104286_lower  = m.rs2104286$auto.lower,
               auto_rs2104286_effect    = m.rs2104286$auto.fixed,
               auto_rs2104286_upper  = m.rs2104286$auto.upper,
               auto_rs2104286_pvalue  = pval.col(m.rs2104286$auto.pvalue),
               #rs11594656
               manual_rs11594656_lower   = m.rs11594656$manual.lower,
               manual_rs11594656_effect   = m.rs11594656$manual.fixed,
               manual_rs11594656_upper   = m.rs11594656$manual.upper,
               manual_rs11594656_pvalue   = pval.col(m.rs11594656$manual.pvalue),
               #
               auto_rs11594656_lower   = m.rs11594656$auto.lower,
               auto_rs11594656_effect   = m.rs11594656$auto.fixed,
               auto_rs11594656_upper   = m.rs11594656$auto.upper,
               auto_rs11594656_pvalue   = pval.col(m.rs11594656$auto.pvalue),
               #age
               manual_age_lower   = m.age$manual.lower,
               manual_age_effect   = m.age$manual.fixed,
               manual_age_upper   = m.age$manual.upper,
               manual_age_pvalue   = pval.col(m.age$manual.pvalue),
               #
               auto_age_lower   = m.age$auto.lower,
               auto_age_effect   = m.age$auto.fixed,
               auto_age_upper   = m.age$auto.upper,
               auto_age_pvalue   = pval.col(m.age$auto.pvalue),
               #sex
               manual_sex_lower   = m.sex$manual.lower,
               manual_sex_effect   = m.sex$manual.fixed,
               manual_sex_upper   = m.sex$manual.upper,
               manual_sex_pvalue   = pval.col(m.sex$manual.pvalue),
               #
               auto_sex_lower   = m.sex$auto.lower,
               auto_sex_effect   = m.sex$auto.fixed,
               auto_sex_upper   = m.sex$auto.upper,
               auto_sex_pvalue   = pval.col(m.sex$auto.pvalue)
               )

library(whisker)
cat(whisker.render(template <- readLines('~/GoogleDrive/PhD/Thesis/IL2RA/cd25pos-assoc-whisker-template.tex'),params))




