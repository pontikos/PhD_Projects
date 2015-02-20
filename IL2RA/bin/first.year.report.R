source('~nikolas/bin/load.cell.phenotypes.R')

manual.beads<-read.csv('~nikolas/Beads.Stats/beads.manual.fcs2.stats')[,-2]
as.character(manual.beads$file.name)->manual.beads$file.name
gsub('_6beads.*$', '', manual.beads$file.name)->manual.beads$file.name

### Recalculate alpha, beta and rse given five bead MFI
log10(c(4100, 10300, 25500, 67300, 139100))->mef
mef.param <- function(beads, i) {
    log10(as.numeric(beads[i,paste("mfi", 2:6, sep=".")]))->mfi
    lm( mef ~ mfi )->m
    return( c(mef.alpha=m$coefficients[[1]], mef.beta=m$coefficients[[2]], mef.rse=as.numeric(summary(m)$sigma)) )
}

t(sapply(1:nrow(manual.beads), function(x) mef.param(manual.beads, x) ))->manual.mef.params
round(manual.mef.params[,'mef.beta'],4)->manual.beads$mef.beta
round(manual.mef.params[,'mef.alpha'],4)->manual.beads$mef.alpha
round(manual.mef.params[,'mef.rse'],4)->manual.beads$mef.rse


auto.beads<-read.csv('~nikolas/Beads.Stats/beads.kmedoids.fcs2.stats')[,-2]
as.character(auto.beads$file.name)->auto.beads$file.name
gsub('_6beads.*$', '', auto.beads$file.name)->auto.beads$file.name

t(sapply(1:nrow(auto.beads), function(x) mef.param(auto.beads, x) ))->auto.mef.params
round(auto.mef.params[,'mef.beta'],4)->auto.beads$mef.beta
round(auto.mef.params[,'mef.alpha'],4)->auto.beads$mef.alpha
round(auto.mef.params[,'mef.rse'],4)->auto.beads$mef.rse



#Comparison of MFI values obtained from auto and manual gating for 57 beads.
#Note that in general the MFI values predicted with auto gating are slightly greater than those predicted
#by manual gating.
plot.manual.auto.mfi <- function(file.name) {
 merge(auto.beads,manual.beads,by="file.name")->m
 pdf(file.name, pointsize=14)
 par( oma = c( 0, 0, 0, 0 ) )
 par(mgp=c(1.75,.5,0), mfrow=c(2,3))
 R.squared <- function(x, y) return( 1-(mean((x-y)**2))/var(c(x,y)) )
 plot(m$mfi.1.x,m$mfi.1.y,xlab="auto MFI",ylab="manual MFI",main="Blank Beads")
 abline(b=1,a=0)
 plot(m$mfi.2.x,m$mfi.2.y,xlab="auto MFI",ylab="manual MFI",main="Beads 2")
 abline(b=1,a=0)
 plot(m$mfi.3.x,m$mfi.3.y,xlab="auto MFI",ylab="manual MFI",main="Beads 3")
 abline(b=1,a=0)
 plot(m$mfi.4.x,m$mfi.4.y,xlab="auto MFI",ylab="manual MFI",main="Beads 4")
 abline(b=1,a=0)
 plot(m$mfi.5.x,m$mfi.5.y,xlab="auto MFI",ylab="manual MFI",main="Beads 5")
 abline(b=1,a=0)
 plot(m$mfi.6.x,m$mfi.6.y,xlab="auto MFI",ylab="manual MFI",main="Beads 6")
 abline(b=1,a=0)
 dev.off()
}
plot.manual.auto.mfi('~nikolas/Plots/FirstYearReport/Agreement/manual-auto-mfi.pdf')



### Agreement of mef parameters estimated using manual gating
plot.auto.manual.agreement <- function(file.name, property, ...) {
    pdf(file.name, pointsize=14)
    par( oma = c( 0, 0, 0, 0 ))
    par(mgp=c(1.75,.5,0))
    range(c(manual.beads[,property], auto.beads[,property]))->xlim
    merge(manual.beads, auto.beads, by=c("file.name"), suffixes=c(".manual", ".auto"))->beads
    plot(beads[,paste(property,"auto",sep=".")], beads[,paste(property,"manual",sep=".")], xlim=xlim, ylim=xlim, ...)
    abline(a=0, b=1)
    R.squared <- function(x, y) return( 1-(mean((x-y)**2))/var(c(x,y)) )
    R.squared(beads[,paste(property,"auto",sep=".")], beads[,paste(property,"manual",sep=".")])->R2
    round(R2, digits=3)->R2
    title( as.expression(bquote( R^2 == .(R2))), line=1 )
    dev.off()
}
plot.auto.manual.agreement('~nikolas/Plots/FirstYearReport/Agreement/mef-beta-agreement.pdf', 'mef.alpha', ylab=expression(paste('manual ', beta)), xlab=expression(paste('auto ', beta)))
plot.auto.manual.agreement('~nikolas/Plots/FirstYearReport/Agreement/mef-alpha-agreement.pdf', 'mef.beta', ylab=expression(paste('manual ', alpha)), xlab=expression(paste('auto ', alpha)))
plot.auto.manual.agreement('~nikolas/Plots/FirstYearReport/Agreement/mef-rse-agreement.pdf', 'mef.rse')



### Effect on repeatability of memory MEF (calli)
calli.memory<-load.individual.cell.phenotypes('~nikolas/CellPhenotypes/calli.memory.stats')$all
as.character(calli.memory$file.name)->calli.memory$file.name
gsub('_i0.*$', '', calli.memory$file.name)->calli.memory$file.name.main

plot.repeatability <- function(method, manual.method, method.name, phenotypes, file.name,...) {
    # 2 measures of repeatability
    # Vincent's calculation
    R.squared <- function(x, y) return( 1-(mean((x-y)**2))/var(c(x,y)) )
    # Pearson's correlation. We use the cor function here.
    r.squared <- function(x, y) return( cov(x,y)/sqrt(var(x)*var(y)) )
    for (p in names(phenotypes)) {
        pdf(file.name, pointsize=14)
        par(mgp=c(1.75,.5,0))
        #par(mfrow=c(1,1), mar=c(2, 4, 2, 2))
        par(mfrow=c(1,1))
        method[,paste(p, "day1", sep=".")]->pheno.day1
        method[,paste(p, "day2", sep=".")]->pheno.day2

        manual.method[,paste(p, "day1", sep=".")]->manual.pheno.day1
        manual.method[,paste(p, "day2", sep=".")]->manual.pheno.day2

        range(manual.pheno.day1, pheno.day1, manual.pheno.day2, pheno.day2)->pheno.range

        plot(pheno.day1, pheno.day2, pch=method$pch,
             xlab=paste(phenotypes[[p]], "Day 1"),
             ylab=paste(phenotypes[[p]], "Day 2"),
             xlim=pheno.range,
             ylim=pheno.range,
             col="red", cex=.75)
        round(cor(pheno.day1, pheno.day2), digits=3)->r2
        print(r2)
        round(R.squared(pheno.day1, pheno.day2), digits=3)->R2
        print(R2)

        points(manual.pheno.day1, manual.pheno.day2, pch=manual.method$pch, col="black", cex=.75)
        round(cor(manual.pheno.day1, manual.pheno.day2), digits=3)->manual.r2
        print(manual.r2)
        round(R.squared(manual.pheno.day1, manual.pheno.day2), digits=3)->manual.R2
        print(manual.R2)
        abline(a=0,b=1)
        dev.off()
    }
}

merge(calli.memory, manual.beads, by.x=c('file.name.main'), by.y=c('file.name'))->d
(10**d$mef.beta)*(d$APC.mean**d$mef.alpha)->d$APC.mef
merge(RECALLED.INDIVIDUALS.PCH, d[,-4], by=c("individual", "fcsFile"))->m
(recalled.phenotypes.table((m[!duplicated(as.character(m$fcsFile)),])))->manual.recalled

merge(calli.memory, auto.beads, by.x=c('file.name.main'), by.y=c('file.name'))->d
(10**d$mef.beta)*(d$APC.mean**d$mef.alpha)->d$APC.mef
merge(RECALLED.INDIVIDUALS.PCH, d[,-4], by=c("individual", "fcsFile"))->m
(recalled.phenotypes.table((m[!duplicated(as.character(m$fcsFile)),])))->auto.recalled

R.squared <- function(x, y) return( 1-(mean((x-y)**2))/var(c(x,y)) )
print( R.squared(manual.recalled$APC.mef.day1, manual.recalled$APC.mef.day2) )
print( R.squared(auto.recalled$APC.mef.day1, auto.recalled$APC.mef.day2) )

plot.repeatability(method=auto.recalled,
manual.method=manual.recalled,
phenotypes=list(APC.mef="Memory CD25 MEF"), file.name="~/Plots/FirstYearReport/Repeatability/mef-parameters-repeatability.pdf")




###

histograms <- function(d, cell.phenotype, file.name) {
    pdf(file.name)
    par(mfrow=c(1,2))
    hist((d[,cell.phenotype]), xlim=c(0,100), breaks=50, main="")
    hist(log10(d[,cell.phenotype]),  breaks=50, main="")
    dev.off()
}

histograms( load.individual.cell.phenotypes("~nikolas/CellPhenotypes/calli.cell.phenotypes")$unique, cell.phenotype="memory.freqpar", file.name="~/Plots/FirstYearReport/Effects/manual-memory-freqpar-histograms.pdf" )
histograms( load.individual.cell.phenotypes("~nikolas/CellPhenotypes/mm.cell.phenotypes")$unique, cell.phenotype="memory.freqpar", file.name="~/Plots/FirstYearReport/Effects/mm-memory-freqpar-histograms.pdf" )
histograms( load.individual.cell.phenotypes("~nikolas/CellPhenotypes/sp.mm.cell.phenotypes")$unique, cell.phenotype="memory.freqpar", file.name="~/Plots/FirstYearReport/Effects/sp-mm-memory-freqpar-histograms.pdf" )
histograms( load.individual.cell.phenotypes("~nikolas/CellPhenotypes/calli.cell.phenotypes")$unique, cell.phenotype="naive.cd25pos.freqpar", file.name="~/Plots/FirstYearReport/Effects/manual-naive-freqpar-histograms.pdf" )
histograms( load.individual.cell.phenotypes("~nikolas/CellPhenotypes/naive.cd25pos.beads.cell.phenotypes")$unique, cell.phenotype="naive.cd25pos.freqpar", file.name="~/Plots/FirstYearReport/Effects/beads-naive-freqpar-histograms.pdf" )

###


### MEF regression
#figure:MEF
plot.MEF <- function(mfi2, file.name) {
    pdf(file.name)
    mef <- c(4100, 10300, 25500, 67300, 139100)
    line(log10(mfi2[-1]),log10(mef))->l2
    coef(l2)->c2
    print(c2)
    MEF<- function(x) return (c2[1]+log10(x)*c2[2])
    plot(log10(mfi2), MEF(mfi2),pch=".",
         xlab=expression(paste(plain(Log[10]),'(MFI)')),
         ylab=expression(paste(plain(Log[10]),'(MEF)'))
         )
         #xlim=c(0,4),ylim=c(0,5))
    points(log10(mfi2[-1]),log10(mef))
    abline(c2)
    abline(h=MEF(mfi2),lty=3)
    draw.line <- function(x, ...) {
        log10(x)->x0
        MEF(x) -> y1
        segments(x0=x0,y0=0,x1=x0,y1=y1,...)
    }
    calli.memory <- load.individual.cell.phenotypes("~nikolas/CellPhenotypes/calli.cell.phenotypes")$all
    draw.line(mfi2[1], col="black", lty=2)
    draw.line(x=min(calli.memory$APC.mfi), col="red")
    draw.line(x=max(calli.memory$APC.mfi), col="green")
    dev.off()
}
beads.manual.fcs2.stats<-manual.beads
plot.MEF(
mfi2=as.numeric(beads.manual.fcs2.stats[which(beads.manual.fcs2.stats$file.name=="cad57_2008feb20_treg"),grep("mfi.*", names(beads.manual.fcs2.stats))]),
file.name="~/Plots/FirstYearReport/BeadNormalisation/MEF.pdf"
)



### cd25pos gate over time
pdf('~/Plots/FirstYearReport/cd25pos-gates.pdf')
read.csv('~/Beads.Stats/beads.fcs2.kmedoids.stats') -> beads.kmedoids.stats
beads.kmedoids.stats[order(as.Date(beads.kmedoids.stats$date)),]->beads.kmedoids.stats
beads.kmedoids.stats$mfi.1->mfi.1
beads.kmedoids.stats$p95fi.1->p95fi.1
### cd25.thresholds obtained from gate.boundaries.R
sapply(cd25.gates, function(x) min((x[,1]))) -> cd25.thresholds
#ylim <- range(mfi.1, p95fi.1, cd25.thresholds) 
ylim <- range(p95fi.1, cd25.thresholds) 
plot(as.Date(d$date), cd25.thresholds, pch="-", ylab="CD25+ Threshold", xlab="", main="CD25+ Gate Over Time", ylim=ylim)
t(sapply(tapply(cd25.thresholds, d$date, range), function(x) x))->s
segments(x0=as.Date(rownames(s),format="%Y-%m-%d"), y0=s[,1], y1=s[,2])
#abline(v=as.Date('2008-03-26'), lty=3)
#abline(v=as.Date('2008-04-16'), lty=3)
lines(as.Date(beads.kmedoids.stats$date), p95fi.1, col="red")
#lines(as.Date(beads.kmedoids.stats$date), mfi.1, col="red")
#lines(as.Date(beads.kmedoids.stats$date), mfi.1+2*beads.kmedoids.stats$sdfi.1, col="red", lty=2)
dev.off()

#confint
load.individual.cell.phenotypes("~nikolas/CellPhenotypes/calli.cell.phenotypes")$unique->d
(round(confint(model(d$memory.freqpar,d)), digits=2))
load.individual.cell.phenotypes("~nikolas/CellPhenotypes/mm.cell.phenotypes")$unique->d
(round(confint(model(d$memory.freqpar,d)), digits=2))
load.individual.cell.phenotypes("~nikolas/CellPhenotypes/sp.mm.cell.phenotypes")$unique->d
(round(confint(model(d$memory.freqpar,d)), digits=2))
load.individual.cell.phenotypes("~nikolas/CellPhenotypes/mm.daily.mean.cell.phenotypes.2")$unique->d
(round(confint(model(d$memory.freqpar,d)), digits=2))

#snp95 table
latex.row <- function(pheno, d, r) {
    model(d[,pheno],d)->m
    round(as.numeric(confint(m)[r,]), digits=3)->conf
    #effect, confint, p-value
    as.numeric(coef(summary(m))[r,]) ->co
    cat( "$", round(co[[1]], digits=3), "$ & $", "[", conf[[1]], ";", conf[[2]], "]", "$ & $", co[[4]], "$\n", sep="" )
}

latex.row("APC.mef", load.individual.cell.phenotypes("~nikolas/CellPhenotypes/calli.cell.phenotypes")$unique)
latex.row("APC.mef", load.individual.cell.phenotypes("~nikolas/CellPhenotypes/mm.cell.phenotypes")$unique)
latex.row("APC.mef", load.individual.cell.phenotypes("~nikolas/CellPhenotypes/sp.mm.cell.phenotypes")$unique)
latex.row("APC.mef", load.individual.cell.phenotypes("~nikolas/CellPhenotypes/mm.daily.mean.cell.phenotypes.2")$unique)
#

latex.row("memory.freqpar", load.individual.cell.phenotypes("~nikolas/CellPhenotypes/calli.cell.phenotypes")$unique)
latex.row("memory.freqpar", load.individual.cell.phenotypes("~nikolas/CellPhenotypes/mm.cell.phenotypes")$unique)
latex.row("memory.freqpar", load.individual.cell.phenotypes("~nikolas/CellPhenotypes/sp.mm.cell.phenotypes")$unique)
latex.row("memory.freqpar", load.individual.cell.phenotypes("~nikolas/CellPhenotypes/mm.daily.mean.cell.phenotypes.2")$unique)
#
latex.row("memory.freqpar", load.individual.cell.phenotypes("~nikolas/CellPhenotypes/fixed.gate.cell.phenotypes")$unique, 2)

latex.row("naive.cd25pos.freqpar", load.individual.cell.phenotypes("~nikolas/CellPhenotypes/calli.cell.phenotypes")$unique, 2)
latex.row("naive.cd25pos.freqpar", load.individual.cell.phenotypes("~nikolas/CellPhenotypes/calli.cell.phenotypes")$unique, 3)
latex.row("naive.cd25pos.freqpar", load.individual.cell.phenotypes("~nikolas/CellPhenotypes/calli.cell.phenotypes")$unique, 4)
latex.row("naive.cd25pos.freqpar", load.individual.cell.phenotypes("~nikolas/CellPhenotypes/calli.cell.phenotypes")$unique, 5)
latex.row("naive.cd25pos.freqpar", load.individual.cell.phenotypes("~nikolas/CellPhenotypes/calli.cell.phenotypes")$unique, 6)

latex.row("naive.cd25pos.freqpar", load.individual.cell.phenotypes("~nikolas/CellPhenotypes/naive.cd25pos.beads.cell.phenotypes")$unique, 2)
latex.row("naive.cd25pos.freqpar", load.individual.cell.phenotypes("~nikolas/CellPhenotypes/naive.cd25pos.beads.cell.phenotypes")$unique, 3)
latex.row("naive.cd25pos.freqpar", load.individual.cell.phenotypes("~nikolas/CellPhenotypes/naive.cd25pos.beads.cell.phenotypes")$unique, 4)
latex.row("naive.cd25pos.freqpar", load.individual.cell.phenotypes("~nikolas/CellPhenotypes/naive.cd25pos.beads.cell.phenotypes")$unique, 5)
latex.row("naive.cd25pos.freqpar", load.individual.cell.phenotypes("~nikolas/CellPhenotypes/naive.cd25pos.beads.cell.phenotypes")$unique, 6)

###

paste(d$snp95, d$snp86, d$snp56,sep=".")->d$genotype

# 7 haplotypes

boxplot(memory.freqpar  ~ genotype, data=calli, at = 0:6*3 + 1, xlim = c(0, 21), ylim = range(calli$memory.freqpar, mm$memory.freqpar), xaxt = "n", col="blue", main="mm")
boxplot(memory.freqpar  ~ genotype, data=mm, at = 0:6*3 + 2, xaxt = "n", add = TRUE, col="red")
axis(1, at = 0:6*3 + 1.5, labels = as.character(levels(as.factor(mm$genotype))), tick = TRUE)

boxplot(memory.freqpar  ~ genotype, data=calli, at = 0:6*3 + 1, xlim = c(0, 21), ylim = range(calli$memory.freqpar, sp.mm$memory.freqpar), xaxt = "n", col="blue", main="sp.mm")
boxplot(memory.freqpar  ~ genotype, data=sp.mm, at = 0:6*3 + 2, xaxt = "n", add = TRUE, col="red")
axis(1, at = 0:6*3 + 1.5, labels = as.character(levels(as.factor(sp.mm$genotype))), tick = TRUE)

boxplot(memory.freqpar  ~ genotype, data=calli, at = 0:6*3 + 1, xlim = c(0, 21), ylim = range(calli$memory.freqpar, learned.mm$memory.freqpar), xaxt = "n", col="blue", main="learned.mm")
boxplot(memory.freqpar  ~ genotype, data=learned.mm, at = 0:6*3 + 2, xaxt = "n", add = TRUE, col="red")
axis(1, at = 0:6*3 + 1.5, labels = as.character(levels(as.factor(sp.mm$genotype))), tick = TRUE)

boxplot(memory.freqpar  ~ genotype, data=mm, at = 0:6*3 + 1, xlim = c(0, 21), ylim = range(mm$memory.freqpar, learned.mm$memory.freqpar), xaxt = "n", col="pink", main="mm vs learned.mm")
boxplot(memory.freqpar  ~ genotype, data=learned.mm, at = 0:6*3 + 2, xaxt = "n", add = TRUE, col="red")
axis(1, at = 0:6*3 + 1.5, labels = as.character(levels(as.factor(sp.mm$genotype))), tick = TRUE)



boxplot(naive.cd25pos.freqpar  ~ genotype, data=calli, at = 0:6*3 + 1, xlim = c(0, 21), ylim = range(calli$naive.cd25pos.freqpar, naive.cd25pos.beads$naive.cd25pos.freqpar), xaxt = "n", col="blue")
boxplot(naive.cd25pos.freqpar  ~ genotype, data=naive.cd25pos.beads, at = 0:6*3 + 2, xaxt = "n", add = TRUE, col="red")
axis(1, at = 0:6*3 + 1.5, labels = as.character(levels(as.factor(calli$genotype))), tick = TRUE)


read.csv('~/CellPhenotypes/info.calli.cell.phenotypes')->d
sapply(1:nrow(d), function(r) mean(read.FCS(sprintf('~/dunwich/FCS.Gated/Calli/NON_T_REGS/%s.fcs', as.character(d[r,'fcsFile'])))@exprs[,5]))->nontreg.CD25.mfi

boxplot(memory.freqpar  ~ snp95, data=calli, at = 0:2*3 + 1, xlim = c(0,9), ylim = range(calli$memory.freqpar, mm$memory.freqpar), xaxt = "n", col="blue", main="mm snp95")
boxplot(memory.freqpar  ~ snp95, data=mm, at = 0:2*3 + 2, xaxt = "n", add = TRUE, col="red")
axis(1, at = 0:2*3 + 1.5, labels = as.character(levels(as.factor(mm$snp95))), tick = TRUE)

boxplot(memory.freqpar  ~ genotype, data=calli, at = 0:6*3 + 1, xlim = c(0, 21), ylim = range(calli$memory.freqpar, sp.mm$memory.freqpar), xaxt = "n", col="blue", main="sp.mm")
boxplot(memory.freqpar  ~ genotype, data=sp.mm, at = 0:6*3 + 2, xaxt = "n", add = TRUE, col="red")
axis(1, at = 0:6*3 + 1.5, labels = as.character(levels(as.factor(sp.mm$genotype))), tick = TRUE)

library(boot)
n<-2
boxplot((logit(naive.cd25pos.freqpar/100))  ~ snp95, data=calli, at = 0:n*3 + 1, xlim = c(0, 3*(n+1)), xaxt = "n", col="blue")
boxplot((logit(naive.cd25pos.freqpar/100))  ~ snp95, data=naive.cd25pos.beads, at = 0:n*3 + 2, xaxt = "n", add = TRUE, col="red")
axis(1, at = 0:n*3 + 1.5, labels = as.character(levels(as.factor(calli$snp95))), tick = TRUE)

library(boot)
n<-2
boxplot((logit(naive.cd25pos.freqpar/100))  ~ snp86, data=calli, at = 0:n*3 + 1, xlim = c(0, 3*(n+1)), xaxt = "n", col="blue")
boxplot((logit(naive.cd25pos.freqpar/100))  ~ snp86, data=naive.cd25pos.beads, at = 0:n*3 + 2, xaxt = "n", add = TRUE, col="red")
axis(1, at = 0:n*3 + 1.5, labels = as.character(levels(as.factor(calli$snp86))), tick = TRUE)

library(boot)
n<-2
boxplot((logit(naive.cd25pos.freqpar/100))  ~ snp56, data=calli, at = 0:n*3 + 1, xlim = c(0, 3*(n+1)), xaxt = "n", col="blue")
boxplot((logit(naive.cd25pos.freqpar/100))  ~ snp56, data=naive.cd25pos.beads, at = 0:n*3 + 2, xaxt = "n", add = TRUE, col="red")
axis(1, at = 0:n*3 + 1.5, labels = as.character(levels(as.factor(calli$snp56))), tick = TRUE)



n<-1
boxplot((logit(naive.cd25pos.freqpar/100))  ~ sex, data=calli, at = 0:n*3 + 1, xlim = c(0, 3*(n+1)), xaxt = "n", col="blue")
boxplot((logit(naive.cd25pos.freqpar/100))  ~ sex, data=naive.cd25pos.beads, at = 0:n*3 + 2, xaxt = "n", add = TRUE, col="red")
axis(1, at = 0:n*3 + 1.5, labels = as.character(levels(as.factor(calli$sex))), tick = TRUE)






read.cell.phenotype <- function(f) {
    read.csv(f)->d
    levels(d$snp95) <- c(0, 1, 2)
    d$snp95<-as.numeric(levels(d$snp95))[d$snp95]
    levels(d$snp86) <- c(0, 1, 2)
    d$snp86<-as.numeric(levels(d$snp86))[d$snp86]
    levels(d$snp56) <- c(2, 1, 0)
    d$snp56<-as.numeric(levels(d$snp56))[d$snp56]
    return(d)
}

#
#calli
read.cell.phenotype('info.calli.cell.phenotypes')->d
#sp.mm
read.cell.phenotype('info.sp.mm.cell.phenotypes')->d
#learned.mm
read.cell.phenotype('info.learned.mm.cell.phenotypes')->d
#mm
read.cell.phenotype('info.mm.cell.phenotypes')->d

#beads.auto
read.csv('info.naive.cd25pos.beads.cell.phenotypes')->d

load.individual.cell.phenotypes("~nikolas/CellPhenotypes/fixed.gate.cell.phenotypes")$unique->d


pheno <- d$memory.freqpar/100

library(boot)
trans <- logit
untrans <- inv.logit


summary(lm(trans(pheno) ~ d$snp95 + d$snp86 + d$snp56 + d$age + d$sex)) -> m


summary(lm(logit(memory.freqpar/100) ~ snp95 + snp86 + snp56 + age + sex, data=read.cell.phenotype('info.calli.cell.phenotypes')))
summary(lm(logit(memory.freqpar/100) ~ snp95 + snp86 + snp56 + age + sex, data=read.cell.phenotype('info.mm.cell.phenotypes')))
summary(lm(logit(memory.freqpar/100) ~ snp95 + snp86 + snp56 + age + sex, data=read.cell.phenotype('info.mm.cell.phenotypes')[c(-61, -72, -79, -102),]))
summary(lm(logit(memory.freqpar/100) ~ snp95 + snp86 + snp56 + age + sex, data=read.cell.phenotype('info.sp.mm.cell.phenotypes')))
summary(lm(logit(memory.freqpar/100) ~ snp95 + snp86 + snp56 + age + sex, data=read.cell.phenotype('info.sp.mm.cell.phenotypes')[c(-35, -56, -102, -42, -72, -77),]))

varnames <- c('intercept', 'rs12722495', 'rs2104286', 'rs11594656', 'Age/10', 'Male')
f <- function(x) return(x/10)
pdf( "~nikolas/thor/Plots/FirstYearReport/Effects/memory-freqpar-association.pdf" )
#library(arm)
coefplot(lm(logit(memory.freqpar/100) ~ snp95 + snp86 + snp56 + f(age) + sex, data=read.cell.phenotype('info.calli.cell.phenotypes')), col.pts='black', main="", vertical=F, var.las=1, CI=2, intercept=F, varnames=varnames)
coefplot(lm(logit(memory.freqpar/100) ~ snp95 + snp86 + snp56 + f(age) + sex, data=read.cell.phenotype('info.mm.cell.phenotypes')), col.pts='blue', add=T, vertical=F, CI=2, varnames=varnames, intercept=F)
coefplot(lm(logit(memory.freqpar/100) ~ snp95 + snp86 + snp56 + f(age) + sex, data=read.cell.phenotype('info.sp.mm.cell.phenotypes')), col.pts='red', add=T, offset=.2, vertical=F, CI=2, varnames=varnames, intercept=F)
dev.off()
pdf( "~nikolas/thor/Plots/FirstYearReport/Effects/memory-freqpar-association-2.pdf" )
library(arm)
coefplot(lm(logit(memory.freqpar/100) ~ snp95 + snp86 + snp56 + f(age) + sex, data=read.cell.phenotype('info.calli.cell.phenotypes')), col.pts='black', vertical=F, var.las=1, main="", CI=2, varnames=varnames, intercept=F)
coefplot(lm(logit(memory.freqpar/100) ~ snp95 + snp86 + snp56 + f(age) + sex, data=read.cell.phenotype('info.learned.mm.cell.phenotypes')), col.pts='red', vertical=F, add=T, CI=2, varnames=varnames, intercept=F)
dev.off()
pdf( "~nikolas/thor/Plots/FirstYearReport/Effects/naive-cd25pos-freqpar-association.pdf" )
library(arm)
coefplot(lm(logit(naive.cd25pos.freqpar/100) ~ snp95 + snp86 + snp56 + f(age) + sex, data=read.cell.phenotype('info.calli.cell.phenotypes')), col.pts='black', vertical=F, var.las=1, main="", CI=2, varnames=varnames, intercept=F)
coefplot(lm(logit(naive.cd25pos.freqpar/100) ~ snp95 + snp86 + snp56 + f(age) + sex, data=read.cell.phenotype('info.naive.cd25pos.beads.cell.phenotypes')), col.pts="red", vertical=F, add=T, CI=2, varnames=varnames, intercept=F)
dev.off()



summary(lm(memory.freqpar ~ snp95, data=read.cell.phenotype('info.calli.cell.phenotypes')))



summary(lm(logit(naive.cd25pos.freqpar/100) ~ snp95 + snp86 + snp56 + age + sex, data=read.cell.phenotype('info.naive.cd25pos.beads.cell.phenotypes')))
summary(lm(logit(naive.cd25pos.freqpar/100) ~ snp95 + snp86 + snp56 + age + sex, data=read.cell.phenotype('info.calli.cell.phenotypes')))


read.cell.phenotype('info.calli.cell.phenotypes')->d



summary(lm(logit(memory.freqpar/100) ~ snp95 + snp86 + snp56 + age + sex, data=d[c(-61, -72, -79, -102),]))


pheno <- d$naive.cd25pos.freqpar/100
pheno <- d$APC.mef




untrans <- inv.logit

m$df[2]->degrees.freedom
print(degrees.freedom)

#throw away the intercept, keep the other coefficients
#we obtain the estimate and the standard error
#estimates
coefficients(m)[-1,"Estimate"]->estimates
effect.names <- c("snp95","snp86","snp56","Age","Male")
names(estimates) <- effect.names
print(estimates)

#standard errors
coefficients(m)[-1,"Std. Error"]->std.errors
names(std.errors) <- effect.names
print(std.errors)

for (e in effect.names) {
    estimates[[e]] -> y.bar
    std.errors[[e]]*qt(.975,df=degrees.freedom)->ci
    untrans(cbind(y.bar-ci, y.bar+ci))->y
    print(e)
    print(y.bar)
    print(y)
}



summary(lm(pheno ~ d$snp95 + d$snp86 + d$snp56 + d$age + d$sex))


summary(lm(log(pheno) ~ d$snp95 + d$snp86 + d$snp56 + d$age + d$sex))



###

plot.agreement(method=list("beads.auto"= load.individual.cell.phenotypes('~nikolas/CellPhenotypes/naive.cd25pos.beads.cell.phenotypes')$all, "manual"=load.individual.cell.phenotypes('~nikolas/CellPhenotypes/calli.cell.phenotypes')$all), phenotypes=list( naive.cd25pos.freqpar="Naive CD25+ % "), file.name="~/Plots/FirstYearReport/Agreement/naive-cd25pos-beads-manual-agreement.pdf", mfrow=c(1,1), width=5, height=5, pheno.range=c(0,60))

plot.agreement(method=list("mm"= load.individual.cell.phenotypes('~nikolas/CellPhenotypes/mm.cell.phenotypes')$all, "manual"=load.individual.cell.phenotypes('~nikolas/CellPhenotypes/calli.cell.phenotypes')$all), phenotypes=list(  memory.freqpar="Memory % "), file.name="~/Plots/FirstYearReport/Agreement/mm-manual-agreement.pdf", mfrow=c(1,1), width=5, height=5, pheno.range=c(0,100))
plot.agreement(method=list("sp.mm"= load.individual.cell.phenotypes('~nikolas/CellPhenotypes/sp.mm.cell.phenotypes')$all, "manual"=load.individual.cell.phenotypes('~nikolas/CellPhenotypes/calli.cell.phenotypes')$all), phenotypes=list(  memory.freqpar="Memory % "), file.name="~/Plots/FirstYearReport/Agreement/sp-mm-manual-agreement.pdf", mfrow=c(1,1), width=5, height=5, pheno.range=c(0,100))
plot.agreement(method=list("learned.mm"= load.individual.cell.phenotypes('~nikolas/CellPhenotypes/learned.mm.cell.phenotypes')$all, "manual"=load.individual.cell.phenotypes('~nikolas/CellPhenotypes/calli.cell.phenotypes')$all), phenotypes=list(  memory.freqpar="Memory % "), file.name="~/Plots/FirstYearReport/Agreement/learned-mm-manual-agreement.pdf", mfrow=c(1,1), width=5, height=5, pheno.range=c(0,100))




####

library(sp)

#read.FCS('~/dunwich/FCS/CAD67_2008apr16_Treg_I007693L_009.fcs')[,-9] -> fcs.data
#c(1,2,3,5,7,8) -> p
#fcs.data[rowSums(fcs.data@exprs[,p]>1)==dim(fcs.data@exprs[,p])[[2]],]->fcs.data

read.FCS('~/dunwich/FCS/CAD67_2008apr16_Treg_I007693L_009.fcs')-> fcs.data

load('~/dunwich/WSP.R.obj/CAD67_160408_Treg.wsp.Rdata')

flowlist@flowJoData[[1]]->f1

exprs <- function(fcs.data) {
    fcs.data@exprs -> x
    log10(x[,3:8]) -> x[,3:8]
    return(x)
}

# FSC-A, SSC-A, Alexa-488-A, PE-Cy7-A, APC-A, PE-A, Alexa-700-A, Pacific Blue-A
exprs(fcs.data)->x

# FSC-A, SSC-A
cbind(x, 'lymph'=0) -> x
f1@filter[grep('CAD67_2008apr16_Treg_I007693L_009.fcs:Lymphocytes$', f1@filterName)][[1]]@boundaries->gate
n<-colnames(gate)
x[,'lymph'] <- as.numeric(point.in.polygon(x[,n[[1]]], x[,n[[2]]], gate[,1], gate[,2]))
100* sum(x[,'lymph']==1)/ (sum(x[,'lymph']==1)+sum(x[,'lymph']==0))

#Alexa-700-A      SSC-A
cbind(x, 'cd4'=0) -> x
trans1 <- log10
trans2 <- identity
f1@filter[grep('CAD67_2008apr16_Treg_I007693L_009.fcs:Lymphocytes:CD4$', f1@filterName)][[1]]@filters[[2]]@boundaries->gate
trans1(gate[,1])->gate[,1]
trans2(gate[,2])->gate[,2]
n<-colnames(gate)
x[,'cd4'] <- as.numeric(point.in.polygon(x[,n[[1]]], x[,n[[2]]], gate[,1], gate[,2]))
100*sum(rowSums(x[,c('cd4','lymph')])==2)/sum((x[,c('lymph')])==1)

#APC-A Alexa-488-A
cbind(x, 'nontregs'=0) -> x
trans1 <- log10
trans2 <- log10
f1@filter[grep("CAD67_2008apr16_Treg_I007693L_009.fcs:Lymphocytes:CD4:CD127hi:127hi$", f1@filterName)][[1]]@filters[[2]]@boundaries -> gate
trans1(gate[,1])->gate[,1]
trans2(gate[,2])->gate[,2]
n<-colnames(gate)
x[,'nontregs'] <- as.numeric(point.in.polygon(x[,n[[1]]], x[,n[[2]]], gate[,1], gate[,2]))
100*sum(rowSums(x[,c('nontregs','cd4','lymph')])==3)/sum(rowSums(x[,c('cd4','lymph')])==2)
100*sum(rowSums(x[,c('nontregs','cd4','lymph')])==3)/sum((x[,c('lymph')])==1)
100*sum(rowSums(x[,c('nontregs','cd4','lymph')])==3)/sum((x[,c('cd4')])==1)

#Pacific Blue-A
cbind(x, 'memory'=0) -> x
log10(f1@filter[grep("CAD67_2008apr16_Treg_I007693L_009.fcs:Lymphocytes:CD4:CD127hi:127hi:CD45RAneg$", f1@filterName)][[1]]@filters[[2]]@max)->memory.gate
x[,'memory'] <- as.numeric(x[,'Pacific Blue-A']<=memory.gate)

#Pacific Blue-A
cbind(x, 'naive'=0) -> x
log10(f1@filter[grep("CAD67_2008apr16_Treg_I007693L_009.fcs:Lymphocytes:CD4:CD127hi:127hi:CD45RApos$", f1@filterName)][[1]]@filters[[2]]@min)->naive.gate
x[,'naive'] <- as.numeric(x[,'Pacific Blue-A']>=naive.gate)

#APC-A
cbind(x, 'cd25pos'=0) -> x
log10(min(f1@filter[grep('CAD67_2008apr16_Treg_I007693L_009.fcs:Lymphocytes:CD4:CD127hiCD25pos$', f1@filterName)][[1]]@filters[[2]]@boundaries[,1])) -> cd25pos.gate
x[,'cd25pos'] <- as.numeric(x[,'APC-A']>=cd25pos.gate)


#how many cd4 over all cells?
100* sum(x[,'cd4']==1)/ length(x[,'lymph'])
#how many cd4+lymph over all cells?
100*sum(rowSums(x[,c('cd4','lymph')])==2)/ length(x[,'lymph'])
#ratio
sum(x[,'cd4']==1)/sum(rowSums(x[,c('cd4','lymph')])==2)


sum( rowSums(x[,c('nontregs','cd4','lymph', 'memory')])==4 )


read.FCS('~/dunwich/FCS.Marcin/CAD67_2008apr16_Treg_I007693L_009.fcs')[,-9] -> fcs3
compensate(fcs3,fcs3@description[["SPILL"]])->fcs3.comp
lgcl <- estimateLogicle(fcs3.comp, channels = c( 'APC-A', 'Pacific Blue-A' ))
fcs3.trans <- transform(fcs3.comp, lgcl)

library(flowViz)
pdf('Plots/rplot.pdf')
plot(fcs3.trans@exprs[,c('APC-A', 'Pacific Blue-A')], pch=".", col=rgb(0,0,0,.25), xlim=c(-.5,3.5),ylim=c(-.25,3))
#plot(fcs3.trans, c('APC-A', 'Pacific Blue-A'))
#contour(fcs3.trans, c('APC-A', 'Pacific Blue-A'), add=T)
points(fcs3.trans@exprs[which(rowSums(x[,c('nontregs','cd4','lymph', 'naive')])==4), c('APC-A', 'Pacific Blue-A')], pch='.', col=rgb(.25,1,.25,.9))
points(fcs3.trans@exprs[which(rowSums(x[,c('nontregs','cd4','lymph', 'naive', 'cd25pos')])==5), c('APC-A', 'Pacific Blue-A')], pch='.', col='red')
points(fcs3.trans@exprs[which(rowSums(x[,c('nontregs','cd4','lymph', 'memory')])==4), c('APC-A', 'Pacific Blue-A')], pch='.', col=rgb(.25,.25,1,.9))
points(fcs3.trans@exprs[which(rowSums(x[,c('cd4','lymph', 'cd25pos')])==3 & x[,'nontregs']==0), c('APC-A', 'Pacific Blue-A')], pch='.', col='pink')
dev.off()



fcs.data-> fcs2
library(flowViz)
pdf('Plots/rplot2.pdf')
plot(log10(fcs2@exprs[,c('APC-A', 'Pacific Blue-A')]), pch=".", col=rgb(0,0,0,.25), xlim=c(0,3.5),ylim=c(0,3))
#contour(fcs3.trans, c('APC-A', 'Pacific Blue-A'), add=T)
points(log10(fcs2@exprs[which(rowSums(x[,c('nontregs','cd4','lymph', 'naive')])==4), c('APC-A', 'Pacific Blue-A')]), pch='.', col=rgb(.25,1,.25,.9))
points(log10(fcs2@exprs[which(rowSums(x[,c('nontregs','cd4','lymph', 'naive', 'cd25pos')])==5), c('APC-A', 'Pacific Blue-A')]), pch='.', col='red')
points(log10(fcs2@exprs[which(rowSums(x[,c('nontregs','cd4','lymph', 'memory')])==4), c('APC-A', 'Pacific Blue-A')]), pch='.', col=rgb(.25,.25,1,.9))
points(log10(fcs2@exprs[which(rowSums(x[,c('cd4','lymph', 'cd25pos')])==3 & x[,'nontregs']==0), c('APC-A', 'Pacific Blue-A')]), pch='.', col='pink')
abline(h=memory.gate, col='red')
abline(h=naive.gate, col='red')
abline(v=cd25pos.gate, col='red')
dev.off()


####
# beads

library(flowCore)
library(flowClust)
library(cluster)

toupper(as.character(read.csv('~nikolas/Beads.Stats/mismatch.fcs', header=FALSE)[,1])) -> fcs

for (f in fcs[13:26]) {

#read.FCS('~/dunwich/FCS.Beads/FCS2.0/CAD111_2008SEP11_TREG_6BEADS_110908.FCS') -> fcs.data
#read.FCS(sprintf('~/dunwich/FCS.Beads/FCS2.0/%s.FCS', f)) -> fcs.data
read.FCS('~/dunwich/FCS.Beads/FCS2.0/CAD89_2008JUL03_TREG_6BEADS_030708.FCS') -> fcs.data

flowClust(fcs.data, varNames=c("FSC-A","SSC-A"), K=1)->res 
fcs.data@exprs[which(Map(res)==1),5]->x

#work with a subset of x which only contains >1
i <- x>1
x <- log10(x[i])


pam.res <- pam(x, k=6)
table(pam.res$clustering)
sort(tapply(x, pam.res$clustering, mean))

res1 <- kmeans(x, centers=6, nstart=2000)
table(res1$cluster)
sort(tapply(x, res1$cluster, mean))

res2 <- kmeans(x, centers=pam.res$medoids)
table(res2$cluster)
sort(tapply(x, res2$cluster, mean))

sort(tapply(x, pam.res$clustering, mean))
sort(tapply(x, res1$cluster, mean))


}

####

library(mixtools)
par(mfrow=c(3,3))
fcs.files <- c('cad64_2008mar26_treg_i007576j_017', 'cad64_2008mar26_treg_i007560r_009', 'cad64_2008mar26_treg_i007584s_021')

for (fcs.file in fcs.files) {

load(sprintf('np.mm.%s.fcs.2.obj', fcs.file))
load(sprintf('mm.%s.fcs.obj', fcs.file))

#mixture
plot(density(np.mm$data, bw=.01), col="gray", main=fcs.file)
hist(np.mm$data, prob=TRUE, add=T, border="gray")
sort(mm$x) -> x
lines(x, mm$lambda[[1]]*dnorm(x, mm$mu[[1]], mm$sigma[[1]])+mm$lambda[[2]]*dnorm(x, mm$mu[[2]], mm$sigma[[2]]), lty=2)
density.npEM(np.mm, component=1, scale=F)->d1
density.npEM(np.mm, component=2, scale=F)->d2
lines(d1$x, d1$y*np.mm$lambdahat[1]+d2$y*np.mm$lambdahat[2], col="black", lty=2)

#individual scaled components
plot(density(np.mm$data, bw=.01), col="gray")
hist(np.mm$data, prob=TRUE, add=T, border="gray")
for (i in 1:2) {
    lines(x, mm$lambda[i] * dnorm(x, mean = mm$mu[i], sd = mm$sigma[i]), col = 1+i)
    density.npEM(np.mm, component=i, scale=T)->d
    lines(d$x, d$y, col=i+1)
}

#posteriors
plot(density(np.mm$data, bw=.01), col="gray", ylim=c(0,1))
hist(np.mm$data, prob=TRUE, add=T, border="gray")
np.mm.density <- (d1$y*np.mm$lambdahat[1]+d2$y*np.mm$lambdahat[2])
lines(d1$x, d1$y*np.mm$lambdahat[1] / np.mm.density, col="red")
lines(d1$x, d2$y*np.mm$lambdahat[2] / np.mm.density, col="green")
mm.density <- mm$lambda[[1]]*dnorm(x, mm$mu[[1]], mm$sigma[[1]])+mm$lambda[[2]]*dnorm(x, mm$mu[[2]], mm$sigma[[2]])
lines(x, mm$lambda[[1]]*dnorm(x, mm$mu[[1]], mm$sigma[[1]])/mm.density, col="red")
lines(x, mm$lambda[[2]]*dnorm(x, mm$mu[[2]], mm$sigma[[2]])/mm.density, col="green")
abline(h=.95, col="gray")

}

###

library(flowCore)
fcs.data <- read.FCS('~nikolas/thor/FCS/cad116_2008oct09_treg_i009546a_013.nontregs.fcs')
x <- fcs.data@exprs[,8]
x <- log10(x[x>1])
fcs.data <- read.FCS('~nikolas/thor/FCS/cad116_2008oct09_treg_i009546a_013.fcs')
y <- fcs.data@exprs[,8]
y <- log10(y[y>1])
hist(y, breaks=100, prob=T, col="blue")
hist(x, breaks=100, add=T, prob=T, col="red")






