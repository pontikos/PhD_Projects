source('~nikolas/bin/load.cell.phenotypes.R')

library(flowCore)


### Repeatability of cell phenotypes per technique
plot.repeatability <- function(method, method.name, phenotypes, file.name, mfrow, ...) {
    # 2 measures of repeatability
    # Vincent's calculation
    R.squared <- function(x, y) return( 1-(mean((x-y)**2))/var(c(x,y)) )
    # Pearson's correlation. We use the cor function here.
    r.squared <- function(x, y) return( cov(x,y)/sqrt(var(x)*var(y)) )
    #pdf(file.name, width=5*mfrow[2], height=5*mfrow[1])
    #par( mfrow=mfrow, oma = c( 0, 0, 2, 0 ) )
    manual.method <- load.individual.cell.phenotypes('~nikolas/CellPhenotypes/calli.cell.phenotypes')$recalled
    if (FALSE)
    method1 <- load.individual.cell.phenotypes('~nikolas/CellPhenotypes/mm.cell.phenotypes')$recalled
    #method1 <- cd25.pos.memory.cell.phenotypes.recalled
    for (p in names(phenotypes)) {
        pdf(sprintf('%s.%s.pdf',file.name,p), pointsize=18)
        par(mgp=c(1.75,.5,0))
        #par(mfrow=c(1,1), mar=c(2, 4, 2, 2))
        par(mfrow=c(1,1))
        method[,paste(p, "day1", sep=".")]->pheno.day1
        method[,paste(p, "day2", sep=".")]->pheno.day2

        if (FALSE) {
        method1[,paste(p, "day1", sep=".")]->method1.pheno.day1
        method1[,paste(p, "day2", sep=".")]->method1.pheno.day2
        }

        manual.method[,paste(p, "day1", sep=".")]->manual.pheno.day1
        manual.method[,paste(p, "day2", sep=".")]->manual.pheno.day2

        if (FALSE)
        range(manual.pheno.day1, method1.pheno.day1, pheno.day1, manual.pheno.day2, method1.pheno.day2, pheno.day2)->pheno.range
        else
        range(manual.pheno.day1, pheno.day1, manual.pheno.day2, pheno.day2)->pheno.range

        plot(pheno.day1, pheno.day2, pch=method$pch,
             xlab=paste(phenotypes[[p]], "Day 1"),
             ylab=paste(phenotypes[[p]], "Day 2"),
             xlim=pheno.range,
             ylim=pheno.range,
             col="red", cex=.75)
             #cex=1.2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
        round(cor(pheno.day1, pheno.day2), digits=3)->r2
        round(R.squared(pheno.day1, pheno.day2), digits=3)->R2

        points(manual.pheno.day1, manual.pheno.day2, pch=manual.method$pch, col="black", cex=.75)
        round(cor(manual.pheno.day1, manual.pheno.day2), digits=3)->manual.r2
        round(R.squared(manual.pheno.day1, manual.pheno.day2), digits=3)->manual.R2

    if (FALSE) {
        points(method1.pheno.day1, method1.pheno.day2, pch=method1$pch, col="blue", cex=.75)
        round(cor(method1.pheno.day1, method1.pheno.day2), digits=3)->method1.r2
        round(R.squared(method1.pheno.day1, method1.pheno.day2), digits=3)->method1.R2
    }

        #title( as.expression(bquote(paste( (R^2 == .(R2)) (r^2 == .(r2)) (R^2 == .(manual.R2)) (r^2 == .(manual.r2)) (N==.(length(pheno.day1))) ))), line=1)
        #title(phenotypes[[p]], line=2)
        abline(a=0,b=1)
        #par(xpd=TRUE)
        #legend(pheno.range[[1]],pheno.range[[2]], c( as.expression(bquote(paste( beads (R^2 == .(R2)) (r^2 == .(r2)) ))), #as.expression(bquote(paste( memory (R^2 == .(method1.R2)) (r^2 == .(method1.r2)) ))), as.expression(bquote(paste( manual (R^2 == .(manual.R2)) (r^2 == .(manual.r2)) )))), border=FALSE, text.col=c('red', #'blue', 'black'), cex=1.1)
        dev.off()
        #par(xpd=FALSE)
    }
    #title(method.name, outer=T)
}

plot.repeatability(method=load.individual.cell.phenotypes('~nikolas/CellPhenotypes/mm.cell.phenotypes')$recalled, method.name="MM", phenotypes=list( APC.mfi="Memory CD25 MFI", APC.mef="Memory CD25 MEF", memory.freqpar="Memory %"), file.name="~/Plots/Repeatability/mm.repeatability", mfrow=c(2,2))

#plot.repeatability(method=recalled.vincent, method.name="Calliope Gating (XML Parsing)", phenotypes=list( APC.mfi="Memory CD25 MFI", APC.mef="Memory CD25 MEF", memory.freqpar="Memory %", naive.cd25pos.freqpar="CD25+ Naive %"), file.name="vincent.repeatability2.pdf", mfrow=c(2,2))
#plot.repeatability(method=recalled.calli, method.name="Calliope Gating (flowFlowJo)", phenotypes=list( APC.mfi="Memory CD25 MFI", APC.mef="Memory CD25 MEF", memory.freqpar="Memory %", naive.cd25pos.freqpar="CD25+ Naive %"), file.name="calli.repeatability.pdf", mfrow=c(2,2))
#plot.repeatability(method=recalled.auto, method.name="CD45RA Thresholding", phenotypes=list( APC.mfi="Memory CD25 MFI", APC.mef="Memory CD25 MEF", memory.freqpar="Memory %", naive.cd25pos.freqpar="CD25+ Naive %"), file.name="auto2.repeatability.pdf", mfrow=c(2,2))
plot.repeatability(method=load.individual.cell.phenotypes('~nikolas/CellPhenotypes/sp.mm.cell.phenotypes')$recalled, method.name="SP.MM", phenotypes=list( APC.mfi="Memory CD25 MFI", APC.mef="Memory CD25 MEF", memory.freqpar="Memory %"), file.name="sp.mm.repeatability.pdf", mfrow=c(2,2))
plot.repeatability(method=load.individual.cell.phenotypes('~nikolas/CellPhenotypes/mm.cell.phenotypes')$recalled, method.name="MM", phenotypes=list( APC.mfi="Memory CD25 MFI", APC.mef="Memory CD25 MEF", memory.freqpar="Memory %"), file.name="mm.repeatability.pdf", mfrow=c(2,2))
#plot.repeatability(method=recalled.manual, method.name="Manual Gates", phenotypes=list( APC.mfi="Memory CD25 MFI", APC.mef="Memory CD25 MEF", memory.freqpar="Memory %"), file.name="manual.repeatability.pdf", mfrow=c(2,2))
plot.repeatability(method=load.individual.cell.phenotypes('~nikolas/CellPhenotypes/gates.6.cell.phenotypes')$recalled, method.name="Learned Gates", phenotypes=list( APC.mfi="Memory CD25 MFI", APC.mef="Memory CD25 MEF", memory.freqpar="Memory %"), file.name="gates.6.repeatability.pdf", mfrow=c(2,2))
plot.repeatability(method=load.individual.cell.phenotypes('~nikolas/CellPhenotypes/fixed.gate.cell.phenotypes')$recalled, method.name="Fixed Gate", phenotypes=list( APC.mfi="Memory CD25 MFI", memory.freqpar="Memory %", naive.cd25pos.freqpar="CD25+ Naive %"), file.name="~/Plots/Repeatability/fixed.repeatability.pdf", mfrow=c(2,2))
#plot.repeatability(method=load.individual.cell.phenotypes('~nikolas/CellPhenotypes/calli.cell.phenotypes')$recalled, method.name="Manual", phenotypes=list( APC.mfi="Memory CD25 MFI", APC.mef="Memory CD25 MEF", memory.freqpar="Memory %", naive.cd25pos.freqpar="CD25+ Naive %"), file.name="~/Plots/Repeatability/manual.repeatability", mfrow=c(2,2))
### 
plot.repeatability.raw.gates <- function(individuals, file.name, gates.file, gates.file2, gate.name) {
    read.csv(gates.file)->gates
    read.csv(gates.file2)->gates2
    as.character(gates$fcsFile)->gates$fcsFile
    pdf(file.name)
    par(mfrow=c(2,2))
    day1 <- list()
    day2 <- list()
    as.character(individuals$fcsFile.day1)->individuals$fcsFile.day1
    as.character(individuals$fcsFile.day2)->individuals$fcsFile.day2
    for (r in 1:nrow(individuals)) {
        individuals[r,] -> d
        read.FCS(sprintf('~nikolas/dunwich/FCS.Gated/Calli/NON_T_REGS/%s.fcs', d$fcsFile.day1))@exprs[,c(5, 8)] -> d.1
        gates[gates$fcsFile==d$fcsFile.day1,gate.name]->gate.day1
        gates2[gates2$fcsFile==d$fcsFile.day1,gate.name]->gate2.day1
        day1[[r]]<-log10(d.1)
        read.FCS(sprintf('~nikolas/dunwich/FCS.Gated/Calli/NON_T_REGS/%s.fcs', d$fcsFile.day2))@exprs[,c(5, 8)] -> d.2
        gates[gates$fcsFile==d$fcsFile.day2,gate.name]->gate.day2
        gates2[gates2$fcsFile==d$fcsFile.day2,gate.name]->gate2.day2
        day2[[r]]<-log10(d.2)
        #
        plot(log10(d.1), main=paste(d$individual, d$pch, d$date.day1, dim(d.1)[[1]]), pch=".", xlim=c(0, 2), ylim=c(0, 2))
        (log10(d.1[,1]))->x
        (log10(d.1[,2]))->y
        abline(v=mean(x))
        abline(h=mean(y))
        mixtools::ellipse(mu=c(mean(x),mean(y)), sigma=cov(cbind(x,y)), alpha=.5)
        abline(h=gate.day1, col="red")
        abline(h=gate2.day1, col="blue")
        #flowMeans(log10(d.1), NumC=2)->res
        #cbind(x,y,res@Label)->r
        #points(r[,1:2], pch=".", col=r[,3])
        #
        plot(log10(d.2), main=paste(d$individual, d$pch, d$date.day2, dim(d.2)[[1]]), pch=".", xlim=c(0, 2), ylim=c(0, 2))
        (log10(d.2[,1]))->x
        (log10(d.2[,2]))->y
        abline(v=mean(x))
        abline(h=mean(y))
        mixtools::ellipse(mu=c(mean(x),mean(y)), sigma=cov(cbind(x,y)), alpha=.5)
        abline(h=gate.day2, col="red")
        abline(h=gate2.day2, col="blue")
        #flowMeans(log10(d.2), NumC=2)->res
        #cbind(x,y,res@Label)->r
        #points(r[,1:2], pch=".", col=r[,3])
    }
    dev.off()
    return(list(day1=day1, day2=day2))
}
plot.repeatability.raw.gates(load.individual.cell.phenotypes('~nikolas/CellPhenotypes/manual.cell.phenotypes')$recalled,
                             'recalled.mm.manual.raw.data.pdf',
                             #'~nikolas/CellPhenotypes/Gates/gates6.csv',
                             '~nikolas/CellPhenotypes/Gates/manual.cd45ra.cd25.csv',
                             #'~nikolas/CellPhenotypes/Gates/gates1.csv',
                             '~nikolas/CellPhenotypes/Gates/mm.cd45ra.csv',
                             'cd45ra.neg') -> res

### 
plot.raw.mef.repeatability <- function(individuals, file.name) {
    pdf(file.name)
    par(mfrow=c(2,2))
    for (r in 1:nrow(individuals)) {
        individuals[r,] -> d
        try(read.csv(sprintf('~nikolas/dunwich/FCS.Gated/Calli/NON_T_REGS/%s.mef.csv', d$fcsFile.day1))[,2]) -> d.1.apc
        if (class(d.1.apc)=='try-error') next
        (read.FCS(sprintf('~nikolas/dunwich/FCS.Gated/Calli/NON_T_REGS/%s.fcs', d$fcsFile.day1))@exprs[,8]) -> d.1.pacific
        try(read.csv(sprintf('~nikolas/dunwich/FCS.Gated/Calli/NON_T_REGS/%s.mef.csv', d$fcsFile.day2))[,2]) -> d.2.apc
        if (class(d.2.apc) == 'try-error') next
        (read.FCS(sprintf('~nikolas/dunwich/FCS.Gated/Calli/NON_T_REGS/%s.fcs', d$fcsFile.day2))@exprs[,8]) -> d.2.pacific
        plot(log10(d.1.apc), log10(d.1.pacific), main=paste(d$individual, d$pch, d$date.day1, length(d.1.apc)), pch=".", xlim=c(1,3.5), ylim=c(0,2))
        plot(log10(d.2.apc), log10(d.2.pacific), main=paste(d$individual, d$pch, d$date.day2, length(d.2.apc)), pch=".", xlim=c(1,3.5), ylim=c(0,2))
    }
    dev.off()
}
plot.raw.mef.repeatability(recalled.calli, 'recalled.calli.raw.data.mef.2.pdf')

### 
plot.raw.overlap.repeatability <- function(individuals, file.name) {
    pdf(file.name)
    par(mfrow=c(2,2))
    for (r in 1:nrow(individuals)) {
        individuals[r,] -> d
        #day1
        try(read.csv(sprintf('~nikolas/dunwich/FCS.Gated/Calli/NON_T_REGS/%s.mef.csv', d$fcsFile.day1))[,2]) -> d.mef.1
        if (class(d.mef.1)=='try-error') next
        (read.FCS(sprintf('~nikolas/dunwich/FCS.Gated/Calli/NON_T_REGS/%s.fcs', d$fcsFile.day1))@exprs[,c(5, 8)]) -> d.fi.1
        #day2
        try(read.csv(sprintf('~nikolas/dunwich/FCS.Gated/Calli/NON_T_REGS/%s.mef.csv', d$fcsFile.day2))[,2]) -> d.mef.2
        if (class(d.mef.2) == 'try-error') next
        (read.FCS(sprintf('~nikolas/dunwich/FCS.Gated/Calli/NON_T_REGS/%s.fcs', d$fcsFile.day2))@exprs[,c(5, 8)]) -> d.fi.2
        #
        log10(d.fi.1[,2])-> pacific.day1
        log10(d.mef.1)->apc.mef.day1
        log10(d.fi.1[,1])-> apc.day1
        #
        log10(d.fi.2[,2])-> pacific.day2
        log10(d.mef.2)->apc.mef.day2
        log10(d.fi.2[,1])-> apc.day2
        #day1
        plot(apc.mef.day1, pacific.day1, main=paste(d$individual, d$pch), pch=".", col='red', ylim=c(0,2))
        points(apc.mef.day2, pacific.day2, pch=".", col='blue')
        #day2
        plot(apc.day1, pacific.day1, main=paste(d$individual, d$pch), pch=".", col='red', ylim=c(0,2))
        points(apc.day2, pacific.day2, pch=".", col='blue')
    }
    dev.off()
}
plot.raw.overlap.repeatability(recalled.calli, 'recalled.calli.overlap.mef.2.pdf')

###
plot.repeatability.channel.1d <- function(individuals, file.name, cell.subset) {
    pdf(file.name)
    par(mfrow=c(2,2))
    for (r in 1:nrow(individuals)) {
        individuals[r,] -> d
        #day1
        try(read.csv(sprintf('~nikolas/dunwich/FCS.Gated/%s/%s.mef.csv', cell.subset, d$fcsFile.day1))[,1]) -> apc.mef.day1
        if (class(apc.mef.day1)=='try-error') next
        (read.FCS(sprintf('~nikolas/dunwich/FCS.Gated/%s/%s.fcs', cell.subset, d$fcsFile.day1))@exprs[,5]) -> apc.day1
        #day2
        try(read.csv(sprintf('~nikolas/dunwich/FCS.Gated/%s/%s.mef.csv', cell.subset, d$fcsFile.day2))[,1]) -> apc.mef.day2
        if (class(apc.mef.day2) == 'try-error') next
        (read.FCS(sprintf('~nikolas/dunwich/FCS.Gated/%s/%s.fcs', cell.subset, d$fcsFile.day2))@exprs[,5]) -> apc.day2
        #MEF
        min(apc.mef.day1)->min.mef.day1
        min(apc.mef.day2)->min.mef.day2
        #
        log10(apc.mef.day1[apc.mef.day1 > min.mef.day1]) -> x.mef.day1
        mean(x.mef.day1) -> mean.mef.day1
        log10(apc.mef.day2[apc.mef.day2 > min.mef.day2]) -> x.mef.day2
        mean(x.mef.day2) -> mean.mef.day2
        #
        density(x.mef.day1) -> density.mef.day1
        density(x.mef.day2) -> density.mef.day2
        plot(density.mef.day1, main=paste("MEF", d$individual, d$pch), pch=".", col='red', ylim=c(c(0,max(density.mef.day1$y, density.mef.day2$y))))
        lines(density.mef.day2, col='blue')
        abline(v=mean.mef.day1, col='red')
        abline(v=mean.mef.day2, col='blue')
        #
        ecdf(x.mef.day1) ->> ecdf.mef.day1
        ecdf(x.mef.day2) ->> ecdf.mef.day2
        sum(sapply(seq(1.5,3.5,.05), function(x) abs(ecdf.mef.day1(x)-ecdf.mef.day2(x)))) -> ecdf.mef.diff
        plot(ecdf.mef.day1, main=paste("MEF", d$individual, d$pch), pch=".", col='red', xlab=round(ecdf.mef.diff,3))
        plot(ecdf.mef.day2, col='blue', pch=".",add=T)
        abline(v=mean.mef.day1, col='red')
        abline(v=mean.mef.day2, col='blue')
        #FI
        log10(apc.day1[apc.day1 > 1]) -> x.day1
        mean(x.day1) -> mean.day1
        log10(apc.day2[apc.day2 > 1]) -> x.day2
        mean(x.day2) -> mean.day2
        #
        density(x.day1) -> density.day1
        density(x.day2) -> density.day2
        plot(density.day1, main=paste("FI", d$individual, d$pch), pch=".", col='red', ylim=c(c(0,max(density.day1$y, density.day2$y))))
        lines(density.day2, col='blue')
        abline(v=mean.day1, col='red')
        abline(v=mean.day2, col='blue')
        #
        ecdf(x.day1) ->> ecdf.day1
        ecdf(x.day2) ->> ecdf.day2
        sum(sapply(seq(0,2,.05), function(x) abs(ecdf.day1(x)-ecdf.day2(x)))) -> ecdf.diff
        plot(ecdf.day1, main=paste("FI", d$individual, d$pch), pch=".", col='red', xlab=round(ecdf.diff,3))
        plot(ecdf.day2, col='blue', pch=".", add=T)
        abline(v=mean.day1, col='red')
        abline(v=mean.day2, col='blue')
    }
    dev.off()
}
#plot.repeatability.channel.1d(recalled.calli, file.name='recalled.calli.nontreg.apc.pdf', cell.subset='Calli/NON_T_REGS')
#plot.repeatability.channel.1d(recalled.calli, file.name='recalled.calli.memory.apc.pdf', cell.subset='Calli/NON_T_REGS/MEMORY')
plot.repeatability.channel.1d(recalled.calli, file.name='recalled.calli.naive.apc.pdf', cell.subset='Calli/NON_T_REGS/NAIVE')




