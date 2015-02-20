#/home/vplagnol/export/FACS/rawData
#beeswarm package
library(gplots)
library(beeswarm)
library(xtable)
library(flowCore)




#"date"
#"individual"
#"fcsFile"              

phenotypes.1 <- c("APC.mfi", "APC.mef", "memory.freqpar", "naive.cd25pos.freqpar" )
#phenotypes.2 <- c("memory.freqpar", "memory.count", "naive.count", "naive.cd25pos.count", "naive.cd25neg.count", "naive.cd25pos.freqpar" )
#phenotypes.log.transformed <- c("memory.freqpar", "naive.cd25pos.freqpar" )
phenotypes.log.transformed <- c()

ci <- function(phenotype, d, file.name, trans=identity) {
    pdf(file.name)
    b<-matrix(nrow=0,ncol=4)
    model(trans(d[,phenotype]), d)->m
    m$coefficients[-1,c(1,2)]->b
    max(abs(b[,1]))+.5->b.max
    plotCI(b[,1], uiw=b[,2]*qt(.975,df=m$df[2]), xlim=c(1,5), ylab=phenotype, xlab="")
    axis(side = 1, at = 1:5, lab=F)
    text(axTicks(1), par("usr")[3] - .5, srt=45, adj=1, labels=rownames(b), xpd=T, cex=0.8)
    #axis(side = 2, at = seq(-b.max,b.max), labels = seq(-b.max,b.max), cex = 0.7, las=1 )
    abline(h=0, lty=2)

    plot.2 <- function(d, ...) {
    model(trans(d[,phenotype]), d)->m
    m$coefficients[-1,c(1,2)]->b
    max(abs(b[,1]))+.5->b.max
    plotCI(b[,1], uiw=b[,2]*qt(.975,df=m$df[2]), xlim=c(1,5), ylab=phenotype, xlab="", axes=F, ...)
    }
    plot.2(vincent, col="red", add=T, pch="o")
    plot.2(auto, col="blue", add=T, pch="x")

    dev.off()
}
ci("memory.freqpar", calli, trans=log, file.name="confidence.intervals.pdf")
ci("memory.freqpar", calli, file.name="confidence.intervals.pdf")

# Calculate Relative Importance for Each Predictor
library(relaimpo)
library(bootstrap)
lm(APC.mfi ~ as.numeric(snp95) + as.numeric(snp86) + as.numeric(snp56) + age + sex, data=vincent)->fit
calc.relimp(fit,type=c("lmg","last","first","pratt"), rela=TRUE)
# Bootstrap Measures of Relative Importance (1000 samples) 
boot <- boot.relimp(fit, b = 1000, type = c("lmg", "last", "first", "pratt"), rank = TRUE, diff = TRUE, rela = TRUE)
booteval.relimp(boot) # print result
pdf('booteval.pdf')
plot(booteval.relimp(boot,sort=TRUE)) # plot result
dev.off()



phenotypes <- function(d, d2) {
    for (p in phenotypes.1)
     print(xtable(cbind(model(d[,p], d), model(d2[,p], d2))[,c(1,3,2,4)], caption=sprintf("%s cell phenotype. Estimates and p-values for both flowFlowJo and Vincent's xml parsing approach for extracting gates from FlowJo workspace.",p)))
    for (p in phenotypes.log.transformed)
     print(xtable(cbind(model(log(d[,p]), d), model(log(d2[,p]), d2))[,c(1,3,2,4)], caption=sprintf("%s cell phenotype (log transformed).  Estimates and p-values for both flowFlowJo and Vincent's xml parsing approach for extracting gates from FlowJo workspace.",p)))
} 
phenotypes(calli, vincent)



#
#read.csv("~nikolas/dunwich/calli.memory.mfi")->phenotype
#na.omit(merge(INDIVIDUAL.INFO, phenotype, by=c("individual")))->d.2

plot.2.df <- function(phenotype, categorical, file.name) {
    pdf(file.name)
    lm(phenotype ~ (categorical))->m
    (summary.aov(m)[[1]])[1,"Pr(>F)"]->p.value
    summary(m)$sigma->rse
    #plot
    beeswarm(phenotype ~ categorical, main=paste("p.value", p.value, "rse", rse))->bs
    y.bar <- tapply(bs$y, FUN = mean, na.rm = TRUE, IND = (bs$x.orig))
    as.matrix(cbind(y.bar, y.bar))->y
    x.range <- tapply(bs$x, FUN = range, na.rm = TRUE, IND = (bs$x.orig))
    as.matrix(t(sapply(x.range, as.numeric)))->x
    segments(x[,1], y[,1], x[,2], y[,2], col="red")
    #abline(coefficients(m), lty=2)
    dev.off()
}
plot.2.df(d$memory.mef, d$snp95, "memory.mef.snp95.2df.pdf")


plot.1.df <- function(phenotype, d, file.name, effect, effect.names=c("snp95","snp86","snp56","Age","Male")) {
    pdf(file.name)
    model(phenotype, d)->m
    #discard intercept
    coefficients(m)[-1,]->coeff
    rownames(coeff)<-effect.names
    as.list((coeff[,"Estimate"]))->estimates
    as.list((coeff[,"Std. Error"]))->std.errors
    as.list((coeff[,"Pr(>|t|)"]))->p.values
    print(coeff)
    #(summary.aov(m)[[1]])[1,"Pr(>F)"]->p.value
    #summary(m)$sigma->rse
    #plot
    beeswarm(phenotype ~ d[,effect], main=paste("estimate", round(estimates[[effect]],2), "p.value", round(p.values[[effect]],2), "rse", round(std.errors[[effect]],2)))->bs
    y.bar <- tapply(bs$y, FUN = mean, na.rm = TRUE, IND = (bs$x.orig))
    as.matrix(cbind(y.bar, y.bar))->y
    x.range <- tapply(bs$x, FUN = range, na.rm = TRUE, IND = (bs$x.orig))
    as.matrix(t(sapply(x.range, as.numeric)))->x
    segments(x[,1], y[,1], x[,2], y[,2], col="red")
    slope <- estimates[[effect]]
    intercept <- mean(phenotype)-2*slope
    abline(intercept, slope, lty=2)
    dev.off()
}
plot.1.df(phenotype=auto$memory.freqpar, d=auto, file.name="auto.memory.freqpar.snp86.pdf", effect="snp86")
plot.1.df(phenotype=log(auto$memory.freqpar), d=auto, file.name="auto.memory.freqpar.snp86.log.pdf", effect="snp86")

plot.1.df(phenotype=auto$naive.cd25pos.freqpar, d=auto, file.name="auto.naive.cd25pos.freqpar.snp86.pdf", effect="snp86")
plot.1.df(phenotype=log(auto$naive.cd25pos.freqpar), d=auto, file.name="auto.naive.cd25pos.freqpar.snp86.log.pdf", effect="snp86")

plot.1.df(phenotype=calli$naive.cd25pos.freqpar, d=calli, file.name="calli.naive.cd25pos.freqpar.snp86.pdf", effect="snp86")
plot.1.df(phenotype=log(calli$naive.cd25pos.freqpar), d=calli, file.name="calli.naive.cd25pos.freqpar.snp86.log.pdf", effect="snp86")
#
#
plot.1.df(calli$memory.freqpar, calli, "calli.memory.freqpar.snp95.pdf")
plot.1.df(log(calli$memory.freqpar), calli, "calli.memory.freqpar.snp95.log.pdf")

plot.1.df(auto$memory.freqpar, auto$snp95, "auto.memory.freqpar.snp95.pdf")
plot.1.df(log(auto$memory.freqpar), auto$snp95, "auto.memory.freqpar.snp95.log.pdf")

#plot.1.df(d$memory.mef, d$snp95, "memory.mef.snp95.pdf")
#plot.1.df(d$memory.mfi, d$snp95, "memory.mfi.snp95.pdf")
#plot.1.df(d$memory.mfi, d$sex, "memory.mfi.sex.pdf")

plot.2.df <- function(phenotype, snp, file.name) {
    pdf(file.name)
    lm(phenotype ~ (snp))->m
    (summary.aov(m)[[1]])[1,"Pr(>F)"]->p.value
    summary(m)$sigma->rse
    #plot
    beeswarm(phenotype ~ snp, main=paste("p.value", p.value, "rse", rse))->bs
    y.bar <- tapply(bs$y, FUN = mean, na.rm = TRUE, IND = (bs$x.orig))
    as.matrix(cbind(y.bar, y.bar))->y
    x.range <- tapply(bs$x, FUN = range, na.rm = TRUE, IND = (bs$x.orig))
    as.matrix(t(sapply(x.range, as.numeric)))->x
    segments(x[,1], y[,1], x[,2], y[,2], col="red")
    abline(coefficients(m), lty=2)
    dev.off()
}
plot.2.df(d$memory.mfi, d$snp95, "memory.mfi.sex.pdf")

plot.age.effect <- function(phenotype, age, file.name) {
    pdf(file.name)
    lm(phenotype ~ age) -> m
    (summary.aov(m)[[1]])[1,"Pr(>F)"]->p.value
    summary(m)$sigma->rse
    plot(age, phenotype, main=paste("p.value", p.value, "rse", rse))
    abline(coefficients(m))
    dev.off()
}
plot.age.effect(d$memory.mfi, d$age, "memory.mfi.age.pdf")

summary(lm(phenotype ~ (snp95), data=d))

#summary(lm(phenotype ~ as.numeric(snp86), data=d))
#summary(lm(phenotype ~ as.numeric(snp56), data=d))


summary(lm(phenotype ~ as.numeric(snp56)*as.numeric(snp86)*as.numeric(snp95)*age*sex,data=genotype.phenotype))
summary(m)

pdf("plot.pdf")
beeswarm(memory.mef ~ rs12722495, data=genotype.phenotype)
boxplot(memory.mef ~ rs12722495, data=genotype.phenotype, add=T,names = c("","", ""), col="#0000ff22") 

   lines(x=x.range
dev.off()

#'Calli_foxp3_bright_CBR200.csv'
#reg
#reg.mef
#reg.percentage

par(mfrow=c(1,1))
scatter.plot(memory.mef,snp86, pch="+")
scatter.plot(memory.mef,snp56, pch="+")

plot.1.df <- function(phenotype, genotype) {
    par(mfrow=c(1,1))
    summary(lm(phenotype ~ genotype))->m
    print(summary(m))
    lm(phenotype ~ as.numeric(genotype))->m.2
    summary.aov(m.2)[[1]][[5]]->p.2
    print(summary(m.2))
    scatter.plot(phenotype,genotype, pch="+",main=p.2)
    abline(coefficients(m.2), lty=2)
}

plot.1.df(memory.mef, snp95)
plot.1.df(memory.percentage, snp95)

plot.1.df(memory.mef, snp86)
plot.1.df(memory.mef, snp56)

summary(lm(memory.mef ~ snp95))
lm(memory.mef ~ as.numeric(snp95))->m
scatter.plot(memory.mef,snp95, pch="+")
abline(coefficients(m), lty=2)


summary(lm(memory.mef ~ snp86))
summary(lm(memory.mef ~ snp56))


par(mfrow=c(3,1))
scatter.plot(memory.percentage,snp95, pch="+")
scatter.plot(memory.percentage,snp86, pch="+")
scatter.plot(memory.percentage,snp56, pch="+")

par(mfrow=c(3,1))
scatter.plot(naive.percentage,snp95, pch="+")
scatter.plot(naive.percentage,snp86, pch="+")
scatter.plot(naive.percentage,snp56, pch="+")


summary(lm(memory.mef ~ snp95))
summary(lm(memory.mef ~ snp86))
summary(lm(memory.mef ~ snp56))

