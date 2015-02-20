library(iterators)
library(plyr)
library(flowCore)
library(lubridate)

REPEATS <- structure(list(individual = c("CB00165D", "CB00366X", "CB00396E", 
"CB00406Q", "CB01484M", "CB01494Y", "CB01495Z", "CB01498C", "CB01503H", 
"CB01504J"), pch = c("a", "b", "c", "d", "e", "f", "g", "h", 
"i", "j"), day1 = c("2012-11-29", "2012-11-07", "2012-09-25", 
"2012-10-16", "2012-09-25", "2012-10-09", "2012-10-09", "2012-10-16", 
"2012-11-07", "2012-11-07"), day2 = c("2013-03-07", "2013-03-27", 
"2013-03-11", "2013-01-22", "2013-03-11", "2013-01-29", "2013-01-29", 
"2013-01-22", "2013-03-07", "2013-03-27"), col = c("#0066FFFF", 
"#FF9900FF", "#00FFFFFF", "#FF0099FF", "#33FF00FF", "#CCFF00FF", 
"#CC00FFFF", "#3300FFFF", "#00FF66FF", "#FF0000FF"), day.diff = structure(c(98, 
140, 167, 98, 167, 112, 112, 98, 120, 140), units = "days", class = "difftime"), 
t1d = c('control', 'case', 'control', 'control', 'case', 'case', 'control', 'case', 'case', 'control')),
.Names = c("individual", "pch", "day1", "day2", "col", "day.diff", "t1d"), row.names = c("CB00165D", 
"CB00366X", "CB00396E", "CB00406Q", "CB01484M", "CB01494Y", "CB01495Z", 
"CB01498C", "CB01503H", "CB01504J"), class = "data.frame") 
REPEATS <- rbind(cbind(REPEATS,date=REPEATS$day1),cbind(REPEATS,date=REPEATS$day2))
DOSES <- c( '0U', '0.1U', '10U', '1000U') 
CLR.CELL.TYPES <- c("Lymphocytes", "Single Cells", "CD4", "Memory", "Memory Eff", "Memory Treg", "Naive", "Naive Eff", "Naive Treg")
CELL.TYPES <- c('Memory Eff', 'Memory Treg', 'Naive Eff', 'Naive Treg') 
blues4 <- blues9[5:9]



BASE.DIR <- '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/Lymphocytes/RData' 
#header
cat('individual', 'norm', 'date.day1', 'date.day2', 'fi.auc.day1', 'pct.auc.day1', 'fi.auc.day2', 'pct.auc.day2', sep=',')
cat('\n')
inverse = function (f, lower = -100, upper = 100) {
   function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
}
getDate <- function(x) format(as.Date(dmy(x@description$`$DATE`, quiet=T)), '%Y-%m-%d') 
f <- function(x, y, qmin=0, qmax=1) (quantile(x,qmax)-quantile(x,qmin)) * (y-quantile(y,qmin))/(quantile(y,qmax)-quantile(y,qmin)) + quantile(x,qmin) 
g <- curv2Filter("FSC-A", "SSC-A", filterId="test filter") 
#head(d[d$individual=='CB01494Y',], n=101)->d2
print(pch <- data.frame(cbind(name=rep.individuals, pch=letters[1:length(rep.individuals)]),stringsAsFactors=F))
d -> d2 
pdf('~nikolas/IL2/Plots/Rplot.pdf')
r <- range(d2[,2],d2[,3])
plot(d2[,2], d2[,3],pch=pch[pch$name==d2$individual,'pch'], xlab='pSTAT5 MFI Day 1', ylab='pSTAT5 MFI Day 2', xlim=r, ylim=r)
abline(b=1,a=0)
#co <- coef(line(d2[,2],d2[,3]))
#abline(co)
d3 <- as.matrix(cbind(1,d2[,2:3]))
#points(d3[,c(1,2)]%*%co, d3[,3], pch=20)
dev.off() 
#
f1 <- rep.fcs$CB01494Y$day1$U1000
f2 <- rep.fcs$CB01494Y$day2$U1000
#
f1 <- read.FCS('~/dunwich/spade.output/units/T1D_TC8_pSTAT5_Treg_I021767J_CB01494Y_0U/T1D_TC8_pSTAT5_Treg_I021767J_CB01494Y_0U.fcs.downsample.fcs')
f2 <- read.FCS('~/dunwich/spade.output/units/T1D_TC8_rep_pSTAT5_Treg_CB01494Y_I023360Q_0U/T1D_TC8_rep_pSTAT5_Treg_CB01494Y_I023360Q_0U.fcs.downsample.fcs')
#
f1 <- read.FCS('~/dunwich/spade.output/units/T1D_TC8_pSTAT5_Treg_I021767J_CB01494Y_1000U/T1D_TC8_pSTAT5_Treg_I021767J_CB01494Y_1000U.fcs.downsample.fcs')
f2 <- read.FCS('~/dunwich/spade.output/units/T1D_TC8_rep_pSTAT5_Treg_CB01494Y_I023360Q_1000U/T1D_TC8_rep_pSTAT5_Treg_CB01494Y_I023360Q_1000U.fcs.downsample.fcs')
d <- data.frame()
for (x in colnames(f1@exprs)) {
    qq<-qqplot(lgcl(f1@exprs[,x]),lgcl(f2@exprs[,x]),plot.it=F)
    l <- line(qq$y, qq$x)
    a <- l$coefficients[1]
    b <- l$coefficients[2]
    d <- rbind(d, data.frame(param=x, a=a, b=b))
}
d1000.down <- d 
source("pstat-rep_functions.R")
trans <- logicleTransform() 
for (individual in rep.individuals) {
  cat(individual,"\t")
  dir.create(individual,showWarnings=FALSE)
  make.plots(individual)
} 
for (individual in rep.individuals) {
  cat(individual,"\t")
  dir.create(individual,showWarnings=FALSE)
  make.plots(individual, cell.type='lymphocytes')
}


individual <- "CB01494Y"
##outdir <- sprintf('/home/nikolas/Tony/pstat5-repeatability/pstat5-norm/%s/', individual)
outdir <- individual
dir.create(outdir)

##pdf(sprintf('/home/nikolas/Tony/pstat5-repeatability/pstat5-norm/%s/norm-none.pdf', individual)) 
myp(day1.1, day1.2, day=1, new=TRUE)
myp(day2.1, day2.2, day=2) 
myp(day1.1, day1.2, day=1, new=TRUE, adjust=TRUE)
myp(day2.1, day2.2, day=2, adjust=TRUE) 
myp(day1.1, day1.2, day=1, q=0.01, new=TRUE)
myp(day2.1, day2.2, day=2, q=0.01) 
baseline.conc <- trans(day2.1@exprs[,7])
baseline.conc.ecdf <- ecdf(baseline.conc)
high.conc <- trans(day2.2@exprs[,7])
high.conc.ecdf <- ecdf(high.conc) 
myp(baseline.conc.ecdf, high.conc.ecdf, day=2)
#dev.off() 
qmin <- 0.01
qmax <- 1-qmin 
pdf(sprintf('/home/nikolas/Tony/pstat5-repeatability/pstat5-norm/%s/norm-%s-%s.pdf', individual, as.character(qmin), as.character(qmax))) 
baseline.conc <- trans(day1.1@exprs[,7])
baseline.conc.ecdf <- ecdf(baseline.conc)
high.conc <- trans(day1.2@exprs[,7])
high.conc <- f(baseline.conc, high.conc, qmin, qmax) 
high.conc.ecdf <- ecdf(high.conc) 
inv.baseline.conc <- function(x) sapply(x, function(x) inverse(baseline.conc.ecdf, -1, 5)(x)$root)
inv.high.conc <- function(x) sapply(x, function(x) inverse(high.conc.ecdf, -1, 5)(x)$root)
inv.ecdf.diff <- function(x) return(inv.high.conc(x)-inv.baseline.conc(x))
x <- seq(0,1,.01) 
plot(baseline.conc.ecdf, xlab='pstat5', main=paste(individual, 'baseline/1000U day1 vs day2 qnorm:', qmin, '-', qmax), col='white')
lines(baseline.conc.ecdf, col='lightblue')
lines(high.conc.ecdf, col='darkblue')
lines(inv.ecdf.diff(x), x, col='blue', lty=2)
inv.ecdf.diff.area <- abs(sum(inv.ecdf.diff(x)))
fi.auc.day1 <- round(inv.ecdf.diff.area,digits=2)
text(0, 1, fi.auc.day1, col='blue') 
x <- seq(min(baseline.conc), max(baseline.conc), .1)
ecdf.diff <- function(x) return(-high.conc.ecdf(x)+baseline.conc.ecdf(x))
lines(x, (ecdf.diff(x)), lty=2, col='blue')
ecdf.diff.area <- abs(sum(ecdf.diff(x)))
pct.auc.day1 <- round(ecdf.diff.area,digits=2)
text(4, 0, pct.auc.day1, col='blue') 
baseline.conc <- trans(day2.1@exprs[,7])
baseline.conc.ecdf <- ecdf(baseline.conc)
high.conc <- trans(day2.2@exprs[,7])
high.conc <- f(baseline.conc, high.conc, qmin, qmax) 
high.conc.ecdf <- ecdf(high.conc) 
inv.baseline.conc <- function(x) sapply(x, function(x) inverse(baseline.conc.ecdf, -1, 5)(x)$root)
inv.high.conc <- function(x) sapply(x, function(x) inverse(high.conc.ecdf, -1, 5)(x)$root)
inv.ecdf.diff <- function(x) return(inv.high.conc(x)-inv.baseline.conc(x))
x <- seq(0,1,.01)
lines(baseline.conc.ecdf, col='lightgreen')
lines(high.conc.ecdf, col='darkgreen')
lines(inv.ecdf.diff(x), x, col='green', lty=2)
inv.ecdf.diff.area <- abs(sum(inv.ecdf.diff(x)))
fi.auc.day2 <- round(inv.ecdf.diff.area,digits=2)
text(median(inv.ecdf.diff(x)), median(x), fi.auc.day2, col='green') 
x <- seq(min(baseline.conc), max(baseline.conc), .1)
ecdf.diff <- function(x) return(-high.conc.ecdf(x)+baseline.conc.ecdf(x))
lines(x, (ecdf.diff(x)), lty=2, col='green')
ecdf.diff.area <- abs(sum(ecdf.diff(x)))
pct.auc.day2 <- round(ecdf.diff.area,digits=2)
text(median(x), median(ecdf.diff(x)), pct.auc.day2, col='green') 
dev.off() 
cat(individual, paste(as.character(qmin),as.character(qmax),sep='-'), getDate(day1.1), getDate(day2.1), fi.auc.day1, pct.auc.day1, fi.auc.day2, pct.auc.day2, sep=',')
cat('\n')




# pSTAT5 MFI
load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/','pstat5.mfi.RData'))
raw.pstat5mfi <- pstat5.mfi 
baseline.raw.pstat5mfi <- cbind(pstat5.mfi[,c('individual','date','cell.type','0U')],pstat5.mfi[,c('01U','10U','1000U')]-pstat5.mfi[,'0U'])
pdf('~nikolas/Thesis/figures/pstat5-mfi-cellsubsets-repeatability.pdf')
par(mfrow=c(2,2))
# Memory Eff
dim(raw.repeats <- merge(raw.pstat5mfi[which(raw.pstat5mfi$cell.type=="Memory Eff"),], REPEATS))
raw.pstat5 <- cbind(raw.repeats[c(TRUE,FALSE),'10U'],raw.repeats[c(FALSE,TRUE),'10U'])
dim(baseline.raw.repeats <- merge(baseline.raw.pstat5mfi[which(baseline.raw.pstat5mfi$cell.type=="Memory Eff"),], REPEATS))
baseline.raw.pstat5 <- cbind(baseline.raw.repeats[c(TRUE,FALSE),'10U'],baseline.raw.repeats[c(FALSE,TRUE),'10U'])
xlim <- range(raw.pstat5,baseline.raw.pstat5)
plot(raw.pstat5, xlim=xlim, ylim=xlim, pch=raw.repeats[c(TRUE,FALSE),'pch'], xlab='day 1: pSTAT5 MFI', ylab='day 2: pSTAT5 MFI', col='black', main='Memory Eff 10U')
points(baseline.raw.pstat5, pch=baseline.raw.repeats[c(TRUE,FALSE),'pch'], col='red')
rp <- vector('expression',2)
rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(raw.pstat5)[1,2]**2,3))**2,3))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(baseline.raw.pstat5)[1,2]**2,3))**2,3))[2]
legend('topleft', legend=rp, text.col=c('black','red'), bty='n')
abline(b=1, a=0)
# Memory Treg
dim(raw.repeats <- merge(raw.pstat5mfi[which(raw.pstat5mfi$cell.type=="Memory Treg"),], REPEATS))
raw.pstat5 <- cbind(raw.repeats[c(TRUE,FALSE),'01U'],raw.repeats[c(FALSE,TRUE),'01U'])
dim(baseline.raw.repeats <- merge(baseline.raw.pstat5mfi[which(baseline.raw.pstat5mfi$cell.type=="Memory Treg"),], REPEATS))
baseline.raw.pstat5 <- cbind(baseline.raw.repeats[c(TRUE,FALSE),'01U'],baseline.raw.repeats[c(FALSE,TRUE),'01U'])
xlim <- range(raw.pstat5,baseline.raw.pstat5)
plot(raw.pstat5, xlim=xlim, ylim=xlim, pch=raw.repeats[c(TRUE,FALSE),'pch'], xlab='day 1: pSTAT5 MFI', ylab='day 2: pSTAT5 MFI', col='black', main='Memory Treg 0.1U')
points(baseline.raw.pstat5, pch=baseline.raw.repeats[c(TRUE,FALSE),'pch'], col='red')
rp <- vector('expression',2)
rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(raw.pstat5)[1,2]**2,3))**2,3))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(baseline.raw.pstat5)[1,2]**2,3))**2,3))[2]
legend('topleft', legend=rp, text.col=c('black','red'), bty='n')
abline(b=1, a=0)
# Naive Eff
dim(raw.repeats <- merge(raw.pstat5mfi[which(raw.pstat5mfi$cell.type=="Naive Eff"),], REPEATS))
raw.pstat5 <- cbind(raw.repeats[c(TRUE,FALSE),'1000U'],raw.repeats[c(FALSE,TRUE),'1000U'])
dim(baseline.raw.repeats <- merge(baseline.raw.pstat5mfi[which(baseline.raw.pstat5mfi$cell.type=="Naive Eff"),], REPEATS))
baseline.raw.pstat5 <- cbind(baseline.raw.repeats[c(TRUE,FALSE),'1000U'],baseline.raw.repeats[c(FALSE,TRUE),'1000U'])
xlim <- range(raw.pstat5,baseline.raw.pstat5)
plot(raw.pstat5, xlim=xlim, ylim=xlim, pch=raw.repeats[c(TRUE,FALSE),'pch'], xlab='day 1: pSTAT5 MFI', ylab='day 2: pSTAT5 MFI', col='black', main='Naive Eff 1000U')
points(baseline.raw.pstat5, pch=baseline.raw.repeats[c(TRUE,FALSE),'pch'], col='red')
rp <- vector('expression',2)
rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(raw.pstat5)[1,2]**2,3))**2,3))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(baseline.raw.pstat5)[1,2]**2,3))**2,3))[2]
legend('topleft', legend=rp, text.col=c('black','red'), bty='n')
abline(b=1, a=0)
# Naive Treg
dim(raw.repeats <- merge(raw.pstat5mfi[which(raw.pstat5mfi$cell.type=="Naive Treg"),], REPEATS))
raw.pstat5 <- cbind(raw.repeats[c(TRUE,FALSE),'01U'],raw.repeats[c(FALSE,TRUE),'01U'])
dim(baseline.raw.repeats <- merge(baseline.raw.pstat5mfi[which(baseline.raw.pstat5mfi$cell.type=="Naive Treg"),], REPEATS))
baseline.raw.pstat5 <- cbind(baseline.raw.repeats[c(TRUE,FALSE),'01U'],baseline.raw.repeats[c(FALSE,TRUE),'01U'])
xlim <- range(raw.pstat5,baseline.raw.pstat5)
plot(raw.pstat5, xlim=xlim, ylim=xlim, pch=raw.repeats[c(TRUE,FALSE),'pch'], xlab='day 1: pSTAT5 MFI', ylab='day 2: pSTAT5 MFI', col='black', main='Naive Treg 0.1U')
points(baseline.raw.pstat5, pch=baseline.raw.repeats[c(TRUE,FALSE),'pch'], col='red')
rp <- vector('expression',2)
rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(raw.pstat5)[1,2]**2,3))**2,3))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(baseline.raw.pstat5)[1,2]**2,3))**2,3))[2]
legend('topleft', legend=rp, text.col=c('black','red'), bty='n')
abline(b=1, a=0)
dev.off()



# nn pSTAT5 MFI
load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/','nn.pstat5.mfi.RData'))
nn.pstat5mfi <- pstat5.mfi[,c('individual','date','cell.type',grep('^X',colnames(pstat5.mfi),value=TRUE))]
colnames(nn.pstat5mfi) <- c('individual','date','cell.type',paste('PSTAT5',1:4,sep='.'))
baseline.nn.pstat5mfi <- pstat5.mfi[,c('individual','date','cell.type',grep('^base.*U',colnames(pstat5.mfi),value=TRUE))]
colnames(baseline.nn.pstat5mfi) <- c('individual','date','cell.type',paste('diff','PSTAT5',1:4,sep='.'))
#
pdf('~nikolas/Thesis/figures/nn-pstat5-mfi-cellsubsets-repeatability.pdf')
par(mfrow=c(2,2))
# Memory Eff
dim(nn.repeats <- merge(nn.pstat5mfi[which(nn.pstat5mfi$cell.type=="Memory Eff"),], REPEATS))
nn.pstat5 <- cbind(nn.repeats[c(TRUE,FALSE),'PSTAT5.3'],nn.repeats[c(FALSE,TRUE),'PSTAT5.3'])
dim(baseline.nn.repeats <- merge(baseline.nn.pstat5mfi[which(baseline.nn.pstat5mfi$cell.type=="Memory Eff"),], REPEATS))
baseline.nn.pstat5 <- cbind(baseline.nn.repeats[c(TRUE,FALSE),'diff.PSTAT5.3'],baseline.nn.repeats[c(FALSE,TRUE),'diff.PSTAT5.3'])
xlim <- range(nn.pstat5,baseline.nn.pstat5)
plot(nn.pstat5, xlim=xlim, ylim=c(0,1), pch=nn.repeats[c(TRUE,FALSE),'pch'], xlab='day 1: pSTAT5 MFI', ylab='day 2: pSTAT5 MFI', col='black', main='Memory Eff 10U')
points(baseline.nn.pstat5, pch=baseline.nn.repeats[c(TRUE,FALSE),'pch'], col='red')
rp <- vector('expression',2)
rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(nn.pstat5)[1,2]**2,3))**2,3))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(baseline.nn.pstat5)[1,2]**2,3))**2,3))[2]
legend('topleft', legend=rp, text.col=c('black','red'), bty='n')
abline(b=1, a=0)
# Memory Treg
dim(nn.repeats <- merge(nn.pstat5mfi[which(nn.pstat5mfi$cell.type=="Memory Treg"),], REPEATS))
nn.pstat5 <- cbind(nn.repeats[c(TRUE,FALSE),'PSTAT5.2'],nn.repeats[c(FALSE,TRUE),'PSTAT5.2'])
dim(baseline.nn.repeats <- merge(baseline.nn.pstat5mfi[which(baseline.nn.pstat5mfi$cell.type=="Memory Treg"),], REPEATS))
baseline.nn.pstat5 <- cbind(baseline.nn.repeats[c(TRUE,FALSE),'diff.PSTAT5.2'],baseline.nn.repeats[c(FALSE,TRUE),'diff.PSTAT5.2'])
xlim <- range(nn.pstat5,baseline.nn.pstat5)
plot(nn.pstat5, xlim=xlim, ylim=c(0,1), pch=nn.repeats[c(TRUE,FALSE),'pch'], xlab='day 1: pSTAT5 MFI', ylab='day 2: pSTAT5 MFI', col='black', main='Memory Treg 0.1U')
points(baseline.nn.pstat5, pch=baseline.nn.repeats[c(TRUE,FALSE),'pch'], col='red')
rp <- vector('expression',2)
rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(nn.pstat5)[1,2]**2,3))**2,3))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(baseline.nn.pstat5)[1,2]**2,3))**2,3))[2]
legend('topleft', legend=rp, text.col=c('black','red'), bty='n')
abline(b=1, a=0)
# Naive Eff
dim(nn.repeats <- merge(nn.pstat5mfi[which(nn.pstat5mfi$cell.type=="Naive Eff"),], REPEATS))
nn.pstat5 <- cbind(nn.repeats[c(TRUE,FALSE),'PSTAT5.4'],nn.repeats[c(FALSE,TRUE),'PSTAT5.4'])
dim(baseline.nn.repeats <- merge(baseline.nn.pstat5mfi[which(baseline.nn.pstat5mfi$cell.type=="Naive Eff"),], REPEATS))
baseline.nn.pstat5 <- cbind(baseline.nn.repeats[c(TRUE,FALSE),'diff.PSTAT5.4'],baseline.nn.repeats[c(FALSE,TRUE),'diff.PSTAT5.4'])
xlim <- range(nn.pstat5,baseline.nn.pstat5)
plot(nn.pstat5, xlim=xlim, ylim=c(0,1), pch=nn.repeats[c(TRUE,FALSE),'pch'], xlab='day 1: pSTAT5 MFI', ylab='day 2: pSTAT5 MFI', col='black', main='Naive Eff 1000U')
points(baseline.nn.pstat5, pch=baseline.nn.repeats[c(TRUE,FALSE),'pch'], col='red')
rp <- vector('expression',2)
rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(nn.pstat5)[1,2]**2,3))**2,3))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(baseline.nn.pstat5)[1,2]**2,3))**2,3))[2]
legend('topleft', legend=rp, text.col=c('black','red'), bty='n')
abline(b=1, a=0)
# Naive Treg
dim(nn.repeats <- merge(nn.pstat5mfi[which(nn.pstat5mfi$cell.type=="Naive Treg"),], REPEATS))
nn.pstat5 <- cbind(nn.repeats[c(TRUE,FALSE),'PSTAT5.2'],nn.repeats[c(FALSE,TRUE),'PSTAT5.2'])
dim(baseline.nn.repeats <- merge(baseline.nn.pstat5mfi[which(baseline.nn.pstat5mfi$cell.type=="Naive Treg"),], REPEATS))
baseline.nn.pstat5 <- cbind(baseline.nn.repeats[c(TRUE,FALSE),'diff.PSTAT5.2'],baseline.nn.repeats[c(FALSE,TRUE),'diff.PSTAT5.2'])
xlim <- range(nn.pstat5,baseline.nn.pstat5)
plot(nn.pstat5, xlim=xlim, ylim=c(0,1), pch=nn.repeats[c(TRUE,FALSE),'pch'], xlab='day 1: pSTAT5 MFI', ylab='day 2: pSTAT5 MFI', col='black', main='Naive Treg 0.1U')
points(baseline.nn.pstat5, pch=baseline.nn.repeats[c(TRUE,FALSE),'pch'], col='red')
rp <- vector('expression',2)
rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(nn.pstat5)[1,2]**2,3))**2,3))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(baseline.nn.pstat5)[1,2]**2,3))**2,3))[2]
legend('topleft', legend=rp, text.col=c('black','red'), bty='n')
abline(b=1, a=0)
dev.off()


# % PSTAT5+
print(load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/','pstat5.pos.RData')))
raw.pstat5pos <- pstat5.pos
print(load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/','nn.pstat5.pos.RData')))
nn.pstat5pos <- pstat5.pos
print(load(file.path('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/','nn.base.pstat5.pos.RData')))
head(baseline.nn.pstat5pos <- pstat5.pos)
#
pdf('~nikolas/Thesis/figures/pstat5-pos-cellsubsets-repeatability.pdf')
par(mfrow=c(2,2))
# Memory Eff
dim(raw.repeats <- merge(raw.pstat5pos[which(raw.pstat5pos$cell.type=="Memory Eff"),], REPEATS))
raw.pstat5 <- cbind(raw.repeats[c(TRUE,FALSE),'10U'],raw.repeats[c(FALSE,TRUE),'10U'])
dim(nn.repeats <- merge(nn.pstat5pos[which(nn.pstat5pos$cell.type=="Memory Eff"),], REPEATS))
nn.pstat5 <- cbind(nn.repeats[c(TRUE,FALSE),'10U'],nn.repeats[c(FALSE,TRUE),'10U'])
dim(baseline.nn.repeats <- merge(baseline.nn.pstat5pos[which(baseline.nn.pstat5pos$cell.type=="Memory Eff"),], REPEATS))
baseline.nn.pstat5 <- cbind(baseline.nn.repeats[c(TRUE,FALSE),'10U'],baseline.nn.repeats[c(FALSE,TRUE),'10U'])
#xlim <- range(raw.pstat5, nn.pstat5,baseline.nn.pstat5)
xlim <- c(0,100)
plot(raw.pstat5, xlim=xlim, ylim=c(0,1), pch=raw.repeats[c(TRUE,FALSE),'pch'], xlab='day 1: % pSTAT5+', ylab='day 2: % pSTAT5+', col='black', main='Memory Eff 10U')
points(nn.pstat5, pch=nn.repeats[c(TRUE,FALSE),'pch'], col='red')
points(baseline.nn.pstat5, pch=baseline.nn.repeats[c(TRUE,FALSE),'pch'], col='blue')
rp <- vector('expression',3)
rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(raw.pstat5)[1,2]**2,3))**2,3))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(nn.pstat5)[1,2]**2,3))**2,3))[2]
rp[3] <- substitute(expression(r^2 == r2),list(r2=(round(cor(baseline.nn.pstat5)[1,2]**2,3))**2,3))[2]
legend('topleft', legend=rp, text.col=c('black','red','blue'), bty='n')
abline(b=1, a=0)
# Memory Treg
dim(raw.repeats <- merge(raw.pstat5pos[which(raw.pstat5pos$cell.type=="Memory Treg"),], REPEATS))
raw.pstat5 <- cbind(raw.repeats[c(TRUE,FALSE),'01U'],raw.repeats[c(FALSE,TRUE),'01U'])
dim(nn.repeats <- merge(nn.pstat5pos[which(nn.pstat5pos$cell.type=="Memory Treg"),], REPEATS))
nn.pstat5 <- cbind(nn.repeats[c(TRUE,FALSE),'01U'],nn.repeats[c(FALSE,TRUE),'01U'])
dim(baseline.nn.repeats <- merge(baseline.nn.pstat5pos[which(baseline.nn.pstat5pos$cell.type=="Memory Treg"),], REPEATS))
baseline.nn.pstat5 <- cbind(baseline.nn.repeats[c(TRUE,FALSE),'01U'],baseline.nn.repeats[c(FALSE,TRUE),'01U'])
#xlim <- range(nn.pstat5,baseline.nn.pstat5)
plot(raw.pstat5, xlim=xlim, ylim=c(0,1), pch=raw.repeats[c(TRUE,FALSE),'pch'], xlab='day 1: % pSTAT5+', ylab='day 2: % pSTAT5+', col='black', main='Memory Treg 0.1U')
points(nn.pstat5, pch=nn.repeats[c(TRUE,FALSE),'pch'], col='red')
points(baseline.nn.pstat5, pch=baseline.nn.repeats[c(TRUE,FALSE),'pch'], col='blue')
rp <- vector('expression',3)
rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(raw.pstat5)[1,2]**2,3))**2,3))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(nn.pstat5)[1,2]**2,3))**2,3))[2]
rp[3] <- substitute(expression(r^2 == r2),list(r2=(round(cor(baseline.nn.pstat5)[1,2]**2,3))**2,3))[2]
legend('topleft', legend=rp, text.col=c('black','red','blue'), bty='n')
abline(b=1, a=0)
# Naive Eff
dim(raw.repeats <- merge(raw.pstat5pos[which(raw.pstat5pos$cell.type=="Naive Eff"),], REPEATS))
raw.pstat5 <- cbind(raw.repeats[c(TRUE,FALSE),'1000U'],raw.repeats[c(FALSE,TRUE),'1000U'])
dim(nn.repeats <- merge(nn.pstat5pos[which(nn.pstat5pos$cell.type=="Naive Eff"),], REPEATS))
nn.pstat5 <- cbind(nn.repeats[c(TRUE,FALSE),'1000U'],nn.repeats[c(FALSE,TRUE),'1000U'])
dim(baseline.nn.repeats <- merge(baseline.nn.pstat5pos[which(baseline.nn.pstat5pos$cell.type=="Naive Eff"),], REPEATS))
baseline.nn.pstat5 <- cbind(baseline.nn.repeats[c(TRUE,FALSE),'1000U'],baseline.nn.repeats[c(FALSE,TRUE),'1000U'])
#xlim <- range(nn.pstat5,baseline.nn.pstat5)
plot(raw.pstat5, xlim=xlim, ylim=c(0,1), pch=raw.repeats[c(TRUE,FALSE),'pch'], xlab='day 1: % pSTAT5+', ylab='day 2: % pSTAT5+', col='black', main='Naive Eff 1000U')
points(nn.pstat5, pch=nn.repeats[c(TRUE,FALSE),'pch'], col='red')
points(baseline.nn.pstat5, pch=baseline.nn.repeats[c(TRUE,FALSE),'pch'], col='blue')
rp <- vector('expression',3)
rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(raw.pstat5)[1,2]**2,3))**2,3))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(nn.pstat5)[1,2]**2,3))**2,3))[2]
rp[3] <- substitute(expression(r^2 == r2),list(r2=(round(cor(baseline.nn.pstat5)[1,2]**2,3))**2,3))[2]
legend('topleft', legend=rp, text.col=c('black','red','blue'), bty='n')
abline(b=1, a=0)
# Naive Treg
dim(raw.repeats <- merge(raw.pstat5pos[which(raw.pstat5pos$cell.type=="Naive Treg"),], REPEATS))
raw.pstat5 <- cbind(raw.repeats[c(TRUE,FALSE),'01U'],raw.repeats[c(FALSE,TRUE),'01U'])
dim(nn.repeats <- merge(nn.pstat5pos[which(nn.pstat5pos$cell.type=="Naive Treg"),], REPEATS))
nn.pstat5 <- cbind(nn.repeats[c(TRUE,FALSE),'01U'],nn.repeats[c(FALSE,TRUE),'01U'])
dim(baseline.nn.repeats <- merge(baseline.nn.pstat5pos[which(baseline.nn.pstat5pos$cell.type=="Naive Treg"),], REPEATS))
baseline.nn.pstat5 <- cbind(baseline.nn.repeats[c(TRUE,FALSE),'01U'],baseline.nn.repeats[c(FALSE,TRUE),'01U'])
#xlim <- range(nn.pstat5,baseline.nn.pstat5)
plot(raw.pstat5, xlim=xlim, ylim=c(0,1), pch=raw.repeats[c(TRUE,FALSE),'pch'], xlab='day 1: % pSTAT5+', ylab='day 2: % pSTAT5+', col='black', main='Naive Treg 0.1U')
points(nn.pstat5, pch=nn.repeats[c(TRUE,FALSE),'pch'], col='red')
points(baseline.nn.pstat5, pch=baseline.nn.repeats[c(TRUE,FALSE),'pch'], col='blue')
rp <- vector('expression',2)
rp[1] <- substitute(expression(r^2 == r2),list(r2=(round(cor(raw.pstat5)[1,2]**2,3))**2,3))[2]
rp[2] <- substitute(expression(r^2 == r2),list(r2=(round(cor(nn.pstat5)[1,2]**2,3))**2,3))[2]
rp[3] <- substitute(expression(r^2 == r2),list(r2=(round(cor(baseline.nn.pstat5)[1,2]**2,3))**2,3))[2]
legend('topleft', legend=rp, text.col=c('black','red','blue'), bty='n')
abline(b=1, a=0)
dev.off()




fun.plot <- function(d,main) {
    #d[is.na(d)] <- 0
    s <- 0
    if (sum(is.na(d[,'0U']))>1) {
        d <- d[,-which('0U'==colnames(d))]
        s <- 1
    }
    #plot(NULL, xlim=c(0,3), ylim=range(d[,grep('U',colnames(d))]), xaxt='n', xlab='dose', ylab=expression(r^2), main=main)
    plot(NULL, xlim=c(0,3), ylim=range(d[,grep('U',colnames(d))]), xaxt='n', xlab='dose', ylab='RMSD', main=main)
    title(nextElem(figure.labels), adj=0)
    axis(1, at=0:3, labels=DOSES)
    sapply( 1:nrow(d), function(i) lines(s:3, d[i,grep('U',colnames(d))], col=i, lwd=2) )
    legend('topleft', CELL.TYPES, text.col=1:4, bty='n')
}

# correlation
fun <- function(x) cor(x[c(TRUE,FALSE)],x[c(FALSE,TRUE)])**2
#100 - nrmsd
fun <- function(x) 100-100*sqrt(mean((x[c(TRUE,FALSE)]-x[c(FALSE,TRUE)])**2))/(max(x)-min(x))
#RMSD
fun <- function(x) log10(sqrt(mean((x[c(TRUE,FALSE)]-x[c(FALSE,TRUE)])**2)))
#Rc
fun <- function(x) {
    x1 <- x[c(TRUE,FALSE)]
    x2 <- x[c(FALSE,TRUE)]
    return( ((sum((x1-mean(x))*(x2-mean(x)))/(length(x1)-1))/(sd(x1)*sd(x2)))**2 )
}
#
fun <- function(x) {
    x1 <- x[c(TRUE,FALSE)]
    x2 <- x[c(FALSE,TRUE)]
    #x1 <- x1/sum(x1)
    #x2 <- x2/sum(x2)
    ( sum((x1-mean(x))*(x2-mean(x))) / ((length(x1)-1)*sd(x1)*sd(x2)) )**2
}


### repeatability pSTAT5 MFI
pdf('~nikolas/Thesis/figures/repeatability-PSTAT5-MFI.pdf',width=10, height=8)
figure.labels <- iter(paste(letters,')',sep=''))
par(mfrow=c(2,2))
# nn.peak.pstat5mfi
d <- ddply( nn.peak.pstat5mfi,c('cell.type'), function(x) {
      print(x$cell.type)
      x <- x[which(x$individual %in% REPEATS$individual),]
      x <- x[which(x$date %in% c(REPEATS$day1, REPEATS$day2)),]
      x <- x[order(x$individual),]
      print(x)
      y <- sapply(grep('^norm2\\.',colnames(x)), function(i) fun(x[,i]))
      names(y) <- DOSES
      y
}) 
fun.plot(d,'NN.peak pSTAT5 MFI') 
# raw.pstat5 mfi
d <- ddply( raw.pstat5mfi,c('cell.type'), function(x) {
      print(x$cell.type)
      x <- x[order(x$individual),]
      print(x)
      y <- sapply(grep('^X',colnames(x)), function(i) fun(x[,i]))
      names(y) <- DOSES
      y
}) 
fun.plot(d,'pSTAT5 MFI')
#p1 <- ggplot(melt(d),aes(x=cell.type,y=value,fill=variable))+geom_bar(position='dodge')+ggtitle('pSTAT5 MFI')+scale_fill_manual(values=blues4)+ylim(0,1)+ylab(expression(r^2))
# baseline raw.pstat5 mfi
d <- ddply( baseline.raw.pstat5mfi,c('cell.type'), function(x) {
      x <- x[order(x$individual),]
      y <- sapply(grep('^X',colnames(x)), function(i) fun(x[,i]))
      names(y) <- DOSES
      y
}) 
fun.plot(d,'base pSTAT5 MFI')
#p2 <- ggplot(melt(d),aes(x=cell.type,y=value,fill=variable))+geom_bar(position='dodge')+ggtitle('pSTAT5 MFI base')+scale_fill_manual(values=blues4)+ylim(0,1)+ylab(expression(r^2))
# pstat5 mfi
d <- ddply( nn.pstat5mfi,c('cell.type'), function(x) {
      x <- x[order(x$individual),]
      y <- sapply(grep('PSTAT5',colnames(x)), function(i) fun(x[,i]))
      names(y) <- DOSES
      y
}) 
fun.plot(d,'NN pSTAT5 MFI')
#p3 <- ggplot(melt(d),aes(x=cell.type,y=value,fill=variable))+geom_bar(position='dodge')+ggtitle('pSTAT5 MFI NN')+scale_fill_manual(values=blues4)+ylim(0,1)+ylab(expression(r^2))
# baseline pstat5 mfi
d <- ddply( baseline.nn.pstat5mfi,c('cell.type'), function(x) {
      x <- x[order(x$individual),]
      y <- sapply(grep('PSTAT5',colnames(x)), function(i) fun(x[,i]))
      y <- c(NA,y)
      names(y) <- DOSES
      y
}) 
fun.plot(d,'NN base pSTAT5 MFI')
#p4 <- ggplot(melt(d),aes(x=cell.type,y=value,fill=variable))+geom_bar(position='dodge')+ggtitle('pSTAT5 MFI NN base')+scale_fill_manual(values=blues4)+ylim(0,1)+ylab(expression(r^2))
#multiplot(p1, p2, p3, p4, cols=2)
dev.off()


### repeatability % pSTAT5+
pdf('~nikolas/Thesis/figures/repeatability-PSTAT5-pos.pdf',width=10,height=8)
figure.labels <- iter(paste(letters,')',sep=''))
par(mfrow=c(2,2))
# raw.pstat5 pos
d <- ddply( raw.pstat5pos,c('cell.type'), function(x) {
      x <- x[order(x$individual),]
      y <- sapply(grep('^X',colnames(x)), function(i) fun(x[,i]))
      names(y) <- DOSES
      y
}) 
fun.plot(d,'% pSTAT5+')
#p1 <- ggplot(melt(d),aes(x=cell.type,y=value,fill=variable))+geom_bar(position='dodge')+ggtitle('% pSTAT5+')+scale_fill_manual(values=blues4)+ylim(0,1)+ylab(expression(r^2))
#
# pstat5 pos
d <- ddply( nn.pstat5pos,c('cell.type'), function(x) {
      x <- x[order(x$individual),]
      y <- sapply(grep('PSTAT5',colnames(x)), function(i) fun(x[,i]))
      names(y) <- DOSES
      y
}) 
fun.plot(d,'NN % pSTAT5+')
#p2 <- ggplot(melt(d),aes(x=cell.type,y=value,fill=variable))+geom_bar(position='dodge')+ggtitle('% pSTAT5+ NN')+scale_fill_manual(values=blues4)+ylim(0,1)+ylab(expression(r^2))
#d$pheno <- 'pstat5+'
# baseline pstat5 pos
d <- ddply( baseline.nn.pstat5pos,c('cell.type'), function(x) {
      x <- x[order(x$individual),]
      y <- c(NA, sapply(grep('PSTAT5',colnames(x)), function(i) fun(x[,i])))
      names(y) <- DOSES
      y
}) 
fun.plot(d,'NN base % pSTAT5+')
#p3 <- ggplot(melt(d),aes(x=cell.type,y=value,fill=variable))+geom_bar(position='dodge')+ggtitle('% pSTAT5+ NN base')+scale_fill_manual(values=blues4)+ylim(0,1)+ylab(expression(r^2))
#multiplot(p1, p2, p3,cols=3)
dev.off()



# RMSD root mean square deviation
# pstat mfi 
cbind(
aggregate(PSTAT5.1 ~ cell.type, aggregate(PSTAT5.1 ~ individual + cell.type, pstat5mfi, function(x) (x[1]-x[2])**2), function(x) sqrt(mean(x))),
'PSTAT5.2'=aggregate(PSTAT5.2 ~ cell.type, aggregate(PSTAT5.2 ~ individual + cell.type, pstat5mfi, function(x) (x[1]-x[2])**2),  function(x) sqrt(mean(x)))[,'PSTAT5.2'],
'PSTAT5.3'=aggregate(PSTAT5.3 ~ cell.type, aggregate(PSTAT5.3 ~ individual + cell.type, pstat5mfi, function(x) (x[1]-x[2])**2),  function(x) sqrt(mean(x)))[,'PSTAT5.3'],
'PSTAT5.4'=aggregate(PSTAT5.4 ~ cell.type, aggregate(PSTAT5.4 ~ individual + cell.type, pstat5mfi, function(x) (x[1]-x[2])**2),  function(x) sqrt(mean(x)))[,'PSTAT5.4']
) 
# baseline pstat mfi 
cbind(
aggregate(diff.PSTAT5.2 ~ cell.type, aggregate(diff.PSTAT5.2 ~ individual + cell.type, baseline.pstat5mfi, function(x) (x[1]-x[2])**2), function(x) sqrt(mean(x))),
'diff.PSTAT5.3'=aggregate(diff.PSTAT5.3 ~ cell.type, aggregate(diff.PSTAT5.3 ~ individual + cell.type, baseline.pstat5mfi, function(x) (x[1]-x[2])**2),  function(x) sqrt(mean(x)))[,'diff.PSTAT5.3'],
'diff.PSTAT5.4'=aggregate(diff.PSTAT5.4 ~ cell.type, aggregate(diff.PSTAT5.4 ~ individual + cell.type, baseline.pstat5mfi, function(x) (x[1]-x[2])**2),  function(x) sqrt(mean(x)))[,'diff.PSTAT5.4']
) 
# pstat pos 
cbind(
aggregate(PSTAT5.1 ~ cell.type, aggregate(PSTAT5.1 ~ individual + cell.type, pstat5pos, function(x) (x[1]-x[2])**2), function(x) sqrt(mean(x))),
'PSTAT5.2'=aggregate(PSTAT5.2 ~ cell.type, aggregate(PSTAT5.2 ~ individual + cell.type, pstat5pos, function(x) (x[1]-x[2])**2),  function(x) sqrt(mean(x)))[,'PSTAT5.2'],
'PSTAT5.3'=aggregate(PSTAT5.3 ~ cell.type, aggregate(PSTAT5.3 ~ individual + cell.type, pstat5pos, function(x) (x[1]-x[2])**2),  function(x) sqrt(mean(x)))[,'PSTAT5.3'],
'PSTAT5.4'=aggregate(PSTAT5.4 ~ cell.type, aggregate(PSTAT5.4 ~ individual + cell.type, pstat5pos, function(x) (x[1]-x[2])**2),  function(x) sqrt(mean(x)))[,'PSTAT5.4']
) 
# pstat pos baseline corrected
cbind(
aggregate(diff.PSTAT5.2 ~ cell.type, aggregate(diff.PSTAT5.2 ~ individual + cell.type, baseline.pstat5pos, function(x) (x[1]-x[2])**2), function(x) sqrt(mean(x))),
'diff.PSTAT5.3'=aggregate(diff.PSTAT5.3 ~ cell.type, aggregate(diff.PSTAT5.3 ~ individual + cell.type, baseline.pstat5pos, function(x) (x[1]-x[2])**2),  function(x) sqrt(mean(x)))[,'diff.PSTAT5.3'],
'diff.PSTAT5.4'=aggregate(diff.PSTAT5.4 ~ cell.type, aggregate(diff.PSTAT5.4 ~ individual + cell.type, baseline.pstat5pos, function(x) (x[1]-x[2])**2),  function(x) sqrt(mean(x)))[,'diff.PSTAT5.4']
)


