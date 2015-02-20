library(flowCore)
library(MASS)


read.FCS('~/dunwich/FCS/cad100_2008aug05_treg_i009121n_013.fcs')->fcs2

pdf('fcs2.pdf')
plot(fcs2, c('APC-A', 'Pacific Blue-A'), xlim=c(1,100), ylim=c(1,100))
contour(fcs2, c('APC-A', 'Pacific Blue-A'), add=T)
dev.off()


pdf('fcs2.trans.pdf')
fcs2.trans <- transform(fcs2, `Pacific Blue-A`=logTrans(`Pacific Blue-A`), `APC-A`=logTrans(`APC-A`))
plot(fcs2.trans, c('APC-A', 'Pacific Blue-A'), xlim=c(0,1), ylim=c(0,1))
contour(fcs2.trans, c('APC-A', 'Pacific Blue-A'), xlim=c(0,1), ylim=c(0,1), add=T)
dev.off()



read.FCS('~/dunwich/FCS.Marcin/cad100_2008aug05_treg_i009121n_013.fcs')->fcs3

pdf('fcs3.asinh.pdf')
#plot(fcs3, c('APC-A', 'Pacific Blue-A'), xlim=c(min(fcs3@exprs[,5]),10), ylim=c(min(fcs3@exprs[,-8]),10))
#plot(fcs3, c('APC-A', 'Pacific Blue-A'), xlim=c(-100,1000), ylim=c(-100,1000))
#kde2d(fcs3.trans@exprs[,5], fcs3.trans@exprs[,8])->fcs3.dens
#contour(fcs3, c('APC-A', 'Pacific Blue-A'), add=T)
asinh(fcs3@exprs[,c('APC-A', 'Pacific Blue-A')])->trans
densCols(trans)->cols
plot(trans, pch=".", col=cols)
dev.off()

pdf('fcs3.density.pdf')
#plot(fcs3, c('APC-A', 'Pacific Blue-A'), xlim=c(min(fcs3@exprs[,5]),10), ylim=c(min(fcs3@exprs[,-8]),10))
#plot(fcs3, c('APC-A', 'Pacific Blue-A'), xlim=c(-100,1000), ylim=c(-100,1000))
#kde2d(fcs3.trans@exprs[,5], fcs3.trans@exprs[,8])->fcs3.dens
#contour(fcs3, c('APC-A', 'Pacific Blue-A'), add=T)
par(mfrow=c(2,1))
plot(density(fcs3@exprs[,c('APC-A')]), xlim=c(-100,10000))
plot(density(fcs3@exprs[,c('Pacific Blue-A')]), xlim=c(-100,10000))
dev.off()

pdf('fcs3.nocomp.pdf')
#plot(fcs3, c('APC-A', 'Pacific Blue-A'), xlim=c(min(fcs3@exprs[,5]),10), ylim=c(min(fcs3@exprs[,-8]),10))
plot(fcs3, c('APC-A', 'Pacific Blue-A'), xlim=c(-100,1000), ylim=c(-100,1000))
kde2d(fcs3@exprs[,5], fcs3.trans@exprs[,8], lims=c(c(-100,1000), c(-100,1000)) )->fcs3.dens
#contour(fcs3, c('APC-A', 'Pacific Blue-A'), add=T)
contour(fcs3.dens, xlim=c(-100,1000), ylim=c(-100,1000), add=T)
#points(fcs3@exprs[,c('APC-A', 'Pacific Blue-A')], pch=".")
dev.off()


pdf('fcs3.comp.pdf')
compensate(fcs3,fcs3@description[["SPILL"]])->fcs3.comp
plot(fcs3.comp, c('APC-A', 'Pacific Blue-A'), xlim=c(-100,1000), ylim=c(-100,1000))
contour(fcs3.comp, c('APC-A', 'Pacific Blue-A'), add=T)
#points(fcs3.comp@exprs[,c('APC-A', 'Pacific Blue-A')], pch=".")
dev.off()


pdf('fcs3.comp.trans.pdf')
lgcl <- estimateLogicle(fcs3.comp, channels = c( 'APC-A', 'Pacific Blue-A' ))
fcs3.trans <- transform(fcs3.comp, lgcl)
plot(fcs3.trans, c('APC-A', 'Pacific Blue-A'), xlim=c(0,1.5), ylim=c(0,3))
contour(fcs3.trans, c('APC-A', 'Pacific Blue-A'), xlim=c(0,1.5), ylim=c(0,3), add=T)
dev.off()




pdf('fcs2.full.trans.pdf')
fcs2.trans <- transform(fcs2, `Pacific Blue-A`=logTrans(`Pacific Blue-A`), `APC-A`=logTrans(`APC-A`))
plot(fcs2.trans, c('APC-A', 'Pacific Blue-A'), xlim=c(0,2), ylim=c(0,2))
#kde2d(fcs2.trans@exprs[,5], fcs2.trans@exprs[,8])->fcs2.dens
contour(fcs2.trans, c('APC-A', 'Pacific Blue-A'), xlim=c(0,2), ylim=c(0,2), add=T)
dev.off()

pdf('fcs3.full.trans.pdf')
lgcl <- estimateLogicle(fcs3.comp, channels = c( 'APC-A', 'Pacific Blue-A' ))
fcs3.trans <- transform(fcs3.comp, lgcl)
plot(fcs3.trans, c('APC-A', 'Pacific Blue-A'), xlim=c(0,3), ylim=c(0,3))
#kde2d(fcs3.trans@exprs[,5], fcs3.trans@exprs[,8])->fcs3.dens
contour(fcs3.trans, c('APC-A', 'Pacific Blue-A'), xlim=c(0,3), ylim=c(0,3), add=T)
dev.off()

pdf('fcs3.subset.pdf')
plot(fcs3.comp@exprs[,5], fcs3.trans@exprs[,5], pch=".", xlim=c(-3000,3000))
#plot(fcs3.trans[1:1000,], c('APC-A', 'Pacific Blue-A'))
lines(sort(fcs3.comp@exprs[,5]), sort(fcs3.trans@exprs[,5]), col="blue")
sinh.inv <- function(x, a, b) a*log(x+sqrt(x**2+1), base=10)+b
lines(0:3000, log10(0:3000))
lines(-3000:0, -log10(abs(-3000:0)))
#lines(0:3000, .5*log10(0:3000))
dev.off()

pdf('fcs3.subset.pdf')
plot(fcs3.comp@exprs[,5], fcs3.trans@exprs[,5], pch=".", xlim=c(-10,10))
lines(sort(fcs3.comp@exprs[,5]), sort(fcs3.trans@exprs[,5]), col="blue")
lines(-1:1, -1:1)
dev.off()


function(X, T, M, W)

logicle <-  function(x, W, M) (10 ** (x-W) - 10 ** -(x-W)) * T * 10 ** -(M-W)
