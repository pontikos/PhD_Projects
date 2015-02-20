setwd('~nikolas/stats-archive/Papers/KIR/') 
source('Code/plate.R') 


d.noqc <- zero.cn.trans(load.experiment(experiment.file='Data/qPCR/processed/experiment.csv', QC=NULL, trans=function(x)-x))
#d <- zero.cn.trans(load.experiment(experiment.file='Data/qPCR/processed/experiment.csv', QC=TRUE, trans=identity))

#d2 <- zero.cn.trans(load.experiment(experiment.file='Data/qPCR/processed/experiment.csv', QC=TRUE, trans=function(x) 2**(-x)))
d <- zero.cn.trans(load.experiment(experiment.file='Data/qPCR/processed/experiment.csv', QC=qc.experiment, trans=function(x)-x))



pdf('bmc_article/figures/KIR3DS1-preQC.pdf', width=6, height=3)
plot.experiment(field='DS1.median', experiment=d.noqc, xlab='', ylab='', panel.extra=panel.extra.peaknorm, strip=strip.custom(bg='lightblue'), scales=list(x=list(draw=F)), ylim=c(-2,max(d.noqc[,'DS1.median'])+1))
dev.off() 
pdf('bmc_article/figures/KIR3DS1-postQC.pdf', width=6, height=3)
plot.experiment(field='DS1.median', experiment=d, xlab='', ylab='', panel.extra=panel.extra.peaknorm, strip=strip.custom(bg='lightblue'), scales=list(x=list(draw=F)), ylim=c(-2,max(d[,'DS1.median'])+1))
dev.off() 
pdf('bmc_article/figures/KIR3DS1-peaknormalised.pdf', width=6, height=3)
plot.experiment(field='DS1.peaknorm.median', experiment=d, plot.fun=stripplot, xlab='', ylab='', strip=strip.custom(bg='lightblue'), scales=list(x=list(draw=F)), panel.extra=panel.extra.abline) 
dev.off()



pdf('bmc_article/figures/KIR3DL1-preQC.pdf', width=6, height=3)
plot.experiment(field='DL1.median', experiment=d.noqc, ylab='', xlab='', panel.extra=panel.extra.peaknorm, strip=strip.custom(bg='lightblue'), scales=list(x=list(draw=F)), ylim=c(-2,max(d.noqc[,'DL1.median'])+1))
dev.off()
pdf('bmc_article/figures/KIR3DL1-postQC.pdf', width=6, height=3)
plot.experiment(field='DL1.median', experiment=d, ylab='', xlab='', panel.extra=panel.extra.peaknorm, strip=strip.custom(bg='lightblue'), scales=list(x=list(draw=F)), ylim=c(-2,max(d[,'DL1.median'])+1))
dev.off()
pdf('bmc_article/figures/KIR3DL1-peaknormalised.pdf', width=6, height=3)
plot.experiment(field='DL1.peaknorm.median', experiment=d, plot.fun=stripplot, xlab='', ylab='', strip=strip.custom(bg='lightblue'), scales=list(x=list(draw=F)), panel.extra=panel.extra.abline) 
dev.off()

# repeats need to be summarised by their median
library(dplyr)
dim( d.norepeats <- ddply(d, .(sample), summarise, DL1.peaknorm.median=median(DL1.peaknorm.median), DS1.peaknorm.median=median(DS1.peaknorm.median)) )

X <- function(d, normtrans='peak') d[, paste(c('DS1', 'DL1'), sprintf('%snorm.median',normtrans), sep=".")]
plot(X(d), pch=20)
plot(X(d.norepeats), pch=20)
plot(d[,c('DS1.median','DL1.median')], pch=20, xlim=c(-2,2), ylim=c(-2,2))

genotypes <- as.matrix( expand.grid(x=0:3,y=0:3) ) [c(2,3,4,5,6,7,9,10),] 
res <- kmeans(X((d)), centers=genotypes)
res <- data.frame(X(d), cluster=res$cluster)
mu <- lapply(sort(unique(res$cluster)), function(i) colMeans(res[res$cluster==i,1:2]))
Sigma <- lapply(sort(unique(res$cluster)), function(i) cov(res[res$cluster==i,1:2]))
lambda <- as.numeric(table(res$cluster))/sum(as.numeric(table(res$cluster)))
plot(X((d)), pch=as.character(as.factor(res$cluster)), main='Both')
library(car)
#for ( i in 1:length(unique((res$cluster))) ) car::ellipse(mu[[i]], Sigma[[i]], 1, col='blue')
#Sigma[[5]] <- Sigma[[7]]
# assign same covariance to 1-2 and 2-1
Sigma[[ which.min(lapply(mu, function(x) sum(abs(x-c(2,1))))) ]] <- Sigma[[ which.min(lapply(mu, function(x) sum(abs(x-c(1,2))))) ]]
for ( i in 1:length(unique((res$cluster))) ) car::ellipse(mu[[i]], Sigma[[i]], 1, col='red')


library(mixtools)

### initialise with k-means then refine parameters with EM
res <- mvnormalmixEM( X(d), mu=mu, lambda=lambda, sigma=Sigma, verb=TRUE )
#pdf('../lc_scale_up/mvnormalmixEM.pdf')
mixtools::plot.mixEM(res, whichplots=2)
#dev.off()

initial.params <- list(mu=mu, lambda=lambda, sigma=Sigma)

save(initial.params, file='initial.params.obj')

load_all( '~nikolas/Projects/KIR/bin/caller' ) 


locs <- clusterdef.ncopies(3, R.) 


library(mixsmsn)
res <- smsn.mmix(X(d), error=10E-10, iter.max=1000, mu=mu, pii=lambda, Sigma=Sigma, shape=lapply(1:8, function(j) mu[[j]]/2), g=8, family='Skew.normal', group=TRUE, get.init=F, obs.prob=TRUE)


ggplot(X(d), aes(x=DS1.peaknorm.median, y=DL1.peaknorm.median)) + geom_point() +
geom_point(data=, col='red')

matrix.sqrt <- function(A) {
  sva <- svd(A)
  if (min(sva$d)>=0)
    Asqrt <- t(sva$v %*% (t(sva$u) * sqrt(sva$d)))
  else
    stop("Matrix square root is not defined")
  return(Asqrt)
}

## Densidade/CDF da SN com locação escala #######
dmvSN <- function(y, mu, Sigma, lambda){
  #y: deve ser uma matrix onde cada linha tem um vetor de dados multivariados de dimensão ncol(y) = p. nrow(y) = tamanho da amostra
  #mu, lambda: devem ser do tipo vetor de mesma dimensão igual a ncol(y) = p
  #Sigma: Matrix p x p
  n <- nrow(y)
  p <- ncol(y)
  dens <- 2*dmvnorm(y, mu, Sigma)*pnorm( apply( matrix(rep(t(lambda)%*%solve(matrix.sqrt(Sigma)),n), n, p, byrow = TRUE)*(y - matrix(rep(mu, n), n, p, byrow = TRUE)), 1, sum  ) )
  return(dens)
}


###########    Densidades das Misturas de SNI   ##################
d.mixedmvSN <- function(y, pi1, mu, Sigma, lambda){
    #y: é a matriz de dados
    #pi1: deve ser do tipo vetor de dimensão g
    #mu: deve ser do tipo list com g entradas. Cada entrada do list deve ser um vetor de dimensão p
    #Sigma: deve ser do tipo list com g entradas. Cada entrada do list deve ser uma matriz p x p
    #lambda: deve ser do tipo list com g entradas. Cada entrada do list deve ser um vetor de dimensão p
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j]*dmvSN(y, mu[[j]], Sigma[[j]], lambda[[j]])
    return(dens)
}

bw=.01
gr <- expand.grid(seq(-2,3.5,bw),seq(-2,3.5,bw))
dens <- sapply(1:length(res$mu), function(j) res$pii[[j]]*dmvSN(gr, res$mu[[j]], res$Sigma[[j]], res$shape[[j]]))
#posteriors <- t(apply(dens, 1, function(x) x/sum(x)))
plot(X(d), pch=20, col=res$group)
points( data.frame(do.call('rbind',mu)), pch='x', col='black', cex=1)
points( data.frame((do.call('rbind',res$mu))), pch='x', col='black', cex=2)
#points(gr, pch=20, cex=.5)

for (j in 1:length(res$mu)) {
gr.sub <- gr[dens[,j]>quantile(dens[,j],probs=(seq(0,1,.01)))[['90%']],]
lines(gr.sub[chull(gr.sub),])
}


cumsum(dens[,1])


points(gr[which(posteriors>.6,arr.ind=T)[,1],], pch=20, cex=.5)
d.mixedmvSN(gr, res$pii, res$mu, res$Sigma, res$shape) 
for ( i in 1:length(unique((res$group))) ) car::ellipse(res$mu[[i]], res$Sigma[[i]], 1, col='red')



library(ggplot2)
pdf('bmc_article/figures/Figure-S2.pdf', width=8, height=8)
ggplot(data=transform(d,t1d=ifelse(t1d==0,'controls','cases'),box384=ifelse(box384=='leftovers',box384,paste('plate',box384,sep=':'))), aes(x=DS1.peaknorm.median, y=DL1.peaknorm.median, group=sample, col=c(blue='control',red='case',black='calibrator')[col])) + facet_wrap(~box384+t1d) + geom_point() + scale_color_manual(values=c(control='blue',case='red',calibrator='black')) + ylab(expression(paste('Peak Normalised ', italic('KIR3DL1 '), Delta, 'Ct'))) + xlab(expression(paste('Peak Normalised ', italic('KIR3DS1 '), Delta, 'Ct'))) + guides(colour=FALSE)
dev.off()

posteriors <- read.csv('Data/qPCR/processed/posteriors.csv')
names(posteriors) <- sub('X(.).(.)', '\\1-\\2', names(posteriors), fixed=F)

posteriors$CN <- apply(posteriors[,grep('\\d-\\d',colnames(posteriors))],1, function(x) names(x)[which.max(x)])

# certainty with which each CN genotype is called
for (geno in grep('\\d-\\d',colnames(posteriors), value=T)) cat(geno, 100*length(which(posteriors[posteriors$CN==geno,geno]>.95)) / length(which(posteriors$CN==geno)), '\n' )


kir.genotypes <- read.csv('Data/qPCR/processed/kir-genotypes.csv')
       
anova(glm(t1d ~ KIR, data=kir.genotypes, family='binomial'),test="Chisq")

