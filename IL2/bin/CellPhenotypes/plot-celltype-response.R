library(scales)
library(ggplot2)
suppressPackageStartupMessages(library("optparse"))

option_list <- list( 
make_option(c("-K",'--clusters'), default='', help = 'number of clusters K')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

K <- opt$clusters
K <- 8

base.dir <- '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3'
print( flowclust.dir <- file.path(base.dir,sprintf('CellPhenotypes/flowClust-%s',K)) )

print(dim(d <- do.call('rbind',lapply(list.files(file.path(flowclust.dir,'CSV'),pattern='.*.csv',full.names=TRUE), read.csv, stringsAsFactors=FALSE))))
#median.pstat5 <- paste('pheno.median.pstat5',1:4,sep='.')
median.pstat5 <- paste('pheno.pstat5.median.diff',1:4,sep='.')

print(dim(d <- merge(d,read.csv(file.path(base.dir,'cb-IL2RA.csv'),stringsAsFactors=FALSE))))

print(load('~/thor/geno-pheno.RData'))
d$rs41295055 <- geno[d$individual,'rs41295055']
d$rs2104286 <- geno[d$individual,'rs2104286']
d$rs12722495 <- geno[d$individual,'rs12722495']

#i<-1
#boxplot(d[which(d$cluster==i),median.pstat5], col=i)
#for (i in 1:5) boxplot(d[which(d$cluster==i),median.pstat5], col=i, add=TRUE)

colMedian <- function(x) apply(x, 2, median) 

pstat5.median <- do.call('rbind', by(d, d$cluster, function(x) colMedian(x[,median.pstat5]))) 
X <- pstat5.median


channels <- c('CD25','CD45RA','CD4','FOXP3','weights')
clusters.mfi <- do.call('rbind', by(d, d$cluster, function(x) colMedian(x[,channels]))) 
Y <- clusters.mfi

library(fmsb)
Y <- data.frame(rbind(apply(clusters.mfi,2,function(x) c(max(x),min(x))),clusters.mfi))
radarchart(Y,axistype=0,seg=5,plty=1,plwd=4) 
#stars(clusters.mfi, locations = c(0, 0), radius = FALSE, key.loc = c(0, 0), main = "", col.lines=1:4, lwd=2)


#days colour by individual
d$date <- as.Date(d$date)
g <- ggplot(d, aes(x=date, y=pheno.median.pstat5.2-pheno.median.pstat5.1,group=individual)) + scale_x_date(labels=date_format( "%b-%Y")) + xlab("") 
g <- g + geom_point()+geom_line(data=d, aes(x=date, y=pheno.median.pstat5.2-pheno.median.pstat5.1, group=individual, col=individual))
g <- g + facet_grid(cluster ~ .) +guides(colour=FALSE)
print(g)


#days colour by genotype
d$date <- as.Date(d$date)
g <- ggplot(d, aes(x=date, y=pheno.median.pstat5.2-pheno.median.pstat5.1,group=individual,col=rs2104286)) + scale_x_date(labels=date_format( "%b-%Y")) + xlab("") 
g <- g + geom_point()
g <- g + facet_grid(cluster ~ .)
print(g)

# transform d to d2 for plotting loess with ggplot2
d2 <- melt(d,grep('pheno.median.pstat5',colnames(d),value=T,invert=T))
d2$dose <- as.character(d2$variable)
d2[which(d2$dose=='pheno.median.pstat5.1'),'dose'] <- '0'
d2[which(d2$dose=='pheno.median.pstat5.2'),'dose'] <- '0.1'
d2[which(d2$dose=='pheno.median.pstat5.3'),'dose'] <- '10'
d2[which(d2$dose=='pheno.median.pstat5.4'),'dose'] <- '1000'
d2$median.pstat5 <- d2$value
d2$dose <- as.numeric(d2$dose)

g <- ggplot(d2, aes(x=as.factor(dose), y=median.pstat5, group=as.factor(cluster), col=as.factor(cluster))) + xlab('dose')
g <- g + geom_point() + geom_smooth()
#g <- g + facet_grid(rs41295055 ~ .)
#g <- g+ facet_grid(rs2104286 ~ .)
g <- g + scale_colour_manual(values = c('black', "red","green", "blue"))
print(g)

# plot dose against median.pstat5 group by cell type and facet by genotype
g <- ggplot(d2, aes(x=as.factor(dose), y=median.pstat5, group=as.factor(cluster), col=as.factor(cluster))) + xlab('dose')
g <- g + geom_point() + geom_smooth()
g <- g + facet_grid(rs41295055 ~ .)
#g <- g+ facet_grid(rs2104286 ~ .)
g <- g + scale_colour_manual(values = c('black', "red","green", "blue"))
print(g)
# plot dose against median.pstat5 group by genotype and facet by cell type
g <- ggplot(d2, aes(x=as.factor(dose), y=median.pstat5, group=rs41295055, col=rs41295055)) + xlab('dose')
g <- g + geom_point() + geom_smooth()
g <- g + facet_grid(cluster ~ .)
print(g)
#
g <- ggplot(d2, aes(x=as.factor(dose), y=median.pstat5, group=rs2104286, col=rs2104286)) + xlab('dose')
g <- g + geom_point() + geom_smooth()
g <- g + facet_grid(cluster ~ .)
print(g)
# plot dose against median.pstat5 group by t1d and facet by cell type
g <- ggplot(d2, aes(x=as.factor(dose), y=median.pstat5, group=as.factor(t1d), col=as.factor(t1d))) + xlab('dose')
g <- g + geom_point() + geom_smooth()
g <- g + facet_grid(cluster ~ .)
print(g) 
# plot dose against median.pstat5 group by sex and facet by cell type
g <- ggplot(d2, aes(x=as.factor(dose), y=median.pstat5, group=as.factor(sex), col=as.factor(sex))) + xlab('dose')
g <- g + geom_point() + geom_smooth()
g <- g + facet_grid(cluster ~ .)
print(g)


g <- ggplot(d, aes(x=rs41295055,y=CD25))
g <- g + geom_boxplot() + geom_jitter() + coord_flip()
g <- g + facet_grid( cluster ~ ., scales='free')
print(g)




r<-(d[which(d$individual %in% names(which(table(d$individual)>4))),])
r <- r[order(r$individual, r$date),]


pdf(file.path(flowclust.dir,'all_pstat5.pdf'))
# all
plot(1:4, X[1,], ylim=range(X), col='white', xlab='IL-2', ylab='median.pstat5', lwd=2, main='All')
sapply(unique(d$cluster), function(i) lines(1:4, X[i,median.pstat5],col=i)) 
pstat5.range <- by(d, d$cluster, function(x) t(sapply(median.pstat5, function(i)quantile(x[,i],probs=c(.1,.90)))))
Y <- pstat5.range 
for (j in sort(unique(d$cluster))) {
  for (i in 1:4) {
    a <- Y[[j]][i,]
    print(a)
    segments(i,a[1],i,a[2],col=j)
  }
} 
#rs41295055
AA.d <- d[which(d$rs41295055=='C/C'),]
AG.d <- d[which(d$rs41295055=='C/T'),]
GG.d <- d[which(d$rs41295055=='T/T'),]
AA.pstat5.median <- do.call('rbind', by(AA.d, AA.d$cluster, function(x) colMedian(x[,median.pstat5]))) 
AA.X <- AA.pstat5.median 
AG.pstat5.median <- do.call('rbind', by(AG.d, AG.d$cluster, function(x) colMedian(x[,median.pstat5]))) 
AG.X <- AG.pstat5.median 
GG.pstat5.median <- do.call('rbind', by(GG.d, GG.d$cluster, function(x) colMedian(x[,median.pstat5]))) 
GG.X <- GG.pstat5.median 
plot(1:4, X[1,], ylim=range(AA.X), col='white', xlab='IL-2', ylab='median.pstat5', lwd=2, main='rs41295055')
sapply(unique(AA.d$cluster), function(i) lines(1:4, AA.X[i,median.pstat5],col=i,lty=2)) 
pstat5.range <- by(AA.d, AA.d$cluster, function(x) t(sapply(median.pstat5, function(i)quantile(x[,i],probs=c(.1,.90)))))
AA.Y <- pstat5.range 
for (j in unique(AA.d$cluster)) {
  for (i in 1:4) {
    a <- AA.Y[[j]][i,]
    print(a)
    segments(i,a[1],i,a[2],col=j,lty=2)
  }
} 
sapply(unique(AG.d$cluster), function(i) lines(1:4, AG.X[i,median.pstat5],col=i)) 
pstat5.range <- by(AG.d, AG.d$cluster, function(x) t(sapply(median.pstat5, function(i)quantile(x[,i],probs=c(.1,.90)))))
AG.Y <- pstat5.range 
for (j in unique(AG.d$cluster)) {
  for (i in 1:4) {
    a <- AG.Y[[j]][i,]
    print(a)
    segments(i,a[1],i,a[2],col=j)
  }
} 
#rs2104286
AA.d <- d[which(d$rs2104286=='A/A'),]
AG.d <- d[which(d$rs2104286=='A/G'),]
GG.d <- d[which(d$rs2104286=='G/G'),]
AA.pstat5.median <- do.call('rbind', by(AA.d, AA.d$cluster, function(x) colMedian(x[,median.pstat5]))) 
AA.X <- AA.pstat5.median 
AG.pstat5.median <- do.call('rbind', by(AG.d, AG.d$cluster, function(x) colMedian(x[,median.pstat5]))) 
AG.X <- AG.pstat5.median 
GG.pstat5.median <- do.call('rbind', by(GG.d, GG.d$cluster, function(x) colMedian(x[,median.pstat5]))) 
GG.X <- GG.pstat5.median 
plot(1:4, X[1,], ylim=range(AA.X), col='white', xlab='IL-2', ylab='median.pstat5', lwd=2, main='rs2104286')
sapply(unique(AA.d$cluster), function(i) lines(1:4, AA.X[i,median.pstat5],col=i,lty=2)) 
pstat5.range <- by(AA.d, AA.d$cluster, function(x) t(sapply(median.pstat5, function(i)quantile(x[,i],probs=c(.1,.90)))))
AA.Y <- pstat5.range 
for (j in unique(AA.d$cluster)) {
  for (i in 1:4) {
    a <- AA.Y[[j]][i,]
    print(a)
    segments(i,a[1],i,a[2],col=j,lty=2)
  }
} 
sapply(unique(AG.d$cluster), function(i) lines(1:4, AG.X[i,median.pstat5],col=i)) 
pstat5.range <- by(AG.d, AG.d$cluster, function(x) t(sapply(median.pstat5, function(i)quantile(x[,i],probs=c(.1,.90)))))
AG.Y <- pstat5.range 
for (j in unique(AG.d$cluster)) {
  for (i in 1:4) {
    a <- AG.Y[[j]][i,]
    print(a)
    segments(i,a[1],i,a[2],col=j)
  }
} 
#rs12722495
AA.d <- d[which(d$rs12722495=='A/A'),]
AG.d <- d[which(d$rs12722495=='A/G'),]
GG.d <- d[which(d$rs12722495=='G/G'),]
AA.pstat5.median <- do.call('rbind', by(AA.d, AA.d$cluster, function(x) colMedian(x[,median.pstat5]))) 
AA.X <- AA.pstat5.median 
AG.pstat5.median <- do.call('rbind', by(AG.d, AG.d$cluster, function(x) colMedian(x[,median.pstat5]))) 
AG.X <- AG.pstat5.median 
GG.pstat5.median <- do.call('rbind', by(GG.d, GG.d$cluster, function(x) colMedian(x[,median.pstat5]))) 
GG.X <- GG.pstat5.median 
plot(1:4, X[1,], ylim=range(AA.X), col='white', xlab='IL-2', ylab='median.pstat5', lwd=2, main='rs12722495')
sapply(unique(AA.d$cluster), function(i) lines(1:4, AA.X[i,median.pstat5],col=i,lty=2)) 
pstat5.range <- by(AA.d, AA.d$cluster, function(x) t(sapply(median.pstat5, function(i)quantile(x[,i],probs=c(.1,.90)))))
AA.Y <- pstat5.range 
for (j in unique(AA.d$cluster)) {
  for (i in 1:4) {
    a <- AA.Y[[j]][i,]
    print(a)
    segments(i,a[1],i,a[2],col=j,lty=2)
  }
} 
sapply(unique(AG.d$cluster), function(i) lines(1:4, AG.X[i,median.pstat5],col=i)) 
pstat5.range <- by(AG.d, AG.d$cluster, function(x) t(sapply(median.pstat5, function(i)quantile(x[,i],probs=c(.1,.90)))))
AG.Y <- pstat5.range 
for (j in unique(AG.d$cluster)) {
  for (i in 1:4) {
    a <- AG.Y[[j]][i,]
    print(a)
    segments(i,a[1],i,a[2],col=j)
  }
} 
# case-control
#case-control
cases.d <- d[which(d$t1d==2),]
cases.pstat5.median <- do.call('rbind', by(cases.d, cases.d$cluster, function(x) colMedian(x[,median.pstat5]))) 
cases.X <- cases.pstat5.median 
controls.d <- d[which(d$t1d==1),]
controls.pstat5.median <- do.call('rbind', by(controls.d, controls.d$cluster, function(x) colMedian(x[,median.pstat5]))) 
controls.X <- controls.pstat5.median 
plot(1:4, X[1,], ylim=range(cases.X), col='white', xlab='IL-2', ylab='median.pstat5', lwd=2, main='case-control')
sapply(unique(cases.d$cluster), function(i) lines(1:4, cases.X[i,median.pstat5],col=i,lty=2)) 
pstat5.range <- by(cases.d, cases.d$cluster, function(x) t(sapply(median.pstat5, function(i)quantile(x[,i],probs=c(.1,.90)))))
cases.Y <- pstat5.range 
for (j in unique(cases.d$cluster)) {
  for (i in 1:4) {
    a <- cases.Y[[j]][i,]
    print(a)
    segments(i,a[1],i,a[2],col=j,lty=2)
  }
} 
sapply(unique(controls.d$cluster), function(i) lines(1:4, controls.X[i,median.pstat5],col=i)) 
pstat5.range <- by(controls.d, controls.d$cluster, function(x) t(sapply(median.pstat5, function(i)quantile(x[,i],probs=c(.1,.90)))))
controls.Y <- pstat5.range 
for (j in unique(controls.d$cluster)) {
  for (i in 1:4) {
    a <- controls.Y[[j]][i,]
    print(a)
    segments(i,a[1],i,a[2],col=j)
  }
} 
# male-female
males.d <- d[which(d$sex==1),]
males.pstat5.median <- do.call('rbind', by(males.d, males.d$cluster, function(x) colMedian(x[,median.pstat5]))) 
males.X <- males.pstat5.median 
females.d <- d[which(d$t1d==2),]
females.pstat5.median <- do.call('rbind', by(females.d, females.d$cluster, function(x) colMedian(x[,median.pstat5]))) 
females.X <- females.pstat5.median 
plot(1:4, X[1,], ylim=range(males.X), col='white', xlab='IL-2', ylab='median.pstat5', lwd=2, main='male-female')
sapply(unique(males.d$cluster), function(i) lines(1:4, males.X[i,median.pstat5],col=i,lty=2)) 
pstat5.range <- by(males.d, males.d$cluster, function(x) t(sapply(median.pstat5, function(i)quantile(x[,i],probs=c(.1,.90)))))
males.Y <- pstat5.range 
for (j in unique(males.d$cluster)) {
  for (i in 1:4) {
    a <- males.Y[[j]][i,]
    print(a)
    segments(i,a[1],i,a[2],col=j,lty=2)
  }
} 
sapply(unique(females.d$cluster), function(i) lines(1:4, females.X[i,median.pstat5],col=i)) 
pstat5.range <- by(females.d, females.d$cluster, function(x) t(sapply(median.pstat5, function(i)quantile(x[,i],probs=c(.1,.90)))))
females.Y <- pstat5.range 
for (j in unique(females.d$cluster)) {
  for (i in 1:4) {
    a <- females.Y[[j]][i,]
    print(a)
    segments(i,a[1],i,a[2],col=j)
  }
} 
#
for (individual in unique(d$individual)) {
  print(individual)
  X2 <- d[which(individual==d$individual),c('individual','date','cluster',median.pstat5)]
  for (Date in unique(X2$date)) {
    print(Date)
    X <- X2[which(Date==X2$date),c('individual','date','cluster',median.pstat5)]
    print(X)
    plot(factor(1:4), X[1,median.pstat5], col='white', ylim=range(X[,median.pstat5]), main=paste(unique(X$individual), unique(X$date)), ylab='median.pstat5', xlab='IL-2')
    sapply(X$cluster,  function(i) lines(factor(1:4),X[i,median.pstat5],col=i))
  }
}
dev.off()


