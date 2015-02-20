setwd('~/Projects/IL2')



#Neil Walker, 4th Jan 2014
#First (?) attempt to make relevant data selection for a subset of
#D-GAP and IP subjects used in Tony's experiments, for Niko
#-----------------------------------------------------------------
#Briefly, I extracted all local (e.g. Taqman) data for CBR-IP and DGAP
#boxsets, for the region (build GRCh37) chromosome Hs10:6037726..6156736
#and then reduced it to the necessary subjects.

#NB [1]: the datasets do not have the same SNPs in them
#NB [2]: one source ID cannot be matched - KM00860R


# 75 snps in IL2RA from Chris
dim(snps <- read.csv('IL2RA-genotype/snps-for-niko-dec-2013.csv'))
length(snps$position.hg19)

#Lookup file to give relationship between Niko's ID - a mixture
#of samples, visits and study subjects, with a few typos -
#a corrected ID, and the subject ID as found on the data files 
dim(lookup <- read.csv('IL2RA-genotype/il2ra-lookup-2014-01-04.csv'))

#Data and support file for genotyping performed in D-GAP subjects (original IDs beginning KM or KA):
dim(dgap <- read.csv('IL2RA-genotype/dgap-IL2RA-2013-12-23.csv'))
dim(dgap.support <- read.csv('IL2RA-genotype/dgap-IL2RA-2013-12-23-support.csv'))

#Data and support file for genotyping performed in IP subjects (original IDs beginning I or CB):
dim(cb <- read.csv('IL2RA-genotype/ip-IL2RA-2013-12-23.csv'))
dim(cb.support <- read.csv('IL2RA-genotype/ip-IL2RA-2013-12-23-support.csv'))


length(intersect(snps$position.hg19, dgap.support$start)) 
length(intersect(snps$snp_id, cb.support$dbSNP))
length(intersect(snps$snp_id, dgap.support$dbSNP))

       
snps.cb <- cb.support[cb.support$start %in% intersect(snps$position.hg19, cb.support$start),] 
snps.dgap <- dgap.support[dgap.support$start %in% intersect(snps$position.hg19, dgap.support$start),]


sapply(snps.cb$variant, function(x) grep(x,colnames(cb)))

cb <- cb[ , c('uniqueID','t1d','sex',sort(c(paste(snps.cb$variant, '1', sep='_'), paste(snps.cb$variant, '2', sep='_'))))] 
dgap <- dgap[ , c('uniqueID','t1d','sex',sort(c(paste(snps.dgap$variant, '1', sep='_'), paste(snps.dgap$variant, '2', sep='_'))))]


dim(cb <- read.csv('IL2RA-genotype/ip-IL2RA-2013-12-23.csv'))
cb<-cbind(cb[,c('uniqueID','t1d','sex')],mapply(paste, cb[,paste(snps.cb$variant, '1', sep='_')], cb[,paste(snps.cb$variant, '2', sep='_')],MoreArgs=list(sep='-')))
colnames(cb)<-gsub('_1', '', colnames(cb))
cb<-merge(lookup,cb)
colnames(cb)[grep('DIL',colnames(cb))] <- as.character(snps.cb[ snps.cb$variant %in% grep('DIL',colnames(cb),value=T), 'dbSNP'])
write.csv(cb,file='IL2RA-genotype/cb-IL2RA.csv', quote=FALSE,row.names=FALSE) 
#cb[which(cb=='N-N',arr.ind=T)] <- NA 
#sort(apply(cb[,grep('DIL',colnames(cb))], 2, function(x) sum(x!='N-N')))


dim(dgap <- read.csv('IL2RA-genotype/dgap-IL2RA-2013-12-23.csv'))
dgap<-cbind(dgap[,c('uniqueID','t1d','sex')],mapply(paste, dgap[,paste(snps.dgap$variant, '1', sep='_')], dgap[,paste(snps.dgap$variant, '2', sep='_')],MoreArgs=list(sep='-')))
colnames(dgap)<-gsub('_1', '', colnames(dgap)) 
dgap<-merge(lookup,dgap)
colnames(dgap)[grep('DIL',colnames(dgap))] <- as.character(snps.dgap[ snps.dgap$variant %in% grep('DIL',colnames(dgap),value=T), 'dbSNP'])
write.csv(dgap,file='IL2RA-genotype/dgap-IL2RA.csv',quote=FALSE,row.names=FALSE)
#dgap[which(dgap=='N-N',arr.ind=T)] <- NA
#sort(apply(dgap[,grep('DIL',colnames(dgap))], 2, function(x) sum(x!='N-N')))



length(snps.cb[grep('DIL', colnames(cb), value=T) %in% snps.cb$variant, 'dbSNP'])
length(snps.dgap[grep('DIL', colnames(dgap), value=T) %in% snps.dgap$variant, 'dbSNP'])


pstat5 <- read.csv('IL2RA-genotype/cb-dgap-pstat5-peaks.csv', sep=',')
#
cb <- read.table('IL2RA-genotype/cb-FCS-IL2RA.csv',sep=';',header=T)
cb$dose <- gsub('U','',as.character(cb$dose))
cb$dose <- gsub('01','0.1',cb$dose)
cb$dose <- factor(cb$dose,levels=c('0','0.1','10','1000'))
dim(cb <- merge(cb,pstat5))

cb$ids.date.panel <- paste(cb$ids, cb$date, cb$panel, sep=',')

#only keep those which have a single dose per individual,date,panel
dim(cb <- cb[cb$ids.date.panel %in% names(which(rowSums(do.call('rbind',by(cb, cb$ids.date.panel, function(x) table(x$dose))))==4)),])


ggplot(data=cb[which(as.numeric(as.character(cb$dose)) > 0),], aes(x=log10(as.numeric(as.character(dose))), y=pstat5.loc.1, group=ids.date.panel, col=ids.date.panel)) + geom_point() + guides(colour=FALSE) + geom_line( data=cb[as.numeric(as.character(cb$dose)) > 0,], aes(x=log10(as.numeric(as.character(dose))), y=pstat5.loc.1, group=ids.date.panel, col=ids.date.panel) )
ggplot(data=cb[which(as.numeric(as.character(cb$dose)) > 0),], aes(x=log10(as.numeric(as.character(dose))), y=pstat5.dens.ratio, group=ids.date.panel, col=ids.date.panel)) + geom_point() + guides(colour=FALSE) + geom_line( data=cb[as.numeric(as.character(cb$dose)) > 0,], aes(x=log10(as.numeric(as.character(dose))), y=pstat5.dens.ratio, group=ids.date.panel, col=ids.date.panel) )

ggplot(data=cb[which(as.numeric(as.character(cb$dose)) > 0),], aes(x=log10(as.numeric(as.character(dose))), y=pstat5.dens.height.2-pstat5.dens.height.1, group=ids.date.panel, col=ids.date.panel)) + geom_point() + guides(colour=FALSE) + geom_line( data=cb[as.numeric(as.character(cb$dose)) > 0,], aes(x=log10(as.numeric(as.character(dose))), y=pstat5.dens.height.2-pstat5.dens.height.1, group=ids.date.panel, col=ids.date.panel) )
ggplot(data=cb[which(as.numeric(as.character(cb$dose)) > 0),], aes(x=log10(as.numeric(as.character(dose))), y=pstat5.loc.2-pstat5.loc.1, group=ids.date.panel, col=ids.date.panel)) + geom_point() + guides(colour=FALSE) + geom_line( data=cb[as.numeric(as.character(cb$dose)) > 0,], aes(x=log10(as.numeric(as.character(dose))), y=pstat5.loc.2-pstat5.loc.1, group=ids.date.panel, col=ids.date.panel) )


#
dgap <- read.table('dgap-FCS-IL2RA.csv',sep=';',header=T)
dgap$dose <- gsub('U','',as.character(dgap$dose))
dgap$dose <- gsub('01','0.1',dgap$dose)
dgap$dose <- factor(dgap$dose,levels=c('0','0.1','10','1000'))
dgap <- merge(dgap,pstat5)

library(reshape)
m <- melt(cb, grep('pstat5', names(cb), invert=T, value=T))

gplotRegression <- function (fit) {
    require(ggplot2)
    ggplot(data=fit$model, aes(x =as.numeric(names(fit$model)[2]), y = as.numeric(names(fit$model)[1]))) + 
     geom_jitter(position=position_jitter(width=.1)) +
     geom_smooth(method = "lm", col = "red") +
     ggtitle(paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                         "; Intercept =",signif(fit$coef[[1]],5 ),
                         "; Slope =",signif(fit$coef[[2]], 5),
                         "; P =",signif(summary(fit)$coef[2,4], 5)))
}
#cb
ggplotRegression(lm(y ~ x, data=transform(cb[as.numeric(as.character(cb$dose)) > 0,], x=log10(as.numeric(as.character(dose))), y=pstat5.loc.1)))+labs(list(x='dose',y='pstat5.loc.1'))
ggplotRegression(lm(y ~ x, data=transform(cb[as.numeric(as.character(cb$dose)) > 0,], x=log10(as.numeric(as.character(dose))), y=pstat5.loc.1)))+labs(list(x='dose',y='pstat5.loc.2'))

ggplotRegression(lm(y ~ x, data=transform(cb[as.numeric(as.character(cb$dose)) > 0,], x=log10(as.numeric(as.character(dose))), y=pstat5.dens.height.1)))+labs(list(x='dose',y='pstat5.dens.height.1'))
ggplotRegression(lm(y ~ x, data=transform(cb[as.numeric(as.character(cb$dose)) > 0,], x=log10(as.numeric(as.character(dose))), y=pstat5.dens.height.2)))+labs(list(x='dose',y='pstat5.dens.height.2'))


#dgap
ggplotRegression(lm(y ~ x, data=transform(dgap[as.numeric(as.character(dgap$dose)) > 0,], x=log10(as.numeric(as.character(dose))), y=pstat5.loc.1)))




par(mfrow=c(2,2))

with( cb[as.numeric(as.character(cb$dose)) > 0,],
     plot(log10(as.numeric(as.character(dose))), pstat5.loc.1) )


lm_eqn = function(df){
  m <- lm(y ~ x, df)
  eq <- substitute(italic(y) == a + b %.% italic(x)*", "~~italic(r)^2~"="~r2, list(a = format(coef(m)[1], digits = 2), b = format(coef(m)[2], digits = 2), r2 = format(summary(m)$r.squared, digits = 3)))
  return(as.character(as.expression(eq)))
}


library(ggplot2)
ggplot(data=cb[as.numeric(as.character(cb$dose)) > 0,], aes(x=log10(as.numeric(as.character(dose))), y=pstat5.loc.1)) + geom_jitter(position=position_jitter(width=.1)) + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)

ggplot(data=cb[as.numeric(as.character(cb$dose)) > 0,], aes(x=log10(as.numeric(as.character(dose))), y=pstat5.loc.1)) + geom_point() + geom_smooth(method = "lm", color="black")


ggplot(data=cb[which(as.numeric(as.character(cb$dose)) > 0 & cb$panel == 'CD25,CD4,CD45RA,FOXP3,PSTAT5'),], aes(x=log10(as.numeric(as.character(dose))), y=pstat5.loc.1, group=id.date, col=id.date)) + geom_point() + guides(colour=FALSE) + geom_line( data=cb[as.numeric(as.character(cb$dose)) > 0,], aes(x=log10(as.numeric(as.character(dose))), y=pstat5.loc.1, group=id.date, col=id.date) )

ggplot(data=cb[which(as.numeric(as.character(cb$dose)) > 0 & cb$panel == 'CD25,CD4,CD45RA,FOXP3,PSTAT5'),], aes(x=log10(as.numeric(as.character(dose))), y=pstat5.dens.height.1, group=id.date, col=id.date)) + geom_point() + guides(colour=FALSE) + geom_line( data=cb[as.numeric(as.character(cb$dose)) > 0,], aes(x=log10(as.numeric(as.character(dose))), y=pstat5.dens.height.1, group=id.date, col=id.date) )


g <- g + geom_point(data=r, aes(x=date, y=(APC.mef), group=individual, col=individual))+geom_line(data=r, aes(x=date, y=(APC.mef), group=individual, col=individual),size=.25)+guides(colour=FALSE)

table(subset(cb, panel='CD25,CD4,CD45RA,FOXP3,PSTAT5')$dose)


with( ,
     {
      m <- lm(y ~ x)
      eq <- substitute(italic(y) == a + b %.% italic(x)*", "~~italic(r)^2~"="~r2, list(a = format(coef(m)[1], digits = 2), b = format(coef(m)[2], digits = 2), r2 = format(summary(m)$r.squared, digits = 3)))
      eq <- as.character(as.expression(eq))
      qplot(x=x, y=y, geom='point')+geom_smooth(method='lm')+xlab(eq)
      annotate("text", x=0.5, y=15000, label=lm_eqn(df), hjust=0, size=8, 
             family="Times", face="italic", parse=TRUE)
     }
)

#+geom_text(aes(x = 25, y = 300, label = lm_eqn(data)), parse = TRUE)

abline(lm(pstat5.loc.1 ~ as.numeric(dose), cb))

plot(cb$dose, cb$pstat5.dens.height.1)
plot(cb$dose, cb$pstat5.loc.2)
plot(cb$dose, cb$pstat5.dens.height.2)

par(mfrow=c(2,1))
plot(cb$date, cb$pstat5.peak1)
plot(cb$date, cb$pstat5.peak2)


#'0U','01U','10U','100U','1000U','10000U'


# Debbie: cb snps
# https://intranet-dil/index.php/en/concrete-proposals/197-genotyping-of-cbr-samples-for-qc
#DIL rs # gene MAF
#DIL90 rs3087243 CTLA4 0.500
#DIL8083 rs2476601 PTPN22 0.092
#
#DIL8103 rs2104286 IL2RA 0.236
#DIL10847 rs11594656 IL2RA 0.287
#DIL9620 rs12722495 IL2RA 0.092
#DIL8092 rs12722489 IL2RA 0.144
#
#DIL14110 rs2187668 DR3 tag 0.103
#DIL14316 rs2844821 A2 tag 0.362
#DIL14385 rs660895 DR4 tag 0.287
#DIL14383 rs9271366 DR15 tag 0.149
#DIL14607 rs1150743 A24 tag 0.098
#DIL14459 rs2394185 A24 tag 0.069
#DIL7550 rs17612648 CD45 0.011
#DIL13098 rs45450798 PTPN2 0.126
#DIL13050 rs478582 PTPN2 0.489
#DIL14022 rs8192284 IL6R 0.356

# Debbie: dgap snps
# https://intranet-dil/index.php/en/completed/137-d-gap
#DIL rs number gene
#DIL14077 rs3024505 IL10
#DIL14110 rs2187668 DR3 tag
#DIL8083 rs2476601 PTPN22
##DIL11843 rs1990760 IFIH1
#DIL90 rs3087243 CTLA4
#
#DIL9620 rs12722495 IL2RA
#DIL8103 rs2104286 IL2RA
#DIL10847 rs11594656 IL2RA
#DIL8092 rs12722489 IL2RA
#
#DIL969 rs689 INS
#DIL13098 rs45450798 PTPN2
#DIL13050 rs478582 PTPN2
#DIL14316 rs2844821 A2 tag
#DIL14385 rs660895 DR4 tag
#DIL7550 rs17612648 CD45
#DIL14383 rs9271366 DR15 tag
#DIL14607 rs1150743 A24 tag
#DIL14459 rs2394185 A24 tag
#DIL12228 rs17884228 HLAB39-in3

library(beeswarm)
library(nlme)
setwd('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/')


#mclust stat5
ext <- 'pstat5.clusters.csv'
individual.date <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/individual-date.csv', col.names=c('individual','date'))
X <- NULL
for (i in 1:nrow(individual.date)) {
    individual <- individual.date[i,'individual']
    date <- individual.date[i,'date']
    f <- file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/Lymphocytes/pstat5-ann-response/RData',paste(individual, date, ext, sep='_'))
    d <- read.csv(f,col.names=c('0U','01U','10U','1000U','cluster'))
    X <- rbind(X, data.frame(individual=individual, date=date, cluster=sort(unique(d$cluster)), pstat5.response=tapply(d[,'X1000U'], d$cluster, mean)))
}
write.csv(X, file='~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',quote=FALSE,row.names=FALSE)



#lymph
head(pstat5 <- read.csv('lymph-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','pstat5.response')))
head(pstat5norm <- read.csv('lymph-pstat5norm-response.csv',stringsAsFactors=F,col.names=c('id','date','pstat5norm.response')))

#naive
dim(pstat5 <- read.csv('naive-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','pstat5.response')))
dim(pstat5norm <- read.csv('naive-pstat5norm-response.csv',stringsAsFactors=F,col.names=c('id','date','pstat5norm.response')))

#memory
head(pstat5 <- read.csv('memory-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','pstat5.response')))
head(pstat5norm <- read.csv('memory-pstat5norm-response.csv',stringsAsFactors=F,col.names=c('id','date','pstat5norm.response')))

#ann mean response 
head(pstat5 <- read.csv('lymph-pstat5-ann-mean-response.csv',stringsAsFactors=F,col.names=c('id','date','pstat5.response')))
head(pstat5norm <- read.csv('lymph-pstat5-ann-mean-response-norm.csv',stringsAsFactors=F,col.names=c('id','date','pstat5norm.response')))

#ann sum response 
head(pstat5 <- read.csv('lymph-pstat5-ann-sum-response.csv',stringsAsFactors=F,col.names=c('id','date','pstat5.response')))
head(pstat5norm <- read.csv('lymph-pstat5-ann-sum-response-norm.csv',stringsAsFactors=F,col.names=c('id','date','pstat5norm.response')))


 
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/lymph_median_pstat5_response.csv',stringsAsFactors=F))
head(pstat5norm <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/lymph_median_pstat5-norm_response.csv',stringsAsFactors=F))
head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))
dim(d <- merge(snps, merge(pstat5, pstat5norm)))
length(dup.individuals <- names(which(table(d$individual)==2)))
length(unique(d$individual)) 
d[d=='N-N'] <- NA 
snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric)))) 
d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric) 
for (p in c('pstat5.response','pstat5.norm.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(individual), data=d)
    m <- lme(fm,random=~ 1|as.factor(individual), data=na.omit(d[,c(p,'individual', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$individual %in% dup.individuals, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]
head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F)) 
#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$individual)==2)))
length(unique(d$individual)) 
#dim(d2 <- d[-which(duplicated(d$id)),]) 
d[d=='N-N'] <- NA 
snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric)))) 
#apply(d[,snps], 2, table) 
#table(d$snp) 
d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric) 
#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(individual), data=na.omit(d[,c(p,'individual', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$individual %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}



#mclust mean response
head(pstat5 <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/mclust-pstat5-response.csv',stringsAsFactors=F,col.names=c('id','date','cluster','pstat5.response')))
pstat5 <- subset(pstat5,cluster==3)
pstat5 <- pstat5[,-3]

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))

#dim(d <- merge(snps, merge(pstat5, pstat5norm)))
dim(d <- merge(snps, pstat5))
length(dup.ids <- names(which(table(d$id)==2)))
length(unique(d$id)) 
#dim(d2 <- d[-which(duplicated(d$id)),])

d[d=='N-N'] <- NA

snps <- grep('^rs', colnames(d), value=T)
snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))

#apply(d[,snps], 2, table) 
#table(d$snp)

d[,snps] <- lapply(d[,snps], factor)
d[,snps] <- lapply(d[,snps], as.numeric)


#for (p in c('pstat5.response','pstat5norm.response')) {
for (p in c('pstat5.response')) {
par(mfrow=c(4,5))
#title(main=p, outer=TRUE)
for (snp in snps) {
    fm <- as.formula(paste(p, snp, sep=' ~ '))
    #m <- lme(fm,random=~ 1|as.factor(id), data=d)
    m <- lme(fm,random=~ 1|as.factor(id), data=na.omit(d[,c(p,'id', snp)]))
    #main <- tryCatch(sprintf('pvalue=%.3f', anova(m)[['Pr(>F)']][1]), error=function(e) {''}, finally='')
    main <- tryCatch(sprintf('pvalue=%.3f', summary(m)$tTable[,'p-value'][[2]]), warning=function(e) {''}, finally='')
    beeswarm(fm,data=d, main=main, ylim=range(d[,p]))
}
x <- d[ d$id %in% dup.ids, p]
x1 <- x[c(TRUE,FALSE)]
x2 <- x[c(FALSE,TRUE)]
plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',cor(x1,x2)))
print(d[9:10,])
points(t(d[9:10,p]),col='red')
abline(b=1,a=0)
}




#library(ggplot2)
#ggplot(data=d, aes(x=pstat5.response, y=pstat5norm.response, group=id)) + geom_point() + geom_line(col=as.factor(d$id))




require('kernlab')
 
kfunction <- function(linear =0, quadratic=0)
{
  k <- function (x,y)
 {
     linear*sum((x)*(y)) + quadratic*sum((x^2)*(y^2))
  }
  class(k) <- "kernel"
  k
}
                       
svp <- ksvm(x,y,type="C-svc",C = 100, kernel=kfunction(1,0),scaled=c())
plot(c(min(x[,1]), max(x[,1])),c(min(x[,2]), max(x[,2])),type='n',xlab='x1',ylab='x2')
title(main='Linear Separable Features')
ymat <- ymatrix(svp)
points(x[-SVindex(svp),1], x[-SVindex(svp),2], pch = ifelse(ymat[-SVindex(svp)] < 0, 2, 1))
points(x[SVindex(svp),1], x[SVindex(svp),2], pch = ifelse(ymat[SVindex(svp)] < 0, 17, 16))
    
# Extract w and b from the model   
w <- colSums(coef(svp)[[1]] * x[SVindex(svp),])
b <- b(svp)
    
# Draw the lines
abline(b/w[2],-w[1]/w[2])
abline((b+1)/w[2],-w[1]/w[2],lty=2)




