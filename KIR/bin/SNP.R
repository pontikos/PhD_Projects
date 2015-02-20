library(ggplot2)
library(mitools)
library(mice)
library(xtable)
setwd('~nikolas/stats-archive/Papers/KIR/') 

dduplicated <- function(df) duplicated(df) | duplicated(df, fromLast = TRUE)
intersection <- function(x, y, ...){if (missing(...)) intersect(x, y) else intersect(x, intersection(y, ...))}

# geno
GENO.COL <- data.frame( geno=c("0-2", "1-1", "2-0", "2-1", "1-2", "0-1", "1-0", "3-0"), col=c('darkgreen', 'pink', 'blue', 'lightblue', 'purple', 'green', 'red', 'grey') )
row.names(GENO.COL) <- GENO.COL$geno
GENO <- GENO.COL$geno

#lookup to revtrieve dil_subjectid from sampleid
head(lookup <- read.csv('Data/ic-cc-clean-lookup-2013-02-27.csv'))

# HLA data
head(hla.data <- read.table('Data/HLA/cc-hla-2013-06-21.tab',header=T))
dim( hla.data )
hla.data$t1d <- hla.data$t1d-1
hla.data$dil_subjectid <- hla.data$dil_subject
dim( hla.data <- subset(hla.data, t1d %in% 0:1) )
hla.data$HLAA <- with(hla.data,paste(HLAA_Bw4_aa80_1, HLAA_Bw4_aa80_2, sep='-'))
hla.data$HLAB <- with(hla.data,paste(HLAB_Bw4_Bw6_aa80_1, HLAB_Bw4_Bw6_aa80_2, sep='-'))
hla.data$HLA_Bw <- ifelse(grepl('6N', with(hla.data, paste(HLAA,HLAB,sep='-'))), '6N', '0')
hla.data$HLA_Bw[grepl('4T', with(hla.data, paste(HLAA,HLAB,sep='-')))] <- '4T'
hla.data$HLA_Bw[grepl('4I', with(hla.data, paste(HLAA,HLAB,sep='-')))] <- '4I'
#drop column collection
head(hla.data <- hla.data[,-which(colnames(hla.data)=='collection')])
# 20,445 = 9,174 : 11,271  


# disagreement of knn predictions between LOO and whole dataset
print(load('~nikolas/stats-archive/Papers/KIR/Data/SNP/loo-knn-calls-for-training-data-6feb2014.RData')) 
loo.knn <- lapply(loo.knn, function(loo.knn) {
 rownames(loo.knn) <- sapply( strsplit(row.names( loo.knn ), '\\.'), '[[', 2)
 loo.knn$sample <- rownames(loo.knn)
 loo.knn
})
print(load('~nikolas/stats-archive/Papers/KIR/Data/SNP/KNN-10-datasets.RData'))
imp.loo <- mapply(function(a, b) merge(a,b), imputed.datasets, loo.knn, SIMPLIFY=FALSE)
sapply(imp.loo , function(x) table(x$knn.call!=x$y.knn))
lapply(imp.loo , function(x) x[x$knn.call!=x$y.knn,])
#do.call('rbind',lapply(imp.loo , function(x) x[x$knn.call!=x$y.knn,]))
m <- do.call('intersection', lapply(imp.loo , function(x) x[x$knn.call!=x$y.knn,'sample']))
dim(unique(do.call('rbind',lapply(imp.loo , function(x) x[x$sample %in% m,]))))


#posteriors
head(qpcr.posteriors <- read.csv('Data/qPCR/processed/posteriors.csv'))
qpcr.posteriors <- merge(qpcr.posteriors, lookup)
colnames(qpcr.posteriors) <- gsub('sampleid', 'sample', colnames(qpcr.posteriors))
qpcr.cn <- grep('^X', colnames(qpcr.posteriors), value=T)
qpcr.posteriors$qpcr.cn <- qpcr.cn[ apply(qpcr.posteriors[,qpcr.cn],1,which.max) ]
qpcr.posteriors$qpcr.cn <- gsub('X(.).(.)', '\\1-\\2', qpcr.posteriors$qpcr.cn)
print( colnames(qpcr.posteriors) <- gsub('X(.).(.)', '\\1-\\2', colnames(qpcr.posteriors)) )
rownames(qpcr.posteriors) <- qpcr.posteriors$sample
head( qpcr.posteriors )


### takes a list of imputed datasets and returns the combined dataset
combine.imputations <- function(imputations) {
    imp <- imputations[[1]]
    imp$KIR <- factor(apply(sapply(imputations, function(x) x$KIR),1,function(y)names(which.max(table(y)))), levels=GENO.COL$geno)
    imp$KIR3DL1 <- factor(gsub('.-(.)', '\\1', imp$KIR), levels=c('2','1','0'))
    imp$KIR3DS1 <- factor(gsub('(.)-.', '\\1', imp$KIR), levels=c('0','1','2','3'))
    return(imp)
}


#these are the LOO datasets, containing the KNN predictions for the qPCR data and the qPCR predictions
print(load('Data/SNP/loo-knn-calls-for-training-data-6feb2014.RData')) 
make.loo.knn <- function(loo.knn) {
#loo.knn
 rownames(loo.knn) <- sapply( strsplit(row.names( loo.knn ), '\\.'), '[[', 2)
 loo.knn$sample <- rownames(loo.knn)
 loo.knn <- merge(loo.knn, lookup, by.x='sample', by.y='sampleid')
 #use the qPCR predictions
 loo.knn$KIR <- factor(gsub('cn(.)(.)', '\\1-\\2', loo.knn$y), levels=GENO)
 loo.knn$KIR3DS1 <- factor(gsub('(.)-.', '\\1', loo.knn$KIR), levels=c('0','1','2','3'))
 loo.knn$KIR3DL1 <- factor(gsub('.-(.)', '\\1', loo.knn$KIR), levels=c('2','1','0'))
 loo.knn$kir3dl1 <- ifelse(as.numeric(gsub('.-(.)','\\1', loo.knn[,'KIR3DL1'])>0), 'KIR3DL1+', 'KIR3DL1-')
 loo.knn$kir3ds1 <- ifelse(as.numeric(gsub('(.)-.','\\1', loo.knn[,'KIR3DS1'])>0), 'KIR3DS1+', 'KIR3DS1-')
 #loo predictions
 loo.knn$loo.KIR <- factor(gsub('cn(.)(.)', '\\1-\\2', loo.knn$y.knn), levels=GENO)
 loo.knn$loo.KIR3DS1 <- factor(gsub('(.)-.', '\\1', loo.knn$loo.KIR), levels=c('0','1','2','3'))
 loo.knn$loo.KIR3DL1 <- factor(gsub('.-(.)', '\\1', loo.knn$loo.KIR), levels=c('2','1','0'))
 loo.knn$t1d <- ifelse(grepl('cases',loo.knn$collection),1,0)
 loo.knn <- loo.knn[,c('dil_subjectid','sample','t1d','KIR','loo.KIR','KIR3DS1','KIR3DL1')]
 loo.knn <- merge(loo.knn, hla.data, all.x=TRUE)
 rownames(loo.knn) <- loo.knn$sample
 return(loo.knn)
}
loo.knn <- lapply(loo.knn, make.loo.knn) 
comb.loo.knn <- combine.imputations(loo.knn)

# check qPCR calls from posteriors match qPCR calls in loo object
if ( !with( merge(qpcr.posteriors[,c('sample','qpcr.cn')], comb.loo.knn), table(qpcr.cn==KIR) ) == dim(comb.loo.knn)[[1]] ) stop('error')


print(load('Data/SNP/KNN-10-datasets.RData'))
#these are the knn imputed datasets which include the qPCR predictions
make.imputed.datasets <- function(imputed.dataset, loo.knn) {
 imputed.dataset <- merge(imputed.dataset, lookup, by.x='sample', by.y='sampleid')
 rownames(imputed.dataset) <- imputed.dataset$sample
 imputed.dataset$knn.call <- factor(gsub('cn(.)(.)', '\\1-\\2', imputed.dataset$knn.call), levels=GENO)
 #check if the loo.knn KNN predictions match the imputed.dataset KNN predictions
 print(with(merge(imputed.dataset,loo.knn, by='sample'), table(knn.call!=loo.KIR)))
 imputed.dataset <- imputed.dataset[!duplicated(imputed.dataset$dil_subjectid),]
 imputed.dataset[rownames(loo.knn), 'knn.call'] <- loo.knn[,'loo.KIR']
 #print(with(merge(imputed.dataset,loo.knn, by='sample'), table(knn.call!=loo.KIR)))
 colnames(imputed.dataset) <- gsub('knn.call', 'loo.KIR', colnames(imputed.dataset))
 #loo predictions
 imputed.dataset$KIR <- imputed.dataset$loo.KIR
 imputed.dataset[rownames(loo.knn), 'KIR'] <- loo.knn[,'KIR']
 #print(with(merge(imputed.dataset,loo.knn) , table(loo.KIR.1 != loo.KIR.2)))
 #correct the datasets to the qPCR cn prediction instead of knn predictions
 imputed.dataset$KIR3DL1 <- factor(gsub('.-(.)', '\\1', imputed.dataset$KIR), levels=c('2','1','0'))
 imputed.dataset$KIR3DS1 <- factor(gsub('(.)-.', '\\1', imputed.dataset$KIR), levels=c('0','1','2','3'))
 imputed.dataset$t1d <- ifelse(grepl('cases',imputed.dataset$collection),1,0)
 imputed.dataset$kir3dl1 <- ifelse(as.numeric(gsub('.-(.)','\\1',imputed.dataset[,'KIR3DL1'])>0), 'KIR3DL1+', 'KIR3DL1-')
 imputed.dataset$kir3ds1 <- ifelse(as.numeric(gsub('(.)-.','\\1',imputed.dataset[,'KIR3DS1'])>0), 'KIR3DS1+', 'KIR3DS1-')
 # only 12,050 out of 12,106 have HLA information
 imputed.dataset <- merge(imputed.dataset, hla.data, all.x=TRUE)
 rownames(imputed.dataset) <- imputed.dataset$sample
 return(imputed.dataset)
}

length(imputed.datasets <- mapply(make.imputed.datasets, imputed.datasets, loo.knn, SIMPLIFY=FALSE))
comb.imputed.datasets <- combine.imputations(imputed.datasets)

qpcr.imputed.datasets <- lapply(imputed.datasets, function(x) x[rownames(qpcr.posteriors),])


### ODD RATIOS
logistic <- function(geno, imputations) {
    #multiple imputation
    require(mitools)
    #snp or qpcr
    X <- mitools::imputationList(imputations)
    x <- mitools::MIcombine(with(X, glm(as.formula(paste('t1d',geno,sep=' ~ ')), family=binomial)))
    s <- summary(x)
    rownames(s) <- levels(imputations[[1]][,geno])
    s <- data.frame(cbind((exp(cbind(OR=coef(x), confint(x)))), pvalue=2*(1-pnorm(abs(s$results/s$se)))), row.names=gsub('geno', '', rownames(s)))
    colnames(s) <- c('OR', 'Lower', 'Upper', 'p-value')
    s[1,] <- c('OR'=1, 'Lower'=0, 'Upper'=0, 'p-value'=0)
    #order rows by same order as factor
    s <- s[levels(imputations[[1]][,geno]),]
    s$vars <- factor(levels(imputations[[1]][,geno]),levels=levels(imputations[[1]][,geno]))
    s$OR <- round(s$OR, digits=2)
    s$Lower <- round(s$Lower, digits=2)
    s$Upper <- round(s$Upper, digits=2)
    s$`p-value` <- round(s$`p-value`, digits=4)
    return(s)
} 

assoc.pool.compare <- function(geno, imputations) {
    #multiple imputation
    require(mice)
    #snp or qpcr
    fit1 <- as.mira( with(imputationList(imputations), glm(as.formula(paste('t1d',geno,sep=' ~ ')), family=binomial)) )
    fit0 <- as.mira( with(imputationList(imputations), glm(t1d ~ 1), family=binomial) )
    round(mice::pool.compare(fit1=fit1, fit0=fit0)$pvalue, digits=4) 
} 

caco.table <- function(geno, imp) {
    caco <- as.data.frame(table(imp[,geno], imp$t1d)[,c(2,1)])
    colnames(caco) <- c('case','control')
    caco$total <- rowSums(caco)
    caco <- rbind(caco, total=colSums(caco))
    return(caco)
}

generate.assoc.table <- function(geno, imputations) {
    caco <- caco.table(geno, combine.imputations(imputations))
    caco <- caco[-nrow(caco),]
    m <- logistic(geno, imputations)
    write.table(data.frame("case:control"=with(caco,paste(case,control,sep=':')), "total"=caco$total, "OR"=m[,'OR'], "95%CI"=with(m,paste(Lower,Upper,sep='-')), "p-value"=m[,'p-value'], check.names=FALSE), sep='&', quote=FALSE, row.names=rownames(caco))
}

#caco.table('KIR', combine.imputations('KIR', imputed.datasets))
#caco.table('KIR3DL1', combine.imputations('KIR3DL1', imputed.datasets))


### Table-1
#table1.a
#qPCR
generate.assoc.table('KIR', qpcr.imputed.datasets)
assoc.pool.compare('KIR', qpcr.imputed.datasets)
#SNP
generate.assoc.table('KIR', imputed.datasets)
assoc.pool.compare('KIR', imputed.datasets)
#table1.b
#qPCR
generate.assoc.table('KIR', qpcr.imputed.datasets)
assoc.pool.compare('KIR', qpcr.imputed.datasets)
#SNP
generate.assoc.table('KIR3DL1', imputed.datasets)
assoc.pool.compare('KIR3DL1', imputed.datasets)
#table1.c
generate.assoc.table('KIR3DS1', imputed.datasets)
assoc.pool.compare('KIR3DS1', imputed.datasets)
#table2.a
#qPCR
generate.assoc.table('KIR', lapply(qpcr.imputed.datasets, function(x) x[grep('4',x$HLA_Bw),]))
assoc.pool.compare('KIR', lapply(qpcr.imputed.datasets, function(x) x[grep('4',x$HLA_Bw),]))
#SNP
generate.assoc.table('KIR', lapply(imputed.datasets, function(x) x[grep('4',x$HLA_Bw),]))
assoc.pool.compare('KIR', lapply(imputed.datasets, function(x) x[grep('4',x$HLA_Bw),]))
#table2.b
#qPCR
generate.assoc.table('KIR3DL1', lapply(qpcr.imputed.datasets, function(x) x[grep('4',x$HLA_Bw),]))
assoc.pool.compare('KIR3DL1', lapply(qpcr.imputed.datasets, function(x) x[grep('4',x$HLA_Bw),]))
#SNP
generate.assoc.table('KIR3DL1', lapply(imputed.datasets, function(x) x[grep('4',x$HLA_Bw),]))
assoc.pool.compare('KIR3DL1', lapply(imputed.datasets, function(x) x[grep('4',x$HLA_Bw),]))
#table2.c
#qPCR
generate.assoc.table('KIR3DS1', lapply(qpcr.imputed.datasets, function(x) x[grep('4I',x$HLA_Bw),]))
assoc.pool.compare('KIR3DS1', lapply(qpcr.imputed.datasets, function(x) x[grep('4I',x$HLA_Bw),]))
#SNP
generate.assoc.table('KIR3DS1', lapply(imputed.datasets, function(x) x[grep('4I',x$HLA_Bw),]))
assoc.pool.compare('KIR3DS1', lapply(imputed.datasets, function(x) x[grep('4I',x$HLA_Bw),]))

### Table 3: INTERACTION TESTS

make.cases <- function(d, qPCR=FALSE) {
    if (qPCR)
        cases <- subset(d[rownames(comb.loo.knn),], t1d==1 & HLA_Bw!='0')
    else
        cases <- subset(d, t1d==1 & HLA_Bw!='0')
    cases$kir3dl1 <- ifelse(as.numeric(gsub('.-(.)','\\1',cases[,'KIR3DL1'])>0), 'KIR3DL1+', 'KIR3DL1-')
    cases$kir3ds1 <- ifelse(as.numeric(gsub('(.)-.','\\1',cases[,'KIR3DS1'])>0), 'KIR3DS1+', 'KIR3DS1-')
    cases$kir <- with(cases, paste(kir3ds1, kir3dl1))
    #cases$HLABW4 <- c('HLA-Bw4-','HLA-Bw4+')[1+as.numeric(grepl('4',cases$HLA_Bw))]
    cases$HLABW4 <- as.numeric(grepl('4',cases$HLA_Bw))
    #cases$HLABW4I <- c('HLA-Bw4-80I-','HLA-Bw4-80I+')[1+as.numeric(grepl('4I',cases$HLA_Bw))]
    cases$HLABW4I <- as.numeric(grepl('4I',cases$HLA_Bw))
    return(cases)
}

all.cases <- lapply(imputed.datasets, make.cases)
all.comb.cases <- make.cases(combine.imputations(all.cases))

qpcr.cases <- lapply(imputed.datasets, make.cases, qPCR=TRUE)
qpcr.comb.cases <- make.cases(combine.imputations(qpcr.cases), qPCR=TRUE)

#library(Zelig)
#z <- zelig(HLABW4 ~ KIR, data=lapply(imputed.datasets, make.cases), model='ls', cite=FALSE)
#anova(z, test='Chisq')

require(mice)


# table3.a qPCR
with(qpcr.comb.cases, table(HLABW4,kir))
sum(with(qpcr.comb.cases, table(HLABW4,kir)))
fit1 <- as.mira(with(imputationList(qpcr.cases), glm(HLABW4 ~ kir, family='binomial')))
fit0 <- as.mira(with(imputationList(qpcr.cases), glm(HLABW4 ~ 1, family='binomial')))
round(mice::pool.compare(fit1=fit1, fit0=fit0)$pvalue, digits=4)
# table3.a SNP
with(all.comb.cases, table(HLABW4,kir))
sum(with(all.comb.cases, table(HLABW4,kir)))
fit1 <- as.mira(with(imputationList(all.cases), glm(HLABW4 ~ kir, family='binomial')))
fit0 <- as.mira(with(imputationList(all.cases), glm(HLABW4 ~ 1, family='binomial')))
round(mice::pool.compare(fit1=fit1, fit0=fit0)$pvalue, digits=4)

# table3.b qPCR
with(qpcr.comb.cases, table(HLABW4,kir3dl1))
sum(with(qpcr.comb.cases, table(HLABW4,kir3dl1)))
fit1 <- as.mira(with(imputationList(qpcr.cases), glm(HLABW4 ~ kir3dl1, family='binomial')))
fit0 <- as.mira(with(imputationList(qpcr.cases), glm(HLABW4 ~ 1, family='binomial')))
round(mice::pool.compare(fit1=fit1, fit0=fit0)$pvalue, digits=4)
# table3.b SNP
with(all.comb.cases, table(HLABW4,kir3dl1))
sum(with(all.comb.cases, table(HLABW4,kir3dl1)))
fit1 <- as.mira(with(imputationList(all.cases), glm(HLABW4 ~ kir3dl1, family='binomial')))
fit0 <- as.mira(with(imputationList(all.cases), glm(HLABW4 ~ 1, family='binomial')))
round(mice::pool.compare(fit1=fit1, fit0=fit0)$pvalue, digits=4)

# table3.c qPCR
with(qpcr.comb.cases, table(HLABW4I,kir3ds1))
sum(with(qpcr.comb.cases, table(HLABW4I,kir3ds1)))
fit1 <- as.mira(with(imputationList(qpcr.cases), glm(HLABW4I ~ kir3ds1, family='binomial')))
fit0 <- as.mira(with(imputationList(qpcr.cases), glm(HLABW4I ~ 1, family='binomial')))
round(mice::pool.compare(fit1=fit1, fit0=fit0)$pvalue, digits=4)
# table3.c SNP
with(comb.cases, table(HLABW4I,kir3ds1))
sum(with(comb.cases, table(HLABW4I,kir3ds1)))
fit1 <- as.mira(with(imputationList(cases), glm(HLABW4I ~ kir3ds1, family='binomial')))
fit0 <- as.mira(with(imputationList(cases), glm(HLABW4I ~ 1, family='binomial')))
round(mice::pool.compare(fit1=fit1, fit0=fit0)$pvalue, digits=4)


#lapply(imputed.datasets, function(d) table(d[grep('4I',d$HLA_Bw),c('t1d','knn.call')]))
#lapply(loo.knn, function(d) table(d[grep('4I',d$HLA_Bw),c('t1d','y')]))

#samples which are consistently misclassified across all 10 multiply imputed datasets
misclass.sampleid <- do.call( 'intersection', lapply(1:10, function(i) row.names(loo.knn[[i]])[which(with(loo.knn[[i]], y!=y.knn))]) )
misclass.dil_subjectid <- lookup[lookup$sampleid %in% misclass.sampleid,'dil_subjectid']

# how good are we at classifying the different cn
lapply(loo.knn, function(x) with(x, table(loo.KIR,KIR)))
lapply(loo.knn, function(x) {
       m <- with(x, table(loo.KIR,KIR))
       100*round(t(t(m)/colSums(m)),digits=2)
       })


# qpcr.data
source('Code/plate.R') 
d <- zero.cn.trans(load.experiment(experiment.file='Data/qPCR/processed/experiment.csv', QC=qc.experiment, trans=function(x)-x))
misclass.d <- d[d$sample %in% misclass.dil_subjectid,]

plot(d$DS1.peaknorm.median, d$DL1.peaknorm.median, pch=20)
points(misclass.d$DS1.peaknorm.median, misclass.d$DL1.peaknorm.median, pch=20, col='red')


sapply(imputed.datasets, function(x) with(x, table(loo.KIR!=KIR)))
sapply(loo.knn, function(x) with(x, table(loo.KIR!=KIR)))

lapply(loo.knn, function(x) with(x, table(loo.KIR,KIR)))

misclass.dil_subjectid <- do.call('intersection', sapply(imputed.datasets, function(x) x[which(with(x, loo.KIR!=KIR)),'dil_subjectid']))
misclass.sample <- do.call('intersection', sapply(imputed.datasets, function(x) x[which(with(x, loo.KIR!=KIR)),'sample']))
#of these, which ones are consistently wrong

#LD errors
LD.errors <- data.frame(
KIR=qpcr.posteriors[misclass.sample,'qpcr.cn'],
KIR.post=round(apply( qpcr.posteriors[misclass.sample,as.character(GENO)],1,max ),digits=2),
loo.KIR=apply( sapply(imputed.datasets, function(x) x[misclass.sample,'loo.KIR']), 1, function(x) names(table(x)) )
) 
LD.errors <- LD.errors[order(LD.errors$KIR, LD.errors$loo.KIR),]
rownames(LD.errors) <- NULL
xtable(LD.errors)


lapply(imputed.datasets, function(x) x[which(with(x, loo.KIR!=KIR)),c('sample','KIR')])


X <- as.data.frame(t(sapply(apply(sapply(imputed.datasets, function(x) x$KIR),1,function(y) table(y)),
            function(x) {
                l <- numeric(length(geno))
                names(l) <- geno
                l[geno] <- 0
                l[names(x)]  <- x/10
                return(l)
            })))
imp.data <- cbind(imp.data, X)
Y <- data.frame(dil_subjectid=misclass.dil_subjectid,
           snp.kir=apply(imp.data[imp.data$dil_subjectid %in% misclass.dil_subjectid,geno],1,function(x) names(which.max(x))),
           snp.post=apply(imp.data[imp.data$dil_subjectid %in% misclass.dil_subjectid,geno],1,max),
           qpcr.kir=apply(qpcr.posteriors[qpcr.posteriors$dil_subjectid %in% misclass.dil_subjectid,geno],1,function(x) names(which.max(x))),
           qpcr.post=apply(qpcr.posteriors[qpcr.posteriors$dil_subjectid %in% misclass.dil_subjectid,geno],1,max))

print(LD.errors <- Y[Y$qpcr.kir != Y$snp.kir & Y$qpcr.post > .99 & Y$snp.post ==1,])


misclass.post <- data.frame( qpcr.cn=qpcr.posteriors[qpcr.posteriors$dil_subjectid %in% misclass.dil_subjectid,'qpcr.cn'], knn.cn=apply( sapply(loo.knn, function(x) x[x$dil_subjectid %in% misclass.dil_subjectid,'y.knn']), 1, unique ), max.post=apply(qpcr.posteriors[qpcr.posteriors$dil_subjectid %in% misclass.dil_subjectid,GENO.COL$geno], 1, max) )

apply(qpcr.posteriors[,GENO.COL$geno], 2, quantile)

#uncertainty per CN group
sapply( GENO.COL$geno, function(g) length(which(qpcr.posteriors[,g]>.95))/table(qpcr.posteriors$qpcr.cn)[[g]])

# 99% of all samples are classified with 99% certainty
prop.table(table(apply(qpcr.posteriors[,GENO.COL$geno],1,max)>.99))



### Figure-2
head(snp.data <- read.csv('Data/SNP/rs592645.csv'))
dim(snp.data <- merge(lookup, snp.data, by.x='sampleid', by.y='sample', all.y=TRUE))

head(qpcr.data <- d[,c('sample','DS1.peaknorm.median','DL1.peaknorm.median')])
colnames(qpcr.data) <- c('dil_subjectid','DS1.peaknorm.median','DL1.peaknorm.median')

dim(snp.data2 <- merge(snp.data, qpcr.data, by.x='dil_subjectid', by.y='dil_subjectid'))

dim(misclass.data <- snp.data2[snp.data2$sampleid %in% misclass.sampleid,])

dim(misclass.data <- merge(misclass.data, qpcr.posteriors))


misclass.data <- merge(misclass.data, data.frame(sampleid=misclass.sampleid, knn.cn=apply( sapply(1:10, function(i) loo.knn[[i]][misclass.sampleid,'y.knn']), 1, unique)))

x <- misclass.data[order(misclass.data$sampleid),c('sampleid','qpcr.cn','knn.cn')]
x <- x[!duplicated(x),]
x[order(x$dist.cn,decreasing=T),]

cn.dist <- function(cn1, cn2) as.numeric(dist(rbind(as.numeric(unlist(strsplit(cn1, '-'))),as.numeric(unlist(strsplit(cn2, '-'))))))

x$dist.cn <- as.numeric(apply(x[,c('qpcr.cn','knn.cn')], 1, function(a) cn.dist(a[1],a[2])))



GENO.COL <- data.frame( geno=c("?-?", "0-2", "1-1", "2-0", "2-1", "0-1", "1-2", "3-0", "1-0"), col=c('black', 'darkgreen', 'pink', 'blue', 'lightblue', 'green', 'purple', 'grey', 'red') )
g<-as.vector(GENO.COL[,2])
names(g)<-GENO.COL[,1]
pdf('~nikolas/stats-archive/Papers/KIR/bmc_article/figures/Figure-2.pdf',height=6,width=9)
snp.data <- snp.data[ order(snp.data$Copies,decreasing=TRUE), ]
snp.data$Copies <- as.character(snp.data$Copies)
snp.data$Copies[ is.na(snp.data$Copies) ] <- "?-?"
snp.data$Copies <- factor(snp.data$Copies,levels=c("?-?","0-1","0-2","1-2","1-1","2-1","1-0","2-0","3-0"))
ggplot(snp.data, aes(x=theta,y=R,col=Copies)) +
  geom_point(data=subset(snp.data, is.na(KIR)),size=1.3) +
  geom_point(data=subset(snp.data, !is.na(KIR)),size=1.3) +
  facet_grid(qPCR~cc, margins=TRUE) + theme_bw() +
  scale_colour_manual(values=g, guide=guide_legend(override.aes=list(size=3),label.hjust=.5, label.vjust=.5, label.theme=element_text(angle=-45), label.position='top', title='KIR3DS1-KIR3DL1 Copy Number', title.position='right', title.vjust=.2, title.theme=element_text(face='italic',angle=0))) + theme(legend.position='top') +
  xlab(expression(theta)) + ylab('R')
dev.off()

ggplot(snp.data, aes(x=theta,y=R,col=Copies)) + geom_point(size=1) + geom_point(data=subset(snp.data,sample%in%misclass.sampleid), size=2) +
facet_wrap( ~ Copies) + theme_bw() +
scale_colour_manual(values=g, guide=guide_legend(override.aes=list(size=3),label.hjust=.5, label.vjust=.5, label.theme=element_text(angle=-45), label.position='top', title='KIR3DS1-KIR3DL1 Copy Number', title.position='right', title.vjust=.2, title.theme=element_text(face='italic',angle=0))) + theme(legend.position='top') +
xlab(expression(theta)) + ylab('R')

quartz()
ggplot(data=snp.data2, aes(x=DS1.peaknorm.median,y=DL1.peaknorm.median,col=Copies)) + geom_point(size=1) + geom_point(data=subset(snp.data2,sampleid%in%misclass.sampleid), size=2) + facet_wrap( ~ Copies) + theme_bw() +
scale_colour_manual(values=g, guide=guide_legend(override.aes=list(size=3),label.hjust=.5, label.vjust=.5, label.theme=element_text(angle=-45), label.position='top', title='KIR3DS1-KIR3DL1 Copy Number', title.position='right', title.vjust=.2, title.theme=element_text(face='italic',angle=0))) + theme(legend.position='top')


