#rm(box)
#unloadNamespace('flowClust')
#association tests
library(nlme)
library(lme4)
library(influence.ME)
library(beeswarm)
library(whisker)

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
.Names = c("individual", 
"pch", "day1", "day2", "col", "day.diff", "t1d"), row.names = c("CB00165D", 
"CB00366X", "CB00396E", "CB00406Q", "CB01484M", "CB01494Y", "CB01495Z", 
"CB01498C", "CB01503H", "CB01504J"), class = "data.frame")
REPEATS <- rbind(cbind(REPEATS,date=REPEATS$day1),cbind(REPEATS,date=REPEATS$day2))

# add metadata such as individual, age, sex, disease status
metadata <- function(pstat5.pheno) {
    pstat5.pheno$individual <- as.character( pstat5.pheno$individual )
    pstat5.pheno[which(pstat5.pheno$individual=='KM00685J'),'individual'] <- 'KM00685T'
    pstat5.pheno[which(pstat5.pheno$individual=='KM00746G'),'individual'] <- 'KM00746K'
    pstat5.pheno[which(pstat5.pheno$individual=='KM00861S'),'individual'] <- 'KM00861K'
    # dgap individuals:
    dgap <- read.csv('~nikolas/Projects/IL2/IL2RA-genotype/dgap-IL2RA.csv',stringsAsFactors=FALSE)
    dgap$individual <- dgap$NikoIDcorrected
    # cb individuals:
    ip <- read.csv('~nikolas/Projects/IL2/IL2RA-genotype/ip-IL2RA-2013-12-23.csv')
    ip$individual <- ip$uniqueID
    #
    ip <- rbind(dgap[,c('individual','t1d', 'sex')], ip[,c('individual','t1d', 'sex')])
    dim(pstat5.pheno <- merge(ip, pstat5.pheno))
    pstat5.pheno$t1d <- ifelse(pstat5.pheno$t1d==1, 'control', 'case')
    pstat5.pheno$t1d <- factor(pstat5.pheno$t1d, level=c('case','control'))
    pstat5.pheno$sex <- factor(ifelse(pstat5.pheno$sex==1,'M','F'), level=c('M','F'))
    dim(pstat5.pheno <- merge(pstat5.pheno, REPEATS, by=c('individual','date','t1d'), all.x=TRUE))
    pstat5.pheno$pch[which(is.na(pstat5.pheno$pch))] <- 'x'
    #d <- subset(pstat5.pheno,cell.type=='Memory Eff')
    #library(nlme)
    #summary(m <- lme(auc.pstat5 ~ t1d, random=~1|individual/date, data=d))
    #summary(m <- lme(auc.pstat5 ~ t1d, random=~1|date/individual, data=d))
    #library(lme4)
    #summary(m<-lmer(auc.pstat5 ~ t1d + (1|individual) + (1|date), data=d))
    #plot(as.Date(d$date), d$auc.pstat5, pch=20, col=ifelse(d$t1d=='case','red','black'), xlab='' , ylab='Memory Eff pSTAT5 AUC')
    #d$date <- as.Date(d$date)
    #by(d, d$individual, function(x) { if (nrow(x)==2) segments(x[1,'date'],x[1,'auc.pstat5'],x[2,'date'],x[2,'auc.pstat5'],lwd=.5,col=unique(as.factor(x$individual)))})
    #abline(v=as.Date(d$date),col='gray',lwd=.25,lty=2)
    #do.call('rbind', by(d, d$date, function(x) table(x$t1d)))
    return(pstat5.pheno)
}


format.pval.col <- function(pval) {
    if (as.numeric(pval)<.05)
        sprintf('\\textcolor{red}{%s}',format.pval(pval))
    else
        format.pval(pval)
}


f <- function(x, pheno, covariate, ...) {
    i <- 1
    fm <- formula(paste('pheno ~',covariate,sep=''))
    fm2 <- formula(paste('pheno ~',covariate,'+ (1|individual) + (1|date)',sep=''))
    print(dim(X <- na.omit(x[,c(pheno,covariate,'individual','pch','date')])))
    X$date <- as.factor(X$date)
    X$individual <- as.factor(X$individual)
    print(head(X))
    colnames(X) <- c('pheno',covariate,'individual','pch','date')
    m1<-lme(fm, random=list(~1|individual,~1|date), data=X, method='ML')
    e <- try(intervals(m1), TRUE)
    if (class(e)=='try-error') {
        cat('Using REML instead\n')
        m1<-lme(fm, random=list(~1|individual,~1|date), data=X, method='REML')
    }
    m2<-lmer(fm2, data=X)
    m <- list(
            m=m1,
            m2=m2,
            cd=cooks.distance(influence(m2,obs=TRUE)),
            #lower=round(as.numeric(intervals(m1)$fixed[,1][2]),digits=3),
            #intercept=round(as.numeric(intervals(m1)$fixed[,2][1]),digits=3),
            #fixed=round(as.numeric(intervals(m1)$fixed[,2][2]),digits=3),
            #upper=round(as.numeric(intervals(m1)$fixed[,3][2]),digits=3),
            pvalue=format.pval(as.numeric(summary(m1)$tTable[,'p-value'][[2]]))
        )
    #par(oma = c(4, .5, 1, 1), mfrow=c(1,1))
    #par(mfrow=c(1,1))
    #if (covariate=='age') plot(fm, data=X, col=X$col, pch=X$pch, ...)
    #else beeswarm(fm, data=X, pwcol=X$col, pwpch=X$pch, cex=1.5, corral='wrap', ...)
    #else plot(as.numeric(X$t1d), X[,'pheno'])
    stripchart(list('Healthy control'=X[which(X$t1d=='control'),'pheno'],'T1D'=X[which(X$t1d=='case'),'pheno']), vertical=T, col=c('black','red'),  ...)
    by(X, X$date, function(x) {
       if (nrow(x)>1) {
           controls <- x[which(x$t1d=='control'),'pheno']
           cases <- x[which(x$t1d=='case'),'pheno']
           if (length(controls)>0 && length(cases)>0) segments(1,controls,2,cases,lwd=.25)
       }
    })
    #for (m1 in m) {
        #abline(b=m1$fixed, a=m1$intercept, col=cols[[i]], lwd=1.5)
        #i <- i+1
    #}
    #par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new=TRUE)
    mtext(paste('p =',as.character(m[['pvalue']])))
    m
} 


plot.t1d.assoc <- function(pstat5.pheno, pdf.file,ylab) {
    pdf(pdf.file,width=10,height=10)
    par(mfrow=c(2,2))
    m <- list()
    CELL.TYPES <- c("Memory Eff", "Memory Treg", "Naive Eff", "Naive Treg")
    #pSTAT5 pheno tested at different doses
    celltype.dose <- list('Memory Eff'='10U', 'Memory Treg'='01U', 'Naive Eff'='1000U', 'Naive Treg'='01U')
    for (cell.type in CELL.TYPES) {
        print(cell.type)
        print(dose <- celltype.dose[[cell.type]])
        d <- pstat5.pheno[which(pstat5.pheno$cell.type==cell.type),]
        m[[cell.type]] <- f(x=d, pheno=dose,covariate='t1d',main=cell.type,ylab=paste(ylab,dose))
    }
    dev.off()
}



# pSTAT5 MFI
print(load('/dunwich/scratch/nikolas/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/pstat5.mfi.RData'))
pstat5.pheno <- metadata(pstat5.mfi)
plot.t1d.assoc(pstat5.pheno, ylab='pSTAT5 MFI', pdf.file='~nikolas/Thesis/figures/pstat5-mfi-t1d-celltypes.pdf')

# NN pSTAT5 MFI
print(load('/dunwich/scratch/nikolas/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/nn.pstat5.mfi.RData'))
pstat5.pheno <- metadata(pstat5.mfi)
colnames(pstat5.pheno) <- gsub('X','',colnames(pstat5.pheno))
plot.t1d.assoc(pstat5.pheno, ylab='pSTAT5 MFI', pdf.file='~nikolas/Thesis/figures/nn-pstat5-mfi-t1d-celltypes.pdf')

# % pSTAT5+
print(load('/dunwich/scratch/nikolas/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/pstat5.pos.RData'))
pstat5.pheno <- metadata(pstat5.pos)
plot.t1d.assoc(pstat5.pheno, ylab='% pSTAT5+', pdf.file='~nikolas/Thesis/figures/pstat5-pos-t1d-celltypes.pdf')

# NN % pSTAT5+
print(load('/dunwich/scratch/nikolas/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/nn.pstat5.pos.RData'))
pstat5.pheno <- metadata(pstat5.pos)
plot.t1d.assoc(pstat5.pheno, ylab='NN % pSTAT5+', pdf.file='~nikolas/Thesis/figures/nn-pstat5-pos-t1d-celltypes.pdf')

# base NN % pSTAT5+
print(load('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/nn.base.pstat5.pos.RData'))
pstat5.pheno <- metadata(pstat5.pos)
plot.t1d.assoc(pstat5.pheno, ylab='base NN % pSTAT5+', pdf.file='~nikolas/Thesis/figures/nn-base-pstat5-pos-t1d-celltypes.pdf')





cases <- subset(d,t1d=='case') 
controls <- subset(d,t1d=='control')

ladderplot( cbind(cases[order(cases$date)[1:50],'auc.norm.pstat5'],controls[order(controls$date)[1:50],'auc.norm.pstat5']) )
matplot( cbind(cases[order(cases$date)[1:50],'auc.norm.pstat5'],controls[order(controls$date)[1:50],'auc.norm.pstat5']) )

d$t1d <- relevel(d$t1d, 'control')
cases[order(cases$date)[1:50],'auc.norm.pstat5'],controls[order(controls$date)[1:50],'auc.norm.pstat5']

print(summary(lmer( t1d ~ auc.norm.pstat5 + (1|as.factor(uniqueID)) + (1|as.factor(cell.type)), data=pstat5.mfi )))


# Hui's analysis

#d = read.csv('/chiswick/data/hui/hui_archive/flow/FACS/Tony_ls_dgap_pairs/LS_DGAP_paired_data_CW_copy_tony_csv', header=T, as.is=T)
d = read.csv('~/dunwich/Projects/IL2/LS_DGAP_paired_data_CW_copy_tony_csv', header=T, as.is=T)

d = d[1:51, -c(5,6)]
d$HC[7] = 'pair 40'
colnames(d) = c('pairs', 'naive', 'memory', 'actv', 'naive.1', 'memory.1', 'actv.1')

## only for the paired data, do paired t-test
d.p = na.omit(d)
d.p$dif.naive = with(d.p, naive.1-naive)
d.p$dif.memory = with(d.p, memory.1 - memory)
d.p$dif.actv = with(d.p, actv.1 - actv) 

plot(d.p$dif.naive)
plot(d.p$dif.memory)
plot(d.p$dif.actv)

t.test(d.p$dif.naive, alternative = "two.sided")
# t = 1.6879, df = 39, p-value = 0.09941

t.test(d.p$dif.memory, alternative = "two.sided")
# t = -0.0077, df = 39, p-value = 0.9939

t.test(d.p$dif.actv, alternative = "two.sided")
# t = -0.7752, df = 39, p-value = 0.4429


## assign Row i with empty pair info with the pair info of Row i-1
## for each trio group, if phenotype measured from 2 IDs, use the average to form 1 pair
  
n = dim(d)[1]
for (i in 1:n) {
  d$pairs[i] = ifelse(d$pairs[i] != '', d$pairs[i], d$pairs[i-1])
}  

d.dup = d[c(which(duplicated(d$pairs)), which(duplicated(d$pairs))-1), ]

d.dup.split = split(d.dup, d.dup$pairs)

d.dup.avrg = vector('list', length(d.dup.split))
for (j in names(d.dup.split)) {
  d.dup.avrg[[j]]$pairs = names(d.dup.split)[j]

  d.dup.avrg[[j]]$naive = mean(na.omit(d.dup.split[[j]]$naive))
  d.dup.avrg[[j]]$memory = mean(na.omit(d.dup.split[[j]]$memory)) 
  d.dup.avrg[[j]]$actv = mean(na.omit(d.dup.split[[j]]$actv)) 

  d.dup.avrg[[j]]$naive.1 = mean(na.omit(d.dup.split[[j]]$naive.1))
  d.dup.avrg[[j]]$memory.1 = mean(na.omit(d.dup.split[[j]]$memory.1)) 
  d.dup.avrg[[j]]$actv.1 = mean(na.omit(d.dup.split[[j]]$actv.1)) 
}


## combine with orignal paired data, do paired t-test  
d.dup.avrg = do.call('rbind.data.frame', d.dup.avrg)
d.dup.avrg$pairs = rownames(d.dup.avrg)
d.dup.avrg$dif.naive = with(d.dup.avrg, naive.1-naive)
d.dup.avrg$dif.memory = with(d.dup.avrg, memory.1 - memory)
d.dup.avrg$dif.actv = with(d.dup.avrg, actv.1 - actv) 

d.p.dup.avrg = rbind(subset(d.p, !pairs %in% d.dup.avrg$pairs), d.dup.avrg)

t.test(d.p.dup.avrg$dif.naive, alternative = "two.sided")
# t = 2.197, df = 39, p-value = 0.03403

t.test(d.p.dup.avrg$dif.memory, alternative = "two.sided")
# t = 0.4235, df = 39, p-value = 0.6743

t.test(d.p.dup.avrg$dif.actv, alternative = "two.sided")
# t = -0.5241, df = 39, p-value = 0.6032


## for each trio group, use data of 1 ID twice to form two pairs

d.dup.pair = vector('list', length(d.dup.split))
for (j in names(d.dup.split)) {
  d.dup.pair[[j]] = d.dup.split[[j]]
  if (is.na(d.dup.split[[j]]$naive[1])) {
    d.dup.pair[[j]][1, 2:4] = d.dup.pair[[j]][2, 2:4]
  }
  if (is.na(d.dup.split[[j]]$naive.1[1])) {
    d.dup.pair[[j]][1, 5:7] = d.dup.pair[[j]][2, 5:7]
  }
}

d.dup.pair = do.call('rbind', d.dup.pair)
d.dup.pair$dif.naive = with(d.dup.pair, naive.1-naive)
d.dup.pair$dif.memory = with(d.dup.pair, memory.1 - memory)
d.dup.pair$dif.actv = with(d.dup.pair, actv.1 - actv)

d.dup.pair.comb = rbind(subset(d.p, !pairs %in% d.dup.avrg$pairs), d.dup.pair) 
d.dup.pair.comb$pairs = as.numeric(gsub('pair ', '', d.dup.pair.comb$pairs))

 



## bootstapping in 2 strata (pairs, trios), use weight 3/4 for trios

library(boot)
set.seed(10000)

d.dup.pair.comb$group = c(rep(1, 29), rep(2, 22))

d.pair = data.frame(d.dup.pair.comb[1:29, -c(2:7,11)], dif.naive.1 = NA, dif.memory.1 = NA, dif.actv.1 = NA, group = d.dup.pair.comb[1:29, 11])

d.trio = d.dup.pair.comb[-(1:29), ]
trio.for.boot = cbind(d.trio[seq(1,dim(d.trio)[1],by=2), c(1,8:10)],d.trio[seq(2,dim(d.trio)[1],by=2), 8:11])
colnames(trio.for.boot)[5:7] = colnames(d.pair)[5:7]

d.boot = rbind(d.pair, trio.for.boot)

w = c(rep(1, 29), rep(3/4, 11))/sum(c(rep(1, 29), rep(3/4, 11)))
grp1 = 1:table(d.boot$group)[1]


avrg = function(d, weight=w) {
     mn.n <- sum(d$dif.naive * weight, d$dif.naive.1[-grp1] * weight[-grp1])  
     mn.m <- sum(d$dif.memory * weight, d$dif.memory.1[-grp1] * weight[-grp1]) 
     mn.a <- sum(d$dif.actv * weight, d$dif.actv.1[-grp1] * weight[-grp1])
     #c(mn.n, sd.n, mn.m)
     c(mn.n, mn.a, mn.m)
}

bstp = boot(data=d.boot, statistic=avrg, R=100000, weights = w, stype = "w", strata = d.boot$group)
m = bstp$t0
s = sqrt(diag(var(bstp$t)))

# naive
t.n = m[1]/s[1]
pval.n = 2*pt(abs(t.n), 50, lower.tail=FALSE)
# t = 1.586395, df = 50, p-value = 0.1189545
 
# memory
t.m = m[2]/s[2]
pval.m = 2*pt(abs(t.m), 50, lower.tail=FALSE)
# t = 1.577787, df = 50, p-value = 0.1209212

# activated
t.a = m[3]/s[3]
pval.a = 2*pt(abs(t.a), 50, lower.tail=FALSE)
# t = -0.9558598, df = 50, p-value = 0.343742



