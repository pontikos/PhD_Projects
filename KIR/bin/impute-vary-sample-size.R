setwd('/chiswick/data/stats-archive/Papers/KIR/')

library(ggplot2)

source("impute_functions.R")

## qPCR calls
print(load("Data/qPCR/MI-10-datasets.RData"))

## add ichip sample id
lookup <- read.table("Data/SNP/KIR-lookup-2013-05-07.tab",sep="\t",header=TRUE,as.is=TRUE)
lookup.clean <- read.csv("Data/ic-cc-clean-lookup-2013-02-27.csv")
lookup <- subset(lookup, lookup$immunochipsample %in% lookup.clean$sampleid)

geno2 <- lapply(imputations, function(geno) {
  tmp <- merge(lookup, geno, by.y="dil_subject", by.x="sample", all.x=TRUE)
  tmp <- subset(tmp,!is.na(immunochipsample))
  tmp <- tmp[!duplicated(tmp$sample),]
  rownames(tmp) <- tmp$immunochipsample
  tmp
})

str(imputations[[1]])
str(geno2[[1]])
geno <- geno2

################################################################################
 ## prep data

## Illumina files for local/predictive SNPs
cs.cols <- c( "Sample_ID", "SNP_Name", "Chr", "Position", "GC_Score",
"Genotype_Forward", "Theta", "R", "X", "Y", "B_Allele_Freq",
"Log_R_Ratio", "GT_Score")

ctr.cols <- c("Sample_ID", "SNP_Name", "Chr", "Position",
"Genotype_Forward", "Theta", "R", "X", "Y", "B_Allele_Freq",
"Log_R_Ratio", "GC_Score", "GT_Score")

case.files <- list.files("lrr-data", full=TRUE, pattern="case")
control.files <- list.files("lrr-data", full=TRUE, pattern="control")
case.data <- lapply(case.files, read.snp.mi, cs.cols)
control.data <- lapply(control.files, read.snp.mi, ctr.cols)
names(case.data) <- sub("lrr-data/cases-","case-",case.files)
names(control.data) <- sub("lrr-data/controls-","control-",control.files)
data <- rbind(do.call("rbind",case.data),
              do.call("rbind",control.data)[,colnames(case.data[[1]])])
data$sample <- sub(".*\\.","",rownames(data))

################################################################################

## summary(m.3ds1 <- glm(KIR.3DS1 ~ (Theta+Log_R_Ratio+GC_Score):SNP_Name, data=data))
## summary(m.3dl1 <- glm(KIR.3DL1 ~ (Theta+Log_R_Ratio+GC_Score):SNP_Name, data=data))

##  glm - what is uninformative?
## seq-rs1654644
## rs581623
## seq-rs55761930
## seq-rs592645
## seq-t1d-19-60034052-C-T
## seq-t1d-19-60056721-C-T

snps.drop <- c("rs10422740","rs3826878","rs62122181","rs12461010","seq-t1d-19-60056605-A-T")
data <- subset(data, !(SNP_Name %in% snps.drop))

################################################################################
## rs10422740              0.7241379 0.2413793 0.7241379 0.2413793
## rs10423751              0.4482759 0.5172414 0.2413793 0.7241379
## rs2569678               0.7241379 0.3103448 0.6896552 0.3448276
## seq-rs11669355          0.7241379 0.3793103 0.6206897 0.3448276
## seq-rs12982080          0.5517241 0.5517241 0.2413793 0.7241379
## seq-rs270775            0.3103448 0.7241379 0.3793103 0.6551724
## seq-rs34078188          0.6896552 0.2758621 0.7241379 0.2413793
## seq-rs604999            0.6551724 0.3103448 0.7241379 0.2413793
## seq-t1d-19-60054973-T-C 0.7241379 0.2413793 0.6206897 0.3448276
## seq-rs2295805           0.6896552 0.2758621 0.7586207 0.2068966
## seq-rs62122181          0.7586207 0.2068966 0.7586207 0.2068966
## seq-rs62123050          0.7586207 0.3448276 0.3103448 0.6551724
## 1kg_19_18101748         0.3448276 0.7586207 0.3103448 0.7931034
## rs3865507               0.9655172 0.1379310 0.7241379 0.2413793
## seq-rs592645            0.9655172 0.1379310 0.1034483 0.8620690


# X,Y,R not useful
#good.snps <- c("seq-rs12976350","seq-rs592645","seq-rs4806568","seq-t1d-19-60034052-C-T","seq-rs674268","seq-rs3865510","rs12461010")

################################################################################

## make data for prediction

## explore different variables for prediction
table(cc <- grepl("cases",data$f))
## datasets <- list(XY.case = make.wide(data[cc,], cols=c("X","Y")),
##                  RT.case = make.wide(data[cc,], cols=c("R","Theta")),
##                  LB.case = make.wide(data[cc,], cols=c("B_Allele_Freq","Log_R_Ratio")),
##                  XY.ctrl = make.wide(data[!cc,], cols=c("X","Y")),
##                  RT.ctrl = make.wide(data[!cc,], cols=c("R","Theta")),
##                  LB.ctrl = make.wide(data[!cc,], cols=c("B_Allele_Freq","Log_R_Ratio")))
## datasets <- list(case=make.wide(data[cc , ],
##                    cols=c("Theta","Log_R_Ratio", "GT_Score","GC_Score"),
##                    mi=TRUE),
##                  ctrl=make.wide(data[!cc, ],
##                    cols=c("Theta","Log_R_Ratio", "GT_Score","GC_Score"),
##                    mi=TRUE),
##                  both=make.wide(data,
##                    cols=c("Theta","Log_R_Ratio", "GT_Score","GC_Score"),
##                    mi=TRUE))
datasets <- list(case=make.wide(data[cc , ],
                   cols=c("Theta","R","Log_R_Ratio", "GT_Score","GC_Score"),
                   mi=TRUE),
                 ctrl=make.wide(data[!cc, ],
                   cols=c("Theta","R","Log_R_Ratio", "GT_Score","GC_Score"),
                   mi=TRUE),
                 both=make.wide(data,
                   cols=c("Theta","R","Log_R_Ratio", "GT_Score","GC_Score"),
                   mi=TRUE))

## drop constant variables
v <- lapply(datasets, function(x) apply(x[-1], 2, var, na.rm=TRUE))
to.drop <- which(v[[3]]==0)
datasets <- lapply(datasets, function(d) d[, setdiff(colnames(d), names(to.drop))])

## drop one of any pair of highly correlated variables
cr <- lapply(datasets, function(x) cor(x[,-1], use="pair"))
wh <- which(abs(cr[[3]])>0.98, arr.ind=TRUE)
wh <- wh[wh[,1]<wh[,2],]
(to.drop <- cbind(wh,cr[[3]][wh])) # don't need theta & BAF. 
##                    row col           
## Theta.seq-rs648689  63  67  0.9949678
## Theta.seq-rs598452  49  71 -0.9860313
## Theta.seq-rs597598  45  75 -0.9906294

datasets <- lapply(datasets, function(d) d[, setdiff(colnames(d), rownames(to.drop))])
lapply(datasets, function(d) table(complete.cases(d)))

################################################################################

## what is best single SNP?
d <- cbind(datasets$both, KIR=geno[[1]][as.character(datasets$both$sample),"geno"])
#d <- subset(d,!is.na(KIR))

d$k1 <- as.numeric(substr(as.character(d$KIR),3,3))
d$k2 <- as.numeric(substr(as.character(d$KIR),4,4))
d$cc <- sub("-.*","",rownames(d))
d$n<- sapply(lapply(strsplit(sub("cn","",d$KIR),"-"),as.numeric),sum)

summary(m.n <- glm(n ~ ., data=d[,c("n",grep("Log_R",colnames(d),value=TRUE))]))
snps.lrr <- c("rs581623","rs10407958","rs3865510","rs604999","seq.t1d.19.60034052.C.T")
summary(m.R <- glm(n ~ ., data=d[,c("n",grep("^R",colnames(d),value=TRUE))]))


summary(m.1 <- glm(k1 ~ ., data=d[,c("k1",grep("Theta|Log_R",colnames(d),value=TRUE))]))
summary(m.theta <- glm(k1 ~ ., data=d[,c("k1",grep("Theta",colnames(d),value=TRUE))]))
## Theta.seq-rs592645                   1.096674   0.126509   8.669  < 2e-16 ***
## `Theta.seq-rs4806568`            1.23140    0.15057   8.178 6.22e-16 ***
## `Theta.seq-rs604999`            -0.54581    0.07344  -7.432 1.82e-13 ***
summary(m.1 <- glm(k1 ~ ., data=d[,c("k1",grep("Log_R",colnames(d),value=TRUE))]))
## `Log_R_Ratio.seq-rs4806568`           -0.84352    0.06709 -12.574  < 2e-16 ***
## `Log_R_Ratio.seq-rs62122181`           0.75434    0.04099  18.403  < 2e-16 ***
## `Log_R_Ratio.seq-rs674268`             0.17392    0.01630  10.673  < 2e-16 ***
## `Log_R_Ratio.seq-t1d-19-60056721-C-T` -1.46477    0.10035 -14.597  < 2e-16 ***


## Log_R_Ratio.rs581623                   0.183401   0.037259   4.922 9.55e-07 ***
## Theta.seq-rs10500318                 0.647370   0.114375   5.660 1.83e-08 ***
## Theta.seq-rs3865510                 0.674448   0.101103   6.671 3.63e-11 ***

summary(m.2 <- glm(k2 ~ ., data=d[,c("k2",grep("Theta|Log_R",colnames(d),value=TRUE))]))
summary(m.2 <- glm(k2 ~ ., data=d[,c("k2",grep("Theta",colnames(d),value=TRUE))]))
## Theta.rs4806585                  0.56055    0.08674   6.462 1.40e-10 ***
## `Theta.seq-rs10500318`          -0.79358    0.09870  -8.040 1.85e-15 ***
## `Theta.seq-rs1654644`           -0.14540    0.01803  -8.063 1.54e-15 ***
## `Theta.seq-rs2295805`           -0.76907    0.14314  -5.373 9.02e-08 ***
## `Theta.seq-rs3865510`           -0.80352    0.09308  -8.633  < 2e-16 ***
## `Theta.seq-rs4806568`            1.19658    0.13113   9.125  < 2e-16 ***
## `Theta.seq-rs604999`            -0.48908    0.06396  -7.647 3.73e-14 ***
## `Theta.seq-rs62122181`          -0.73499    0.07798  -9.425  < 2e-16 ***
summary(m.2 <- glm(k2 ~ ., data=d[,c("k2",grep("Log_R",colnames(d),value=TRUE))]))
## `Log_R_Ratio.seq-rs4806568`            0.74327    0.05705  13.028  < 2e-16 ***
## `Log_R_Ratio.seq-rs62122181`          -0.75438    0.03486 -21.640  < 2e-16 ***
## `Log_R_Ratio.seq-rs674268`            -0.18051    0.01386 -13.026  < 2e-16 ***
## `Log_R_Ratio.seq-t1d-19-60056721-C-T`  1.24747    0.08534  14.618  < 2e-16 ***

## `Theta.seq-rs592645`                  -0.782991   0.098425  -7.955 3.61e-15 ***
## `Theta.seq-rs649216`                  -0.504138   0.066929  -7.532 8.80e-14 ***
## `Log_R_Ratio.seq-rs3865510`            0.246832   0.041358   5.968 3.02e-09 ***

snps.r <- make.names(c("rs581623","seq-rs10407958","seq-rs10500318","seq-rs2295805","rs3865510",
                       "seq-rs592645","seq-rs649216","seq-t1d-19-60007809-C-G",
                       "seq-t1d-19-60034052-C-T","seq.t1d.19.60034052.C.T"))

snps.lrr <- make.names(c( "seq-rs4806568", "seq-rs62122181", "seq-rs674268",
              "seq-t1d-19-60056721-C-T"))

snps.theta <- make.names(c( "seq-rs592645", "seq-rs10500318", "seq-rs1654644",
                "seq-rs2295805", "seq-rs3865510", "seq-rs4806568",
                "seq-rs604999", "seq-rs62122181"))

library(ggplot2)
library(reshape)
d$LRR <- apply(d[,grep("Log_R",colnames(d))],1,mean)
df <- melt(d[,c("sample","cc","k1","k2","KIR","LRR",grep("^R|Theta|Log_R",colnames(d),value=TRUE))],
           id.vars=c("sample","cc","k1","k2","KIR"))
df$variable <- factor(make.names(as.character(df$variable)))
df <- df[order(df$sample,df$variable),]
ggplot(subset(df,variable %in% c("Log_R_Ratio.seq.rs3865510",
                                 "Theta.seq.rs1654644",
                                 "Theta.seq.rs649216",
                                 "Theta.seq.rs62122181",
                                 "Log_R_Ratio.rs581623",
                                 "Theta.seq.rs10500318",
                                 "Theta.seq.rs3865510",
                                 "Theta.seq.rs592645")),
              aes(x=variable,y=value,col=as.factor(KIR),group=sample)) +
  geom_point() + geom_path() + facet_grid(cc ~ .)

ggplot(subset(df,variable %in% c("LRR",
                                 "Log_R_Ratio.seq.rs3865510",
                                 "Theta.seq.rs1654644",
                                 "Theta.seq.rs649216",
                                 "Theta.seq.rs62122181",
                                 "Log_R_Ratio.rs581623",
                                 "Theta.seq.rs10500318",
                                 "Theta.seq.rs3865510",
                                 "Theta.seq.rs592645")),
              aes(x=variable,y=value,col=as.factor(KIR),group=sample)) +
  geom_point() + geom_path() + scale_colour_brewer(type="qual",palette="Dark2")

x <- d
colnames(x) <- make.names(colnames(x))

snps <- c("seq.rs3865510","seq.rs1654644","seq.rs649216","seq.rs62122181","seq.rs10500318","seq.rs592645")
theta <- df[grep("Theta",df$variable),]
lrr <- df[grep("Log_R",df$variable),]
r <- df[grep("^R",df$variable),]
theta$snp <- sub("Theta.","",theta$variable)
lrr$snp <- sub("Log_R_Ratio.","",lrr$variable)
r$snp <- sub("R.","",r$variable)
colnames(theta)[7] <- "theta"
colnames(lrr)[7] <- "lrr"
colnames(r)[7] <- "R"
wdf <- merge(theta[,-6],lrr[,-6],all=TRUE)
wdf <- merge(wdf, r[,-6], all=TRUE)
ggplot(wdf, aes(x=theta,y=lrr,col=KIR)) + geom_point() + facet_grid(snp~cc)

wsub <- subset(wdf,snp %in% snps & !is.na(KIR))
ggplot(wsub, aes(x=theta,y=R,col=KIR)) + geom_point() + facet_grid(snp~cc) + theme_bw()

wbest <- subset(wdf, snp=="seq.rs592645" & R>0.4)
wbest$qPCR <- ifelse(!is.na(wbest$KIR),"qPCR","without qPCR")
wbest$Copies <- wbest$KIR
ggplot(wbest, aes(x=theta,y=R,col=Copies)) +
  geom_point(data=subset(wbest, is.na(KIR))) +
  geom_point(data=subset(wbest, !is.na(KIR))) +
  facet_grid(qPCR~cc, margins=TRUE) + theme_bw()

GENO.COL <- data.frame( geno=c("?-?", "0-2", "1-1", "2-0", "2-1", "0-1", "1-2", "3-0", "1-0"), col=c('black', 'darkgreen', 'pink', 'blue', 'lightblue', 'green', 'purple', 'grey', 'red') )
g<-as.vector(GENO.COL[,2])
names(g)<-GENO.COL[,1]

## FIGURE
library(RColorBrewer)

#pdf("~/Words/papers/KIR/figures/fig-rs592645.pdf",height=6,width=8)
pdf('~nikolas/stats-archive/Papers/KIR/bmc_article/figures/Figure-2.pdf',height=6,width=9)
wbest <- wbest[ order(wbest$Copies,decreasing=TRUE), ]
wbest$Copies <- as.character(wbest$Copies)
wbest$Copies[ is.na(wbest$Copies) ] <- "?-?"
wbest$Copies <- factor(wbest$Copies,levels=c("?-?","0-1","0-2","1-2","1-1","2-1","1-0","2-0","3-0"))
ggplot(wbest, aes(x=theta,y=R,col=Copies)) +
  geom_point(data=subset(wbest, is.na(KIR)),size=1.3) +
  geom_point(data=subset(wbest, !is.na(KIR)),size=1.3) +
  facet_grid(qPCR~cc, margins=TRUE) + theme_bw() +
  #scale_colour_manual(breaks=c("?-?","0-1","0-2","1-2","1-1","2-1","1-0","2-0","3-0"), values=c("grey",rev(brewer.pal(n=8,name="Paired")[order(levels(wbest$Copies)[-1])]))) +
  scale_colour_manual(values=g, guide=guide_legend(override.aes=list(size=3),label.hjust=.5, label.vjust=.5, label.theme=element_text(angle=-45), label.position='top', title='KIR3DS1-KIR3DL1 Copy Number', title.position='right', title.vjust=.2, title.theme=element_text(face='italic',angle=0))) + theme(legend.position='top') +
  xlab(expression(theta)) + ylab('R')
dev.off()

vars <- c(paste0("Log_R_Ratio.","seq.rs592645"),
          paste0("R.","seq.rs592645"),
          paste0("Theta.","seq.rs592645"))
dfsub <- subset(df, variable %in% vars & !is.na(KIR))
ggplot(dfsub,
       aes(x=variable,y=value,col=KIR,group=sample)) + geom_point() + geom_path(alpha=0.2)
wsub <- subset(wdf,snp %in% "seq.rs592645" & !is.na(KIR))
ggplot(wsub,
       aes(x=theta,ymin=lrr,ymax=R,col=KIR,group=sample)) + geom_linerange() + geom_path(alpha=0.2)
ggplot(wsub,
       aes(x=lrr,y=R,col=KIR,group=sample)) + geom_point()
ggplot(wsub,
       aes(x=theta,y=lrr,col=KIR,group=sample)) + geom_point() + facet_grid(.~cc) + geom_hline(yintercept=c(-0.6,0,0.3))

################################################################################

## SUPP TABLE: informative snps
st<-data.frame(snp=unique(c(snps.theta, snps.r)))
st<-cbind(st, theta=st$snp %in% snps.theta, R=st$snp %in% snps.r)
ic <- read.table("/dunwich/scratch/chrisw/immunochip/final_chip/ImmunoChip_final_markerlist.txt.gz",
                 header=TRUE, as.is=TRUE)
ic$snp <- make.names(ic$Name)
st <- merge(st,ic,by.x="snp",by.y="snp")
st <- st[order(st$MapInfo),]
library(Hmisc)
#latex(st[,c(4:6,2:3)], file="~/Words/papers/KIR/figures/informative-snps.tex", booktabs=TRUE,rownames=NULL)

################################################################################

datasets.bak <- datasets

## k nearest neighbour 
library(class)
datasets <- lapply(datasets.bak, function(d) {
  colnames(d) <- make.names(colnames(d));
  wh <- grep("Log_R_Ratio",colnames(d))
  if(length(wh))
    d <- d[, -wh]
  return(d)
})
geno <- lapply(geno, function(x) {
  x$KIR <- paste0("cn",sub("-","",x$geno))
   rownames(x) <- as.character(geno[[1]]$immunochipsample)
  x
})
sapply(geno, function(x) table(x$KIR))
##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
## cn01   24   24   24   24   24   24   24   24   24    24
## cn02  890  890  890  890  890  890  890  890  890   890
## cn10    7    7    7    7    7    7    7    7    7     7
## cn11  436  436  436  436  436  436  436  436  436   436
## cn12   27   27   27   27   27   27   27   27   27    27
## cn20   52   55   53   54   53   53   52   54   53    54
## cn21   31   31   31   31   31   31   31   31   31    31
## cn30    7    4    6    5    6    6    7    5    6     5

## knn with best SNP
library(parallel)
options(mc.cores=18)
loo.all <- lapply(datasets, function(d) {
  d <- d[,-grep("GT_Score|GC_Score",colnames(d))]
  loo <- mclapply(geno, function(g) {
    d <- cbind(d, KIR=g[as.character(d$sample),"KIR"])
    train <- complete.cases(d)  
    mycv.knn(d[train,],k=1:20,trim.pattern="rs592645", do.plot=FALSE)
  })
  loo <- do.call("cbind",lapply(loo, function(x) x[,2]))
  mn <- rowMeans(loo)
  sd <- apply(loo,1,sd)
  xmin <- apply(loo,1,min)
  xmax <- apply(loo,1,max)
  data.frame(k=1:20,mn=mn,sd=sd,min=xmin,max=xmax)
})
for(n in names(datasets))
  loo.all[[n]]$group <- n
loo.best <- loo <- do.call("rbind",loo.all)
library(ggplot2)

pdf("~/Words/papers/KIR/figures/fig-knn-loo-rs592645.pdf",height=6,width=8)
ggplot(loo, aes(col=group,x=k+0.1*as.numeric(as.factor(group))-0.2,y=mn,ymin=min,ymax=max)) + 
  geom_pointrange() + xlab("k") + ylab("LOO error rate")
## mean + min/max error rate in LOO CV, averages/min/max taken over 10 MI datasets for the qPCR copy number calls
dev.off()

## knn with more snps
trim1 <- "rs592645"
trim3 <- "seq.rs592645|seq.t1d.19.60034052.C.T|R.seq.t1d.19.60056721.C.T"
snps.r
snps.theta
patt <- paste(unique(c(snps.r,snps.theta)),collapse="|")

trimpred <- paste(make.names(setdiff(c( "rs10422740", "rs10423751", "rs2569678",
"seq-rs11669355", "seq-rs12982080", "seq-rs270775", "seq-rs34078188",
"seq-rs604999", "seq-t1d-19-60054973-T-C", "seq-rs2295805",
"seq-rs62122181", "seq-rs62123050", "1kg_19_18101748", "seq-rs592645"),snps.drop)),collapse="|")


library(parallel)
options(mc.cores=20)
loo.all <- lapply(datasets, function(d) {
#  d <- d[,-grep("GT_Score|GC_Score",colnames(d))]
  loo <- mclapply(geno, function(g) {
    d <- cbind(d, KIR=g[as.character(d$sample),"KIR"])
    train <- complete.cases(d)  
    mycv.knn(d[train,],k=1:20,trim.pattern=trim3, do.plot=FALSE)
  })
  loo <- do.call("cbind",lapply(loo, function(x) x[,2]))
  mn <- rowMeans(loo)
  sd <- apply(loo,1,sd)
  xmin <- apply(loo,1,min)
  xmax <- apply(loo,1,max)
  data.frame(k=1:20,mn=mn,sd=sd,min=xmin,max=xmax)
})
for(n in names(datasets))
  loo.all[[n]]$group <- n
loo.lots <- loo <- do.call("rbind",loo.all)
loo[ which.min(loo$mn), ]

library(ggplot2)
pdf("~/Words/papers/KIR/figures/fig-knn-loo-thirteen.pdf",height=6,width=8)
ggplot(loo, aes(col=group,x=k+0.1*as.numeric(as.factor(group))-0.2,y=mn,ymin=min,ymax=max)) + 
  geom_pointrange()
dev.off()


## check min error in LOO
pdf("loo-cv-error-rate.pdf",height=6,width=6)
p2
dev.off()


## minimum at 3, sometimes 4.

save(loo, file="knn-mi-loo-030913.RData")


################################################################################

## vary size of training set, random sampling

d2 <- lapply(geno, function(g) {
  d <- cbind(datasets$both,
             KIR=g[as.character(datasets$both$sample),"KIR"])
#  d <- d[,-grep("GT_Score|GC_Score|Log_R",colnames(d))]
  train <- complete.cases(d)  
  d[train,]
})

library(parallel)
options(mc.cores=25)
nrep <- 25
results <- vector("list",9*nrep)
i <- 1
for(rep in 1:nrep) {
  for(p.train in c(seq(0.1,0.9,0.1),0.99)) {
    tmp <- do.call("cbind",mclapply(d2,function(d) {
      d <- d[,c(1,grep("Theta.seq.rs592645|R.seq.rs592645|KIR",colnames(d)))]
      varycv.knn(d,train=p.train,k=1:10)
    }))
    results[[i]] <- cbind(tmp[,c(1,2)], rowMeans(tmp[, colnames(tmp)=="e"]))
    i <- i+1
  }
}
results <- do.call("rbind",results)
colnames(results)[3] <- "e"
table(results$sample.size <- sprintf("%s%% (%s)",100*results$train,round(1474*results$train)))
results$fraction <- paste("Fraction",results$train)

pdf("~/Words/papers/KIR/figures/fig-knn-vary-sample-size.pdf",height=6,width=6)
ggplot(subset(results,fraction!="Fraction 0.99"),aes(x=k,y=e,col=sample.size)) +
  geom_hline(yintercept=0.02, linetype="dashed") +
  geom_point() + geom_smooth(se=FALSE) + facet_wrap(~sample.size) +
  ylab("Error rate") + theme(legend.position="none") + scale_x_continuous(breaks=seq(2,10,2)) +
  scale_y_continuous(breaks=seq(0,0.1,0.02)) + ylim(0,0.1)
dev.off()



colnames(datasets$both) <- make.names(colnames(datasets$both))
loo <- lapply(geno, function(g) {
  d <- cbind(datasets$both,
             KIR=g[as.character(datasets$both$sample),"KIR"])
  train <- complete.cases(d)  
  mycv.knn(d[train,],k=1:20,trim=FALSE)
})


################################################################################

## now call out of sample
data <- datasets$both

## drop samples excluded from Immunochip
dim(data) # 12269    90
data <- subset(data, as.character(sample) %in% as.character(lookup.clean$V2))
dim(data) # 12144    90

## nearest neighbour imputation
## use 3 nn
imputed.datasets <- lapply(geno, function(g) {
  d <- cbind(data, KIR=g[as.character(data$sample),"KIR"])
  d <- d[,-grep("GT|GC",colnames(d))]
  d <- knn.impute(d,trim.pattern="seq.rs592645")
  d$sample <- sub(".*\\.","",d$sample)
  return(d)
})

save(imputed.datasets, file="knn-mi-imputed-data-211113.RData")



################################################################################


## redo signal cloud figure
library(RColorBrewer)
wbest <- wbest[ sample(1:nrow(wbest)), ]
wbest$Copies <- as.character(wbest$Copies)
wbest$Copies[ is.na(wbest$Copies) ] <- "?-?"
wbest$Copies <- factor(wbest$Copies,levels=c("?-?","0-1","0-2","1-2","1-1","2-1","1-0","2-0","3-0"))

imp <- do.call("cbind",imputed.datasets)
imp <- data.frame(sample=imp[,1],
                  knn.best=apply(imp[,colnames(imp)=="knn.call"],1,function(x) {
                    tt <- table(x)
                    names(tt)[ which.max(tt) ]}),
                  knn.prob=rowMeans(imp[,colnames(imp)=="knn.prob"]))
imp <- merge(imp,wbest)
imp$Copies <- sub("cn","",imp$knn.best)
imp$Copies <- paste(substr(imp$Copies,1,1),substr(imp$Copies,2,2),sep="-")
head(imp)
imp$type <- "qPCR + Imputation"
wbest$type <- ifelse(wbest$Copies=="?-?","Extended sample","qPCR sample")
wbest2 <- rbind(wbest,imp[,colnames(wbest)])
wbest2$type <- factor(wbest2$type,levels=c("qPCR sample","Extended sample","qPCR + Imputation"))


pdf("fig-rs592645-v2.pdf",height=6,width=8)
ggplot(wbest2, aes(x=theta,y=R,col=Copies)) +
  geom_point(data=subset(wbest2, Copies=="?-?"),size=1) +
  geom_point(data=subset(wbest2, Copies!="?-?"),size=1) +
  facet_grid(type~cc, margins="cc") + theme_bw() +
  scale_colour_manual("3DS1-3DL1",
                      breaks=c("?-?","0-1","0-2","1-2","1-1","2-1","1-0","2-0","3-0"),
                      values= c("grey",rev(brewer.pal(n=8,name="Paired"))))
dev.off()

## disagreements
loo.knn <- function(wdata, KIR="KIR", k=8,trim=TRUE,trim.pattern="rs592645") {
  kcols <- grep("KIR",colnames(wdata))
  x <- wdata[,-c(1,kcols)]
  if(trim)
    x <- x[,grep(trim.pattern,colnames(x))]
  message("predicting using ",ncol(x)," columns:\n",paste(colnames(x),collapse=" "))
  y <- wdata[,KIR]
  m <- knn.cv(train=x, cl=y, k=k, prob=TRUE)
  tt <- table(m,y)
  x$y.pred <- m
  x$y.obs <- y
  x$sample <- wdata$sample
  return(x)
}

library(parallel)
options(mc.cores=20)
d <- datasets$both
d <- d[,-grep("GT_Score|GC_Score",colnames(d))]
loo.all <- lapply(geno, function(g) {
  d <- cbind(d, KIR=g[as.character(d$sample),"KIR"])
  train <- complete.cases(d)  
  loo.knn(d[train,],k=8,trim.pattern="rs592645")
})


loo <- as.data.frame(do.call("cbind",lapply(loo.all, function(x) as.character(x[,3]))))
loo$y.obs <- loo.all[[1]][,4]
loo$sample <- loo.all[[1]][,5]

mismatch <- loo[ apply(loo,1,function(x) any(x[1:10]!=x[[11]])),][,10:12]
colnames(mismatch)[1] <- "y.knn"
colnames(mismatch)[2] <- "y.qpcr"
write.table(mismatch,file="mismatch.csv",col.names=TRUE,row.names=FALSE,sep="\t")


mn <- rowMeans(loo)
  sd <- apply(loo,1,sd)
  xmin <- apply(loo,1,min)
  xmax <- apply(loo,1,max)
  data.frame(k=1:20,mn=mn,sd=sd,min=xmin,max=xmax)
})
for(n in names(datasets))
  loo.all[[n]]$group <- n
loo.best <- loo <- do.call("rbind",loo.all)
