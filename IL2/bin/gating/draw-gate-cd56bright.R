#!/usr/bin/env Rscript
source('~nikolas/bin/FCS/fcs.R')
source('~nikolas/Projects/IL2/bin/common.R')
source('~nikolas/Projects/IL2/bin/gating/flowclust-functions.R')
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(sp))
suppressPackageStartupMessages(library(flowClust))

### interactive mode does not work :(
#option_list <- list( 
#make_option(c("--individual"), default=NULL, help = "individual"),
#make_option(c("--date"), default=NULL, help = "date")
#)
#OptionParser(option_list=option_list) -> option.parser
#parse_args(option.parser) -> opt

#
save.gate <- function(f) {
    print( f <- file.path(BASE.DIR, 'CLR', f) )
    save(CLR, file=f)
}

BASE.DIR <- '~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/'
#for (f in list.files(file.path(BASE.DIR,'FCS','pstat5-join'), pattern='.*.RData')) { 
load('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/transform-w1.RData')

panel <-list("FSC-A"="FSCA", "FSC-H"="FSCH", "FSC-W"="FSCW", "SSC-A"="SSCA", "SSC-H"="SSCH", "SSC-W"="SSCW", "Alexa Fluor 488-A"="PSTAT5", "PerCP-Cy5-5-A"="CD3", "APC-A"="CD25", "Alexa Fluor 700-A"="CD4", "APC-Cy7-A"="APC-Cy7-A", "Pacific Blue-A"="CD56", "Qdot 605-A"="CD8", "PE YG-A"="FOXP3", "PE-Cy7 YG-A"="CD45RA", "Time"="Time") 
#f <- '/dunwich/scratch/nikolas/FCS.Tony/pSTAT5_DGAP_KM00782Z/pSTAT5_DGAP_KM00782Z_0U.fcs'
f <- '/dunwich/scratch/nikolas/FCS.Tony/DGAP_NKCD8_KA866533H/DGAP_NKCD8_KA866533H_0U.fcs'
f1 <- flowCore::read.FCS(f)
d <- f1@exprs
colnames(d) <- panel[colnames(f1)]
fcs.data <- d[,c("FSCA", "FSCH", "FSCW", "SSCA", "SSCH", "SSCW", "PSTAT5", "CD3",  "CD25", "CD4", "CD56", "CD8", "FOXP3", "CD45RA")]
print(colnames(fcs.data))
fcs.data <- applyTransforms(fcs.data, transforms) 
print(names(transforms))

fcs.data<-as.data.frame(fcs.data)
X <- fcs.data[ which( .5 < fcs.data$CD25 & fcs.data$CD25 < 1.5 & 2 < fcs.data$CD3 & fcs.data$CD3 < 2.5 & .5 < fcs.data$CD4 & fcs.data$CD4 < 1.5 & 3.5 < fcs.data$CD45RA &  2 < fcs.data$CD56 & .5 < fcs.data$CD8 & fcs.data$CD8 < 1 ),]

par(mfrow=c(3,3))
for (marker in CORE.MARKERS) smoothPlot1D(X[,marker], outliers=TRUE, main=marker)

gate.f <- file.path(BASE.DIR,'CLR','KM00782Z_0U_2012-07-24.RData')
load(gate.f)

 #Lymphocytes
background.mean <- colMeans(fcs.data) 
prior <- list() 

#Lymphocytes
channels <- c('CD3','CD4')
X <- fcs.data[as.logical(CLR[,'Lymphocytes']),channels]
smoothPlot(X,outliers=FALSE, main='CD3- CD4- ?')
G <- locator(type='l') 
#G <- structure(list(CD3 = c(2.68970800837676, 2.44841544366885, 2.07598561553272,  1.78223814545352, 1.8871479561961, 2.15466797358965, 2.62676212193122,  2.56381623548567), CD4 = c(1.58187141454148, 1.96727880340973,  2.08486071865768, 1.69292100116454, 1.22912566879766, 0.954767866552464,  1.22259334017277, 1.58187141454148)), .Names = c("CD3", "CD4" ))
names(G) <- channels
in.poly <- point.in.polygon( X[,1] , X[,2], G[[1]], G[[2]])
e <- classification.to.ellipse(X,in.poly)
plot.ellipses(list(e))
CLR <- cbind(CLR,'CD3.CD4'=0)
CLR[as.logical(CLR[,"Lymphocytes"]),'CD3.CD4'] <- in.poly

#  
channels <- c('CD25','CD56')
X <- fcs.data[as.logical(CLR[,'Lymphocytes'])&as.logical(CLR[,'CD3.CD4']),channels]
smoothPlot(X,outliers=FALSE, main='CD25- CD56hi?')
#G <- locator(type='l') 
G <- structure(list(CD25 = c(1.30234719696643, 1.04830101072666, 1.07727118985927,  1.29120482037697, 1.5898205129746, 1.70347275418713, 1.51405235216624,  1.30903262292011), CD56 = c(3.05731610501181, 2.96773903242385,  2.78005373747763, 2.74166356351136, 2.77578816259249, 2.95920788265356,  3.03598823058611, 3.05731610501181)), .Names = c("CD25", "CD56" ))
names(G) <- channels
in.poly <- point.in.polygon( X[,1] , X[,2], G[[1]], G[[2]])
e <- classification.to.ellipse(X,in.poly)
plot.ellipses(list(e))
CLR <- cbind(CLR,'CD25.CD56'=0)
CLR[as.logical(CLR[,'Lymphocytes'])&as.logical(CLR[,'CD3.CD4']),'CD25.CD56'] <- in.poly

CLR <- cbind(CLR,'CD56bright'=0)
CLR[as.logical(rowSums(CLR[,c('CD3.CD4','CD25.CD56')])==2),'CD56bright']<-1
apply(CLR,2,table)
save(CLR,file=gate.f)


channels <- c('CD3','CD56')
X <- fcs.data[,channels]
smoothPlot(X,outliers=FALSE, main='CD25- CD56+ ?')
#G <- locator(type='l') 
G <- structure(list(CD25 = c(1.30234719696643, 1.04830101072666, 1.07727118985927,  1.29120482037697, 1.5898205129746, 1.70347275418713, 1.51405235216624,  1.30903262292011), CD56 = c(3.05731610501181, 2.96773903242385,  2.78005373747763, 2.74166356351136, 2.77578816259249, 2.95920788265356,  3.03598823058611, 3.05731610501181)), .Names = c("CD25", "CD56" ))
names(G) <- channels
in.poly <- point.in.polygon( X[,1] , X[,2], G[[1]], G[[2]])
e <- classification.to.ellipse(X,in.poly)
plot.ellipses(list(e))

CLR <- cbind(CLR,'CD25.CD56'=0)
CLR[,'CD25.CD56'] <- in.poly




#
plot.ellipses(list(classification.to.ellipse(X,CLR[,'Lymphocytes'])))
save.gate(f) 
#prior[['Lymphocytes']] <- manual.prior(X,1+CLR[,'Lymphocytes']) 
#Single cells
X <- fcs.data[which(as.logical(CLR[,'Lymphocytes'])),c('SSCH','SSCW')]
smoothPlot(X,outliers=FALSE,main='Single cells?')
#
G <- locator(type='l') 
in.poly <- point.in.polygon( X[,1] , X[,2], G[[1]], G[[2]])
#in.poly <- point.in.polygon( X[,1] , X[,2], G$x, G$y)
e <- classification.to.ellipse(X,in.poly)
plot.ellipses(list(e)) 
CLR[,'Single cells'] <- 0
CLR[which(as.logical(CLR[,'Lymphocytes'])),'Single cells'] <- in.poly
#
plot.ellipses(list(classification.to.ellipse(X,CLR[which(as.logical(CLR[,'Lymphocytes'])),'Single cells'])))
save.gate(f) 
#prior[['Single cells']] <- manual.prior(X,1+CLR[which(as.logical(CLR[,'Lymphocytes'])),'Single cells']) 
#CD4+
X <- fcs.data[which(as.logical(CLR[,'Single cells'])),c('CD4','SSCA')]
smoothPlot(X,outliers=FALSE, main='CD4+ ?')
#
G <- locator(type='l') 
in.poly <- point.in.polygon( X[,1] , X[,2], G$x, G$y)
e <- classification.to.ellipse(X,in.poly)
plot.ellipses(list(e)) 
CLR[,'CD4'] <- 0
CLR[which(as.logical(CLR[,'Single cells'])),'CD4'] <- in.poly
#
plot.ellipses(list(classification.to.ellipse(X,CLR[which(as.logical(CLR[,'Single cells'])),'CD4'])))
save.gate(f) 
#prior[['CD4']] <- manual.prior(X,1+CLR[which(as.logical(CLR[,'Single cells'])),'CD4']) 
# Memory / Naive
X <- fcs.data[which(as.logical(CLR[,'CD4'])),c('CD45RA','SSCA')]
smoothPlot(X,outliers=TRUE, main='Memory?')
# Memory
G <- locator(type='l') 
in.poly <- point.in.polygon( X[,1] , X[,2], G$x, G$y)
e <- classification.to.ellipse(X,in.poly)
plot.ellipses(list(e))
CLR[,'Memory'] <- 0
CLR[which(as.logical(CLR[,'CD4'])),'Memory'] <- in.poly
#
smoothPlot(X,outliers=TRUE, main='Naive?')
plot.ellipses(list(classification.to.ellipse(X,CLR[which(as.logical(CLR[,'CD4'])),'Memory'])))
# Naive 
G <- locator(type='l') 
in.poly <- point.in.polygon( X[,1] , X[,2], G$x, G$y)
e <- classification.to.ellipse(X,in.poly)
plot.ellipses(list(e))
CLR[,'Naive'] <- 0
CLR[which(as.logical(CLR[,'CD4'])),'Naive'] <- in.poly
#
plot.ellipses(list(classification.to.ellipse(X,CLR[which(as.logical(CLR[,'CD4'])),'Naive'])))
save.gate(f) 
#prior[['Memory / Naive']] <- manual.prior(X,CLR[which(as.logical(CLR[,'CD4'])),'Memory']+2*CLR[which(as.logical(CLR[,'CD4'])),'Naive']) 
# Naive
X <- fcs.data[which(as.logical(CLR[,'Naive'])),c('CD25','FOXP3')]
smoothPlot(X,outliers=TRUE, main='Naive Eff?')
plot.gate.chull(fcs.data[,c('CD25','FOXP3')], classification=(CLR[,'Naive Eff']))
plot.gate.chull(fcs.data[,c('CD25','FOXP3')], classification=(CLR[,'Naive Treg']))
# Naive Eff
G <- locator(type='l') 
in.poly <- point.in.polygon( X[,1] , X[,2], G$x, G$y)
e <- classification.to.ellipse(X,in.poly)
plot.ellipses(list(e))
CLR[,'Naive Eff'] <- 0
CLR[which(as.logical(CLR[,'Naive'])),'Naive Eff'] <- in.poly
#
smoothPlot(X,outliers=TRUE, main='Naive Treg?')
plot.ellipses(list(classification.to.ellipse(X,CLR[which(as.logical(CLR[,'Naive'])),'Naive Eff'])))
# Naive Treg
G <- locator(type='l') 
in.poly <- point.in.polygon( X[,1] , X[,2], G$x, G$y)
e <- classification.to.ellipse(X,in.poly)
plot.ellipses(list(e))
CLR[,'Naive Treg'] <- 0
CLR[which(as.logical(CLR[,'Naive'])),'Naive Treg'] <- in.poly
plot.ellipses(list(classification.to.ellipse(X,CLR[which(as.logical(CLR[,'Naive'])),'Naive Treg']))) 
save.gate(f) 
#prior[['Naive Eff / Treg']] <- manual.prior(X,1+CLR[which(as.logical(CLR[,'Naive'])),'Naive Eff']+2*CLR[which(as.logical(CLR[,'Naive'])),'Naive Treg'], outliers=TRUE) 
# Memory
X <- fcs.data[which(as.logical(CLR[,'Memory'])),c('CD25','FOXP3')]
smoothPlot(X,outliers=TRUE, main='Memory Eff?')
plot.gate.chull(fcs.data[,c('CD25','FOXP3')], classification=(CLR[,'Memory Eff']))
plot.gate.chull(fcs.data[,c('CD25','FOXP3')], classification=(CLR[,'Memory Treg']))
# Memory Eff
G <- locator(type='l') 
in.poly <- point.in.polygon( X[,1] , X[,2], G$x, G$y)
e <- classification.to.ellipse(X,in.poly)
plot.ellipses(list(e))
CLR[,'Memory Eff'] <- 0
CLR[which(as.logical(CLR[,'Memory'])),'Memory Eff'] <- in.poly
#
smoothPlot(X,outliers=TRUE, main='Memory Treg?')
plot.ellipses(list(classification.to.ellipse(X,CLR[which(as.logical(CLR[,'Memory'])),'Memory Eff'])))
# Memory Treg
G <- locator(type='l') 
in.poly <- point.in.polygon( X[,1] , X[,2], G$x, G$y)
e <- classification.to.ellipse(X,in.poly)
plot.ellipses(list(e))
CLR[,'Memory Treg'] <- 0
CLR[which(as.logical(CLR[,'Memory'])),'Memory Treg'] <- in.poly
plot.ellipses(list(classification.to.ellipse(X,CLR[which(as.logical(CLR[,'Memory'])),'Memory Treg']))) 
save.gate(f) 
#prior[['Memory Eff / Treg']] <- manual.prior(X,1+CLR[which(as.logical(CLR[,'Memory'])),'Memory Eff']+2*CLR[which(as.logical(CLR[,'Memory'])),'Memory Treg'], outliers=TRUE) 
#print(save(prior,transforms,file=file.path(BASE.DIR,'flowclust-priors.RData'))) 

#}

