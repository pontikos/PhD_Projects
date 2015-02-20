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

# K is the cluster labels, must be numerically ordered starting at 1
# kappa > 1
# smaller kappa means less certainty
# one value of kappa probably doesn't fit all
# Returns two-component with 
# background cluster (position 1)
# and prior for gate (position 2).
# Probably better for the background to be everything.
manual.prior <- function(d, K, kappa=100, Nt=NULL, outliers=TRUE) {
    print(sort(unique(K)))
    p <- ncol(d)
    if (is.null(Nt)) Nt <- nrow(d)
    Mu <- do.call('rbind', by(d, K, colMeans))
    Mu[1,] <- background.mean[colnames(d)]
    Sigma <- array(dim=c(length(unique(K)),p,p))
    for (i in sort(unique(K))) Sigma[i,,] <- cov(d[which(K==i),])
    Sigma[1,,] <- background.cov[colnames(d),colnames(d)]
    w <- as.numeric(prop.table(table(K)))
    w[1] <- .0001
    theta <- list(K=as.numeric(length(unique(K))), mu=Mu, sigma=Sigma, w=w, lambda=1, nu=4)
    print(prior <- flowClust2Prior(theta, kappa=kappa, Nt=Nt))
    level=.9
    u.cutoff=.5
    B=500
    kappa=1
    print(colnames(d))
    res <- flowClust::flowClust(d,varNames=colnames(d),K=prior$K,B=B,level=level,u.cutoff=u.cutoff,trans=0,lambda=1,control=list(B.lambda=0),usePrior='yes',prior=prior)
    smoothPlot(X, classification=MAP(X, res), outliers=outliers)
    return(prior)
}

#
save.gate <- function(f) {
    print( f <- file.path(BASE.DIR, 'CLR', f) )
    save(CLR, file=f)
}


BASE.DIR <- '~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3'
f <- 'CB01494Y_2012-10-09.RData'

BASE.DIR <- '~/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/'
f <- 'CB00573X_2012-06-12.RData'

BASE.DIR <- '~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3'
f <- 'KM00861S_2012-09-20.RData'

BASE.DIR <- '~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3'
f <- 'KM00915T_2012-10-25.RData'

#for (f in list.files(file.path(BASE.DIR,'FCS','pstat5-join'), pattern='.*.RData')) { 
print(f)
print((load(file.path(BASE.DIR,'RData/pstat5-join/',f))))

load('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/transform-w1.RData')

panel <-list("FSC-A"="FSCA", "FSC-H"="FSCH", "FSC-W"="FSCW", "SSC-A"="SSCA", "SSC-H"="SSCH", "SSC-W"="SSCW", "Alexa Fluor 488-A"="PSTAT5", "PerCP-Cy5-5-A"="CD3", "APC-A"="CD25", "Alexa Fluor 700-A"="CD4", "APC-Cy7-A"="APC-Cy7-A", "Pacific Blue-A"="CD56", "Qdot 605-A"="CD8", "PE YG-A"="FOXP3", "PE-Cy7 YG-A"="CD45RA", "Time"="Time") 
f <- '/dunwich/scratch/nikolas/FCS.Tony/pSTAT5_DGAP_KM00782Z/pSTAT5_DGAP_KM00782Z_0U.fcs'
f1 <- flowCore::read.FCS(f)
d <- f1@exprs
colnames(d) <- panel[colnames(f1)]
fcs.data <- d[,c("FSCA", "FSCH", "FSCW", "SSCA", "SSCH", "SSCW", "PSTAT5", "CD3",  "CD25", "CD4", "CD56", "CD8", "FOXP3", "CD45RA")]

print(colnames(fcs.data))
CLR <- matrix(0, nrow=nrow(fcs.data),ncol=length(CLR.CELL.TYPES))
colnames(CLR) <- CLR.CELL.TYPES 
fcs.data <- applyTransforms(fcs.data, transforms) 
#Lymphocytes
print(names(transforms))
#transforms[['CD3']] <- transforms[['CD4']]
#transforms[['CD56']] <- transforms[['CD4']]
#transforms[['CD8']] <- transforms[['CD4']]
background.cov <- cov(fcs.data)
background.mean <- colMeans(fcs.data) 
prior <- list() 
#Lymphocytes
channels <- c('FSCA','SSCA')
X <- fcs.data[,c('FSCA','SSCA')]
smoothPlot(X,outliers=FALSE, main='Lymphocytes?')
#
G <- locator(type='l') 
names(G) <- channels
in.poly <- point.in.polygon( X[,1] , X[,2], G[[1]], G[[2]])
e <- classification.to.ellipse(X,in.poly)
plot.ellipses(list(e))
CLR[,'Lymphocytes'] <- 0
CLR[,'Lymphocytes'] <- in.poly
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

