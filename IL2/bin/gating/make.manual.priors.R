#!/usr/bin/env Rscript
source('~nikolas/bin/FCS/fcs.R',chdir=T)
suppressPackageStartupMessages(library(optparse))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowClust, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(tools, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))


flowClust2Prior <- function (x, kappa, Nt) {
    p <- ncol(x$mu)
    K <- x$K
    nu0 <- Ng <- x$w * Nt
    if (all((nu0 * kappa - p - 1) > 0)) {
        Lambda0 <- x$sigma
        for (i in 1:K) {
            Lambda0[i, , ] <- Lambda0[i, , ] * (kappa * nu0[i] - p - 1)
        }
    }
    else {
        stop("Can't proceed. Prior nu0 is negative for cluster(s) ", paste(which((nu0 - p - 1) > 0), collapse = ","), "\n(p-1) = ", p - 1, ": Try increasing kappa")
    }
    Omega0 <- array(0, c(K, p, p))
    for (i in 1:K) {
        Omega0[i, , ] <- diag(1, p)
        if (p == 1) {
            dS <- x$sigma[i, , ]
            dO <- Omega0[i, , ]
        }
        else {
            dS <- det(x$sigma[i, , ])
            dO <- det(Omega0[i, , ])
        }
        k <- (dO/dS)^(1/p)
        Omega0[i, , ] <- Omega0[i, , ] * k
        Omega0[i, , ] <- solve(Omega0[i, , ] * Ng[i] * kappa)
    }
    nu0 <- nu0 * kappa
    Mu0 <- x$mu
    lambda <- x$lambda
    w0 <- x$w * Nt
    prior <- list(Mu0 = Mu0, Lambda0 = Lambda0, Omega0 = Omega0, w0 = w0, nu0 = nu0, nu = x$nu, lambda = x$lambda, K = K)
    class(prior) <- "flowClustPrior"
    attr(prior, "lambda") <- x$lambda
    #name all array dimensions for later subsetting
    channels <- colnames(prior$Mu0)
    prior$Lambda0 <- array(prior$Lambda0,dim=c(2,2,2),dimnames=list(1:2,channels,channels))
    prior$Omega0 <- array(prior$Omega0,dim=c(2,2,2),dimnames=list(1:2,channels,channels))
    prior
}


# kappa > 1
# smaller kappa means less certainty
# one value of kappa probably doesn't fit all
# Returns two-component with 
# background cluster (position 1)
# and prior for gate (position 2).
manual.prior <- function(K, channels, kappa=100, Nt=NULL) {
    g <- as.matrix(G[,K])
    K <- as.numeric(as.factor(K))
    d <- fcs.data.prior[,channels] 
    if (is.null(Nt)) Nt <- nrow(d)
    p <- length(channels)
    f <- lapply(K, function(i) which(as.logical(g[,i])))
    #prepend the background cluster mean
    Mu <- rbind(colMeans(d), do.call('rbind', lapply(K, function(i) colMeans(d[f[[i]],]))) )
    Sigma <- array(dim=c(length(K)+1,p,p))
    #prepend the background cluster covariance
    Sigma[1,,] <- cov(d)
    for (i in K) Sigma[i+1,,] <- cov(d[f[[i]],])
    w <- sapply(K, function(i) length(f[[i]])/length(d[,1]))
    w <- c(1-sum(w), w) 
    theta <- list(K=length(K)+1, mu=Mu, sigma=Sigma, w=w, lambda=1, nu=4)
    prior <- flowClust2Prior(theta, kappa=kappa, Nt=Nt)
    return(prior)
}

#
manual.prior2 <- function(K, channels, kappa=100, Nt=NULL) {
    g <- as.matrix(G[,K])
    K <- as.numeric(as.factor(K))
    d <- fcs.data.prior[,channels] 
    if (is.null(Nt)) Nt <- nrow(d)
    p <- length(channels)
    f <- lapply(K, function(i) which(as.logical(g[,i])))
    #prepend the background cluster mean
    Mu <- rbind(colMeans(d), do.call('rbind', lapply(K, function(i) colMeans(d[f[[i]],]))) )
    Sigma <- array(dim=c(length(K)+1,p,p))
    #prepend the background cluster covariance
    Sigma[1,,] <- cov(d)
    for (i in K) Sigma[i+1,,] <- cov(d[f[[i]],])
    w <- sapply(K, function(i) length(f[[i]])/length(d[,1]))
    w <- c(1-sum(w), w) 
    theta <- list(K=length(K)+1, mu=Mu, sigma=Sigma, w=w, lambda=1, nu=4)
    prior <- flowClust2Prior(theta, kappa=kappa, Nt=Nt)
    return(prior)
}




option_list <- list( 
make_option(c("--in.file"), default=NULL, help = "RData file (must contain fcs.data, must have same number of rows as CLR file)"),
make_option(c("--clr"), default=NULL, help = "CLR file to use (must have same number of rows as RData file)"),
make_option(c("--out.dir"), default=NULL, help = "dir to send gates (.prior.RData files)")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

#Lymphocytes
#Lymphocytes.Single.Cells
#Lymphocytes.Single.Cells.CD4
#Lymphocytes.Single.Cells.CD4.Memory
#Lymphocytes.Single.Cells.CD4.Naive

#print(channels <- grep('PSTAT5', colnames(fcs.data), invert=TRUE, value=TRUE))
#print(channels <- c('FSCA','SSCA','CD4','CD25','CD45RA','CD3','CD56','FOXP3','CD8'))
#print(channels <- c('FSCA','SSCA','CD4'))


print(channels <- c("FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3"))

#setwd('~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3') 
#print(colnames(G <- read.csv('CLR/I022267C_CB00010K_0U.clr'))) 
#print(load('pstat5-join/All/CB00010K_2012-11-13.RData'))
#print(load('pstat5-join/All/CB01494Y_2012-10-09.RData'))

print(gates <- colnames(G <- read.csv(opt$clr))) 
print(load(opt$in.file))
fcs.data.prior <- fcs.data

#the CD4 lymphocyte prior
for (g in gates) {
    print(g)
    prior <- manual.prior(K=g, channels=channels, kappa=1000, Nt=1) 
    print(out.file <- file.path(opt$out.dir,sprintf("%s.prior.RData",g)))
    save(prior, file=out.file)
}



