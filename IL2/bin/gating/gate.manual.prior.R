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
    prior
}


flowClust <- function(d,B,prior)
    flowClust::flowClust(d[,colnames(prior$Mu0)],K=prior$K,B=B,usePrior='yes',prior=prior,level=.9,u.cutoff=.5,trans=0,lambda=1,control=list(B.lambda=0))


option_list <- list( 
make_option(c("--in.file"), default=NULL, help = ".RData file"),
make_option(c("--prior"), default=NULL, help = ".prior.RData file"),
make_option(c("--channels"), default=NULL, help = ""),
make_option(c("--plot.file"), default='', help = ""),
make_option(c("--out.file"), default='', help = "")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

#Lymphocytes
#Lymphocytes.Single.Cells
#Lymphocytes.Single.Cells.CD4
#Lymphocytes.Single.Cells.CD4.Memory
#Lymphocytes.Single.Cells.CD4.Naive

channels <- unlist(strsplit(opt$channels, ','))

load(opt$prior)
prior$Mu0 <- prior$Mu0[,channels]
prior$Lambda0 <- prior$Lambda0[,channels,channels]
prior$Omega0 <- prior$Omega0[,channels,channels]
print(prior)

load(opt$in.file)
original.fcs.data <- fcs.data
fcs.data <- fcs.data[,channels]

res.1 <- flowClust(fcs.data, 1, prior)
res <- flowClust(fcs.data, 2000, prior)

#postscript('file.eps', colormodel="cmyk")
plotClustRes(fcs.data,res,before.res=res.1, outliers=FALSE,plot.file=opt$plot.file,cex=1.5)

clusters <- Map(res)
fcs.data <- original.fcs.data[which(clusters==2),]

save(fcs.data, res, file=opt$out.file)


#print(channels <- grep('PSTAT5', colnames(fcs.data), invert=TRUE, value=TRUE))
#print(channels <- c('FSCA','SSCA','CD4','CD25','CD45RA','CD3','CD56','FOXP3','CD8'))
#print(channels <- c('FSCA','SSCA','CD4'))
#print(channels <- c("FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3"))


