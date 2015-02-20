#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list( 
make_option(c("-f","--fcs"), help = "fcsfile to parse"),
make_option(c("-p","--parameters"), help="parameters on which to gate")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowCore, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(cluster, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(mixtools, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))


unlist(strsplit(opt$parameters, ",")) -> gating.parameters

suppressMessages(suppressWarnings( read.FCS(opt$fcs, alter.names=TRUE, column.pattern=".A") -> fcs.data ))
dimnames(fcs.data@exprs)[[2]] -> parameters
parameters[!grepl("^Time|^FSC|^SSC",parameters)]->parameters
cat("All Parameters:", parameters, "\n")

parameters[sapply(gating.parameters, function(x) grep(x, parameters))]->gating.parameters
#only gate on fluorescent channels
cat("Gating Parameter:", gating.parameters, "\n")

K <- 2

log10(fcs.data@exprs[,gating.parameters]) -> x


x[rowSums(x>0 & x<2)==2,] -> x

pam(x, k=2)->pam.res

print(pam.res$medoids)

### bivariate mixture of normals
mvnormalmixEM(x, mu=list(pam.res$medoids[1,], pam.res$medoids[2,]), k=2) -> mm.2d
save(mm.2d, file=sprintf('mm.2d.%s.obj', opt$fcs))

exit(0)




#x is the original data
print(length(x))
print(mean(x))
print(quantile(x))

#y is the transformed data for clustering
log10(x[x>1])->y
print(length(y))
print(mean(y))
print(quantile(y))

thresholds <- list()

# most basic split just take the mean
mean(y) -> thresholds$mean

#partition around medoid to identify cluster medoid
#i.e points which minimise the overall inter cluster distance
pam(y, k=K)->res

sort(res$medoids)->medoids
print(medoids)

thresholds$medoids.1 <- medoids[[1]]
thresholds$medoids.2 <- medoids[[2]]

mean(medoids) -> thresholds$medoids.mean


### parametric mixture model
#normalmixEM(y, mu=medoids) -> mm
#save(mm, file=sprintf('mm.%s.obj',opt$fcs))

### semi parametric symmetric mixture model
#spEMsymloc(y, mu=medoids, bw=.01) -> np.mm
#save(np.mm, file=sprintf('np.mm.%s.2.obj', opt$fcs))


### non parameteric mixture model
#npEM(

#fcs.data@exprs[,

