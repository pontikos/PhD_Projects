#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list( 
make_option(c("-f","--fcs"), help = "fcsfile to parse"),
make_option(c("-p","--parameter"), help="parameter on which to gate")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowCore, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(cluster, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(mixtools, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))

suppressMessages(suppressWarnings( read.FCS(opt$fcs, alter.names=TRUE, column.pattern=".A") -> fcs.data ))
dimnames(fcs.data@exprs)[[2]] -> gating.parameters
gating.parameters[!grepl("^Time|^FSC|^SSC",gating.parameters)]->gating.parameters
cat("All Parameters:", gating.parameters, "\n")
gating.parameters[grep(opt$parameter,gating.parameters)]->gating.parameter
#only gate on fluorescent channels
cat("Gating Parameter:", gating.parameter, "\n")


K <- 2

#x is the original data
fcs.data@exprs[,gating.parameter]->x
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

#range between medoids
y[medoids[1] < y & y < medoids[2]]->y.valley
#100 bins good number?
hist(y.valley, breaks=100, plot=F)->h
#smallest bin
h$mid[which.min(h$density)]->thresholds$valley


### parametric mixture model
normalmixEM(y, mu=medoids) -> mm
max(y[mm$posterior[,1]>.95]) -> thresholds$mm.1
mm$mu[[1]]  -> thresholds$mu.1
mm$sigma[[1]]  -> thresholds$sigma.1
min(y[mm$posterior[,2]>.95]) -> thresholds$mm.2
mm$mu[[2]] -> thresholds$mu.2
mm$sigma[[2]] -> thresholds$sigma.2

### non parametric mixture model
spEMsymloc(y, mu=medoids, bw=.1) -> np.mm
np.mm$mu[[1]]  -> thresholds$mu.1
np.mm$sigma[[1]]  -> thresholds$sigma.1

max(y[mm$posterior[,1]>.95]) -> thresholds$np.mm.1
min(y[mm$posterior[,2]>.95]) -> thresholds$np.mm.2



### blank bead threshold

### print out
#header
cat('>>fcsFile')
for (threshold.name in names(thresholds)) cat(',', threshold.name, sep='')
cat("\n")
#values
cat('>>', opt$fcs, sep='')
for (threshold.name in names(thresholds)) {
    cat(',', thresholds[[threshold.name]], sep='')
}
cat("\n")

