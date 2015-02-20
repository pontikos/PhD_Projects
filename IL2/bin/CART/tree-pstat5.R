#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tree))
suppressPackageStartupMessages(library(rpart))
#suppressPackageStartupMessages(library(mvpart))
suppressPackageStartupMessages(library(partykit))
suppressPackageStartupMessages(library(diptest))
source('~nikolas/bin/FCS/fcs.R',chdir=T)

option_list <- list( 
make_option(c("--in.file"), default=NULL, help = ".RData file"),
make_option(c('--channels'), default='FSCW,SSCW,FSCH,SSCH,FSCA,SSCA,CD4,CD25,CD45RA,FOXP3', help=''),
make_option(c("--cell.subset"), default=NULL, help = "any of the subsets defined in the CLR file"),
make_option(c("--clr"), default='magnetic-manual-gates', help = "any of the subsets defined in the CLR file"),
make_option(c("--out.dir"), default='~/Plots-tree-pstat5/', help = "Output directory"),
make_option(c("--splits"), default=NULL, help = "number of splits")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

if (!is.null(opt$channels)) {
    channels <- unlist(strsplit(opt$channels, ","))
} else  {
    channels <- NULL
}

base.dir <- '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All'
#file.name <- 'CB01477E_2012-09-26.RData'
file.name <- opt$in.file
file.name <- file.path(base.dir,file.name)
print(basename(file.name))
clr <- 'magnetic-manual-gates2'
clr <- opt$clr
gate.dir <- sprintf('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/%s/CLR/',clr)

#splits <- 3
splits <- opt$splits
#cell.subset <- 'CD4'
cell.subset <- opt$cell.subset

print(load(file.name))
print(load(file.path(gate.dir,basename(file.name))))

fcs.data <- baseline.relative.pstat5(fcs.data) 
if (!is.null(cell.subset)) fcs.data <- fcs.data[which(as.logical(CLR[,cell.subset])),] 

channels <- c("FSCW", "SSCW", "FSCH", "SSCH", "FSCA", "SSCA", "CD4", "CD25", "CD45RA", "FOXP3")
form <- formula(sprintf('diff.PSTAT5.4 ~ %s', paste(channels, collapse=' + ')))

t <- tree::tree(form, data.frame(fcs.data),split='deviance')
if (!is.null(splits)) t <- tree::prune.tree(t,best=splits)
print(t$frame) 
w <- as.factor(t$where)
levels(w) <- 1:length(unique(w))
table(w)

par(mfrow=c(2,5))
for (i in sort(unique(w))) {
    x <- fcs.data[which(w==i),'diff.PSTAT5.4']
    plot(normalised.density(x), main=mean(x)/sd(x))
}

for (i in sort(unique(w)) {
    t <- tree::tree(form, data.frame(fcs.data[which(w==i),]),split='deviance')
    w <- as.factor(t$where)
    levels(w) <- 1:length(unique(w))
    table(w)
}

channels <- c("CD4", "CD25", "CD45RA", "FOXP3")

recursive.part <- function(d) {
    x <- d[,'diff.PSTAT5.4']
    #not significantly bimodal
    if (diptest::dip.test(x)$p.value > 2.2e-16) {
        cat('unimodal',nrow(d),'\n')
        return(d)
    }
    #t <- tree::tree(form, data.frame(d),split='deviance',tree.control=tree.control(ql)
    t <- rpart::rpart(form, data=data.frame(d),method=list(eval=rpart.eval, split=rpart.split, init=rpart.init),control=rpart.control(minsplit=60,minbucket=20,usesurrogate=0))
    w <- as.factor(t$where)
    levels(w) <- 1:length(unique(w))
    print(table(w))
    splits <- sort(unique(w))
    if (length(splits)==1) {
        cat('cannot split further', nrow(d),'\n')
        return(d)
    }
    l <- list()
    for (i in splits) {
        x <- d[which(w==i),'diff.PSTAT5.4']
        #significantly bimodal
        if (diptest::dip.test(x)$p.value < 2.2e-16) l[[i]] <- recursive.part(d[which(w==i),])
        else l[[i]] <- d
    }
    return(l)
} 
X <- recursive.part(fcs.data) 
#X <- recursive.part(scaled.fcs.data <- apply(fcs.data,2,scale))


t <- rpart::rpart('diff.PSTAT5.4 ~ FSCA + SSCA', data=data.frame(fcs.data),method=list(eval=rpart.eval, split=rpart.split, init=rpart.init),control=rpart.control(minsplit=60,minbucket=20,usesurrogate=0))

###
t <- tree::tree('diff.PSTAT5.4 ~ FSCA + SSCA', data=data.frame(fcs.data))
t <- prune.tree(t, best=3)
w <- as.factor(t$where)
levels(w) <- 1:length(unique(w))
table(w) 
plotClusters(fcs.data[,c('diff.PSTAT5.4','FSCA','SSCA')],classification=w,ellipses=FALSE,chull.lwd=2)

t <- rpart::rpart(form, data=data.frame(fcs.data),method=list(eval=rpart.eval, split=rpart.split, init=rpart.init),control=rpart.control(minsplit=60,minbucket=20,usesurrogate=0))
w <- as.factor(t$where)
levels(w) <- 1:length(unique(w))
print(table(w))
splits <- sort(unique(w))

#t <- mvpart::mvpart( fcs.data[,c('diff.PSTAT5.4','diff.PSTAT5.3','diff.PSTAT5.2')] ~ FSCW + SSCW + FSCH + SSCH + FSCA + SSCA + CD4 + CD25 + CD45RA + FOXP3, data=data.frame(fcs.data),method=list(eval=rpart.eval, split=rpart.split, init=rpart.init),control=rpart.control(minsplit=60,minbucket=20,usesurrogate=0))

par(mfrow=c(3,4))
for (i in sort(unique(t$where))) {
    x <- fcs.data[which(t$where==i),]
    plot(normalised.density(x[,'diff.PSTAT5.2']), main=i, xlim=c(-.5,2), col=1)
    lines(normalised.density(x[,'diff.PSTAT5.3']), col=2)
    lines(normalised.density(x[,'diff.PSTAT5.4']), col=3)
}

par(mfrow=c(1,1))
#t$where <- w
plot(as.party(t))




recursive.part.plot.density <- function(x) {
    if (is.matrix(x)) {
        y <- x[,'diff.PSTAT5.4']
        plot(density(y), main=paste(table(y > .5),collapse=':'))
        abline(v=.5)
    } else if (is.list(x)) {
        for (i in 1:length(x)) recursive.part.plot.density(x[[i]])
    } else {
        print(class(x))
    }
} 
recursive.part.plot.density(X)


# The 'evaluation' function.  Called once per node.
#  Produce a label (1 or more elements long) for labeling each node,
#  and a deviance.  The latter is
#	- of length 1
#       - equal to 0 if the node is "pure" in some sense (unsplittable)
#       - does not need to be a deviance: any measure that gets larger
#            as the node is less acceptable is fine.
#       - the measure underlies cost-complexity pruning, however
rpart.eval <- function(y, wt, parms) {
    #wmean <- sum(y*wt)/sum(wt)
    #rss <- sum(wt*(y-wmean)^2)
    #list(label= wmean, deviance=rss)
    print(t <- table(y < .5))
    list(label=paste(as.character(t),collapse=':'), deviance=min(t))
    #list(label=c(t[[1]],t[[2]]), deviance=min(prop.table(t)))
}

# The split function, where most of the work occurs.
#   Called once per split variable per node.
#   The actual x variable is ordered
#   y is supplied in the sort order of x, with no missings,
#   return two vectors of length (n-1):
#      goodness = goodness of the split, larger numbers are better.
#                 0 = couldn't find any worthwhile split
#        the ith value of goodness evaluates splitting obs 1:i vs (i+1):n
#      direction= -1 = send "y< cutpoint" to the left side of the tree
#                  1 = send "y< cutpoint" to the right
#         this is not a big deal, but making larger "mean y's" move towards
#         the right of the tree, as we do here, seems to make it easier to
#         read
# The reason for returning a vector of goodness is that the C routine
#   enforces the "minbucket" constraint. It selects the best return value
#   that is not too close to an edge.
rpart.split <- function(y, wt, x, parms, continuous) {
    #browser()
    ## Center y
    #browser()
    n <- length(y)
    #y <- y- sum(y*wt)/sum(wt)
	## continuous x variable
	#temp <- cumsum(y*wt)[-n]
	#left.wt  <- cumsum(wt)[-n]
	#right.wt <- sum(wt) - left.wt
	#lmean <- temp/left.wt
	#rmean <- -temp/right.wt
	#goodness <- (left.wt*lmean^2 + right.wt*rmean^2)/sum(wt*y^2)
	#list(goodness= goodness, direction=sign(lmean))
    #v <- cumsum(x**2)/(1:length(x))- (cumsum(x)**2) / (1:length(x))**2
    #v2 <- (sum(x**2)-cumsum(x**2))/(length(x):1) - ((sum(x)-cumsum(x))**2) / (length(x):1)**2
    l <- list(goodness=otsu(x)/(y[1:(n-1)]-.5)**2, direction=rep(-1,n-1))
    cat(max(l$goodness), x[which.max(l$goodness)], y[which.max(l$goodness)], '\n')
    return(l)
}

# The init function:
#   fix up y to deal with offsets
#   return a dummy parms list
#   numresp is the number of values produced by the eval routine's "label"
#   numy is the number of columns for y
#   summary is a function used to print one line in summary.rpart
# In general, this function would also check for bad data, see rpart.poisson
#   for instace.
rpart.init <- function(y, offset, parms, wt) {
    if (!is.null(offset)) y <- y-offset
    list(y=y, parms=0, numresp=1, numy=1,
         summary=function(yval, dev, wt, ylevel, digits ) {
		  paste("  mean=", format(signif(yval, digits)), ", MSE=" , format(signif(dev/wt, digits)), sep='')
	     })
    }


with(as.data.frame(fcs.data), length(which(CD4>=1.849 & CD4< 3.357 & CD25>=0.9651 & CD25< 1.873 & CD45RA>=1.183 & CD45RA< 1.528 & CD4>=2.16 & CD25>=1.378)))


plot(density(fcs.data[,'CD4']))
abline(v=c(1.849, 2.16, 3.357))

plotClusters(fcs.data[,c('SSCA','FSCA','CD4','CD25','CD45RA')], classification=with(as.data.frame(fcs.data),as.numeric(as.logical(CD4>=1.849 & CD4< 3.357 & CD25>=0.9651 & CD25< 1.873 & CD45RA>=1.183 & CD45RA< 1.528 & CD4>=2.16 & CD25>=1.378))))


  CD4>=1.849
   CD4< 3.357
   CD25>=0.9651
   CD25< 1.873
   CD45RA>=1.183
   CD45RA< 1.528
   CD4< 2.16

plotClusters(fcs.data[,c('SSCA','FSCA','CD4','CD25','CD45RA')], classification=with(as.data.frame(fcs.data),as.numeric(as.logical( CD4>=1.849 & CD4< 3.357 & CD25>=0.9651 & CD25< 1.873 & CD45RA>=1.183 & CD45RA< 1.528 & CD4>=2.16 & CD25< 1.378 & SSCA< 3.754))))


f <- as.data.frame(t$frame)
pdf(file.path(opt$out.dir,gsub('.RData','.pdf',basename(file.name))))
plotClusters(fcs.data[,c('diff.PSTAT5.4', as.character(unique(f$var[which(f$splits[,'cutleft']!='')])))], classification=w, ellipses=FALSE, chull.lwd=2)
dev.off()

#plotClusters(d[,c('diff.PSTAT5.4','CD45RA','CD25')],classification=t$where, ellipses=FALSE, chull.lwd=2) 
#plotClusters(d[,c('diff.PSTAT5.4','CD45RA','FOXP3','CD25')],classification=prune.tree(t,best=3)$where, ellipses=FALSE, chull.lwd=2)

w <- as.factor(t$where)
levels(w) <- 1:length(unique(w))
table(w)


par(mfrow=c(2,3))
for (i in sort(unique(w))) {
smoothPlot(fcs.data[which(w==i),c('CD4','SSCA')])
}
for (i in sort(unique(w))) {
plot(normalised.density(fcs.data[which(w==i),'diff.PSTAT5.4']))
}

d <- fcs.data[which(w==2),]

t <- tree::tree('diff.PSTAT5.4 ~ CD4 + CD45RA + CD25 + FOXP3', data.frame(d), split='deviance')
t <- tree::prune.tree(t,best=20)
f <- as.data.frame(t$frame) 
w <- as.factor(t$where)
levels(w) <- 1:length(unique(w))
table(w)

par(mfrow=c(2,4))
for (i in sort(unique(w))) {
    smoothPlot(d[which(w==i),c('CD25','FOXP3')],main=i)
}
for (i in sort(unique(w))) {
plot(normalised.density(d[which(w==i),'diff.PSTAT5.4']), main=i)
}



plotClusters(d[,c('diff.PSTAT5.4', as.character(unique(f$var[which(f$splits[,'cutleft']!='')])))], classification=t$where, ellipses=FALSE, chull.lwd=2)

plotClusters(d[,c('diff.PSTAT5.4', 'CD4', 'CD45RA')], classification=w, ellipses=FALSE, chull.lwd=2)


by(fcs.data, t$where, function
fcs.data[t$where





