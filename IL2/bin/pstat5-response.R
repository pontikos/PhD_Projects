#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("mclust"))

option_list <- list( 
make_option(c("--in.dir"), help = "in.dir"),
make_option(c("--plot.dir"), help = "plot.dir"),
make_option(c("--individual"), help = "individual"),
make_option(c("--date"), help = "date"),
make_option(c("--ext"), default='', help = "pstat5"),
make_option(c("-p","--parameter"), default='pSTAT5', help=" [default %default]"),
make_option(c("--doses"), default='0U,01U,10U,1000U', help=" [default %default]")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

print(in.dir <- opt$in.dir)
print(out.dir <- opt$plot.dir)
print(individual <- opt$individual)
print(date <- opt$date)
print(p <- opt$parameter)
print(ext <- opt$ext)
print(doses <- unlist(strsplit(opt$doses, ",")))

d <- list()
for (dose in doses) {
    print(f <- file.path(in.dir, sprintf("%s.csv",paste(individual, dose, date, ext, sep='_'))))
    x <- read.csv(f)
    print(colnames(x) <- toupper(colnames(x)))
    x <- as.data.frame(x)
    x$pop <- NA
    cd45ra <- x[,'CD45RA']
    m <- mclust::Mclust(cd45ra, G=2, modelNames='V')
    g.cd45ra <- cd45ra[which.max(m$uncertainty)]
    # gate naive
    x <- x[cd45ra > g.cd45ra,]
    # gate memory
    x <- x[cd45ra < g.cd45ra,]
    d[[dose]] <- x[,toupper(p)]
}

cdf <- lapply(d, ecdf) 
style <- data.frame(conc=doses, col=c('black', 'blue', 'green', 'red'), lwd=seq(.5, by=.75, length=4), stringsAsFactors=FALSE)

x <- seq(min(unlist(d)), max(unlist(d)), .1)
ecdf.diff <- function(x) return(cdf[['1000U']](x)-cdf[['0U']](x))
area <- sum(abs(ecdf.diff(x)))
cat('>>', individual, date, area, '\n', sep=',')

print(plot.file <- file.path(out.dir, sprintf("%s.png",paste(individual, 'cdf', date, ext, sep='_'))))
png(plot.file)
plot(cdf[['0U']], main=paste(individual, date, ext), xlab=ext, ylab='', lwd=0, col='white')
sapply(1:length(d), function(i) lines( cdf[[style$conc[i]]], lwd=style$lwd[i], col=style$col[i])) 
lines(x, abs(ecdf.diff(x)), lty=1, lwd=.25) 
legend(x="bottomright", legend=style$conc, col=style$col, bty='n')
dev.off()




