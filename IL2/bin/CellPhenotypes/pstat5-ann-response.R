#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("RANN"))

option_list <- list( 
make_option(c("--in.dir"), help = "in.dir"),
make_option(c("--plot.dir"), help = "plot.dir"),
make_option(c("--out.dir"), help = "out.dir"),
make_option(c("--individual"), help = "individual"),
make_option(c("--date"), help = "date"),
make_option(c("--ext"), default='', help = "pstat5"),
make_option(c("-p","--parameter"), default='pSTAT5', help=" [default %default]"),
make_option(c("--doses"), default='0U,01U,10U,1000U', help=" [default %default]")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

in.dir <- '~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/Lymphocytes/RData'
date <- '2012-11-21'
individual <- 'CB01509P'
doses <- c('0U','01U','10U','1000U')
ext <- 'pstat5'

print(in.dir <- opt$in.dir)
print(out.dir <- opt$out.dir)
print(individual <- opt$individual)
print(date <- opt$date)
print(p <- opt$parameter)
print(ext <- opt$ext)
print(doses <- unlist(strsplit(opt$doses, ",")))

d <- list()
for (dose in doses) {
    print(f <- file.path(in.dir, sprintf("%s.csv",paste(individual, dose, date, ext, sep='_'))))
    x <- read.csv(f)
    print(p <- colnames(x) <- toupper(colnames(x)))
    d[[dose]] <- x[,toupper(p)]
}


# nearest neighbour (excluding pstat5 dimension) from 0U to 0U, 01U, 10U, and 1000U 
# returns dataset pstat5.diff of dim 0U containing pstat5 difference relative to 0U
p <- colnames(d[[1]])[colnames(d[[1]])!='PSTAT5']
pstat5.diff <- numeric(dim(d[['0U']])[[1]])
for (dose in doses) {
    nn <- nn2(d[[dose]][,p],query=d[['0U']][,p],k=1)
    pstat5.diff <- cbind(pstat5.diff, d[[dose]][nn$nn.idx,'PSTAT5']-d[['0U']][,'PSTAT5'])
}
pstat5.diff <- pstat5.diff[,-1]

colnames(pstat5.diff) <- doses

#write file of pstat5.diff to csv file
print(out.file <- file.path(out.dir, sprintf("%s.csv",paste(individual, date, ext, sep='_'))))

write.csv(pstat5.diff, file=out.file, quote=FALSE, row.names=FALSE)



#which(rowSums(pstat5.diff)>6)
#x <- seq(min(unlist(d)), max(unlist(d)), .1)
#ecdf.diff <- function(x) return(cdf[['1000U']](x)-cdf[['0U']](x))
#area <- sum(abs(ecdf.diff(x)))

#cat('>>', individual, date, mean(rowSums(pstat5.diff)), '\n', sep=',')
cat('>>', individual, date, sum(pstat5.diff), '\n', sep=',')


#
#print(plot.file <- file.path(out.dir, sprintf("%s.png",paste(individual, 'pstat5-ann', date, ext, sep='_'))))
#png(plot.file)
#plot([['0U']], main=paste(individual, date, ext), xlab=ext, ylab='', lwd=0, col='white')
#sapply(1:length(d), function(i) lines( [[style$conc[i]]], lwd=style$lwd[i], col=style$col[i])) 
#lines(x, abs(ecdf.diff(x)), lty=1, lwd=.25) 
#legend(x="bottomright", legend=style$conc, col=style$col, bty='n')
#dev.off()
#






