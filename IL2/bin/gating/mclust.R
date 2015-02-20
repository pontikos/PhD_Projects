#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("mclust"))
source('~nikolas/bin/FCS/fcs.R')

option_list <- list( 
make_option(c("--in.dir"), help = "in.dir"),
make_option(c("--out.dir"), help = "out.dir"),
make_option(c("--individual"), help = "individual"),
make_option(c("--date"), help = "date"),
make_option(c("--ext"), default='', help = "pstat5"),
make_option(c("-p","--parameters"), default='CD25,CD45RA,CD4,FOXP3', help=" [default %default]"),
make_option(c("--dose"), default=NULL, help=" "),
make_option(c("--predict"), default=NULL, help="Mclust object used for prediction")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

print(p <- unlist(strsplit(opt$parameters, ",")))
print(ext <- opt$ext)

# read from Lymphocytes/All/
print(f <- file.path(opt$in.dir, sprintf("%s.csv",paste(opt$individual, opt$date, opt$ext, sep='_'))))
x <- read.csv(f)

#print(rdata.file <- file.path(opt$out.dir, 'RData', sprintf("%s.RData",paste(opt$individual, opt$dose, opt$date, opt$ext, sep='_'))))

if (!is.null(opt$predict)) {
    print(load(opt$predict))
    res <- predict(m, x[,p])
    #plot
    #print(plot.file <- file.path(opt$out.dir, 'Plot', sprintf("%s.png",paste(opt$individual, opt$dose, opt$date, opt$ext, sep='_'))))
    #png(plot.file)
    #plot(res$classification)
    #dev.off()
    print(head(X <- cbind(x, cluster=res$classification)))
    # output file
    print(f <- file.path(opt$out.dir, 'Csv', sprintf("%s.csv",paste(opt$individual, opt$date, opt$ext, sep='_'))))
    write.csv(X, file=f, quote=FALSE, row.names=FALSE)
} else {
    #initialisation uses hierarchical clustering so we need to provide a subset of the data
    res <- Mclust(x[,p], initialization=list(subset=1:1000))
    #plot
    print(plot.file <- file.path(opt$out.dir, 'Plot', sprintf("%s.pdf",paste(opt$individual, opt$dose, opt$date, opt$ext, sep='_'))))
    pdf(plot.file)
    plot(res)
    dev.off()
    print(plot.file <- file.path(opt$out.dir, 'Plot', sprintf("%s.png",paste(opt$individual, opt$dose, opt$date, opt$ext, sep='_'))))
    png(plot.file)
    #plot(res)
    plot.clusters(x[,p], clusters=res$classification)
    dev.off()
}

save(res, file=file.path(opt$out.dir, 'RData', sprintf("%s.RData",paste(opt$individual, opt$dose, opt$date, opt$ext, sep='_'))))




