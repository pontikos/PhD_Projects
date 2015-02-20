#!/usr/bin/env Rscript
source('~nikolas/bin/FCS/fcs.R',chdir=T)
suppressPackageStartupMessages(library("optparse"))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(tools, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))


### MAIN
#base.dir <- '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/Lymphocytes-test2/pstat5-ann/RData/'
base.dir <- '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/'

option_list <- list( 
    make_option(c("-f","--in.file"), default=file.path(base.dir,'individual-date.csv'), help = 'File containing list of individuals and dates.'),
    make_option(c('--plot.file'), default=NULL, help = 'File to which to plot the clustering results.'),
    #make_option(c("--ext"), default='pstat5', help = "pstat5"),
    make_option(c('--parameter'), help='Channel on which to plot density')
)

base.dir <- '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/Lymphocytes5/RData/'

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

print(opt$in.file)
# Two modes of operation: either compute prior or run on all files with specified prior
# If no prior file specified
# make a prior from a pool of fcs files
# otherwise load prior from specified file
head(d <- read.csv(opt$in.file, col.names=c('individual','date')))
dens <- list()
for (i in 1:nrow(d)) {
    for (dose in c('0U','01U','10U','1000U')) {
        print(f <- file.path(base.dir, sprintf("%s.fcs.csv", paste(d[i,'individual'], dose, d[i,'date'], sep='_'))))
        if (file.exists(f)) {
            print(head(x <- read.csv(f)))
            dens[[f]] <- normalised.density(as.numeric(x[,opt$parameter]))
        } else {
            warning(f, ' does not exist!')
        }
    }
} 


png(opt$plot.file)
plot(dens[[1]], main=opt$parameter)
for (i in 2:length(dens)) lines(dens[[i]])
dev.off()

