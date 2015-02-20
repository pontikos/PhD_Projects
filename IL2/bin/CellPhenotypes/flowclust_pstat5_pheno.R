#!/usr/bin/env Rscript
source('~nikolas/bin/FCS/fcs.R')
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tools"))


option_list <- list( 
make_option(c("--in.file"), default='', help = 'flowClust result file'),
make_option(c("--out.file"), default='', help = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser) 
in.file <- opt$in.file

l <- unlist(strsplit(basename(file_path_sans_ext(in.file)), '_'))
print(individual <- l[[1]])
print(Date <- as.Date(l[[2]]))

#clusters
cluster.file <- in.file
load(cluster.file)
colnames(result$res@mu) <- result$res@varNames

out.file <- opt$out.file

#fcs.name comes from loading RData
print(result$fcs.name)
print(head(x <- as.data.frame(result$fcs.data)))

x$PSTAT5.2.diff <- x$PSTAT5.2-x$PSTAT5.1
x$PSTAT5.3.diff <- x$PSTAT5.3-x$PSTAT5.1
x$PSTAT5.4.diff <- x$PSTAT5.4-x$PSTAT5.1

print(varNames <- result$prior$res@varNames)
print(sum(table(clusters <- MAP(result$fcs.data[,varNames], result$res, .95))))
print(clust <- sort(unique(clusters)))


print(
d <- cbind(data.frame(result$res@mu[clust,]),
           data.frame(weights=result$res@w[clust]),
           data.frame(
            individual=individual,
            date=Date,
            cluster=clust,
            pheno.pstat5.median.diff.2=sapply(clust, function(pop) median(x[which(clusters==pop),'PSTAT5.2.diff'])),
            pheno.pstat5.median.diff.3=sapply(clust, function(pop) median(x[which(clusters==pop),'PSTAT5.3.diff'])),
            pheno.pstat5.median.diff.4=sapply(clust, function(pop) median(x[which(clusters==pop),'PSTAT5.4.diff'])),
            #pheno.pstat5.mean.diff=sapply(clust, function(pop) mean(x[which(clusters==pop),'PSTAT5.4']-x[which(clusters==pop),'PSTAT5.1'])),
            #pheno.pstat5.mean.diff=sapply(clust, function(pop) mean(x[which(clusters==pop),'PSTAT5.3']-x[which(clusters==pop),'PSTAT5.1'])),
            #pheno.pstat5.var.diff=sapply(clust, function(pop) var(x[which(clusters==pop),'PSTAT5.4']-x[which(clusters==pop),'PSTAT5.1'])),
            #pheno.pstat5.sd.diff=sapply(clust, function(pop) sd(x[which(clusters==pop),'PSTAT5.4']-x[which(clusters==pop),'PSTAT5.1'])),
            #pheno.pstat5.sd.diff=sapply(clust, function(pop) sd(x[which(clusters==pop),'PSTAT5.3']-x[which(clusters==pop),'PSTAT5.1'])),
            #pheno.pstat5.diff.mean=sapply(clust, function(pop) mean(x[which(clusters==pop),'PSTAT5.4'])-mean(x[which(clusters==pop),'PSTAT5.1'])),
            #pheno.pstat5.median.diff=sapply(clust, function(pop) median(x[which(clusters==pop),'PSTAT5.4']-x[which(clusters==pop),'PSTAT5.1'])),
            #pheno.pstat5.diff.median=sapply(clust, function(pop) median(x[which(clusters==pop),'PSTAT5.4'])-median(x[which(clusters==pop),'PSTAT5.1'])),
            pheno.median.pstat5.1=sapply(clust, function(pop) median(x[which(clusters==pop),'PSTAT5.1'])),
            pheno.median.pstat5.2=sapply(clust, function(pop) median(x[which(clusters==pop),'PSTAT5.2'])),
            pheno.median.pstat5.3=sapply(clust, function(pop) median(x[which(clusters==pop),'PSTAT5.3'])),
            pheno.median.pstat5.4=sapply(clust, function(pop) median(x[which(clusters==pop),'PSTAT5.4']))
            )
)
)

write.csv(d, file=out.file, quote=FALSE, row.names=FALSE)


