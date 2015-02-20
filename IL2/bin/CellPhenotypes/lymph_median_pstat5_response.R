#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

option_list <- list( 
make_option(c("--base.dir"), default='~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/', help=''),
make_option(c("--in.dir"), default='Lymphocytes/All', help = "in.dir [default %default]"),
make_option(c("--out.file"), default='', help = "")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

base.dir <- '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/'
base.dir <- opt$base.dir

#join BEFORE CD4 lymph clustering
in.dir <- "~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/Lymphocytes6/RData/"
#join AFTER CD4 lymph clustering
in.dir <- "~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/Lymphocytes/"
in.dir <- opt$in.dir

#print(out.file <- file.path(base.dir,'CellPhenotypes',paste('lymph_median',ext,'response.csv',sep='_')))
print(out.file <- opt$out.file)

individual.date <- read.csv(file.path(base.dir,'individual-date.csv'), col.names=c('individual','date'))
X <- NULL
for (i in 1:nrow(individual.date)) {
    individual <- individual.date[i,'individual']
    date <- individual.date[i,'date']
    print( f <- sprintf('%s.RData',file.path(in.dir, paste(individual, date, sep='_'))) )
    print(load(f))
    print(head(x <- fcs.data))
    pstat5 <- paste('PSTAT5',1:4,sep='.')
    #X <- rbind( X, data.frame( individual=individual, date=date, var=pstat5, mean=colMeans(x[,pstat5]), row.names=NULL ) )
    X <- rbind( X, data.frame( individual=individual, date=date, var=paste('pSTAT5.diff',1:4,sep='.'), mean=sapply(1:4, function(i) mean(x[,paste('PSTAT5',i,sep='.')]-x[,'PSTAT5.1'])), row.names=NULL ) )
}
X$day <- NA
for (ind in X$individual) {
  cat(ind, X[which(X$individual==ind),'day'] <- as.factor(as.character(X[which(X$individual==ind),'date'])), '\n')
}
X.rep <- X[which(X$individual %in% unique(X[which(X$day==2),'individual'])),] 
X.rep <- cbind(
X.rep[order(X.rep$individual,X.rep$var),][c(TRUE,FALSE),c('individual','var')],
day.1=X.rep[order(X.rep$individual,X.rep$var),][c(TRUE,FALSE),'mean'], day.2=X.rep[order(X.rep$individual,X.rep$var),][c(FALSE,TRUE),'mean']
)
print(X.rep)

df.cor <- as.data.frame(do.call('rbind',as.list(by(X.rep,X.rep$var,function(x) sprintf("r^2==%.2f", cor(x$day.1,x$day.2)**2)))))
df.cor$var <- rownames(df.cor)
print(df.cor)

library(ggplot2)
ggplot(X.rep, aes(x=day.1, y=day.2)) + facet_wrap( ~ var) + geom_point() + xlab('Day 1') + ylab('Day 2') + geom_abline(intercept=0,slope=1) + geom_text(data=df.cor, aes(x=0, y=3, label=V1), parse=TRUE)




colnames(X) <- c('individual','date',paste(gsub('-','.',ext),'response',sep='.'))

write.csv(X,file=out.file,quote=FALSE,row.names=FALSE)



