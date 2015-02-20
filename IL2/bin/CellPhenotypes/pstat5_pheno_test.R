#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(beeswarm))
suppressPackageStartupMessages(library(nlme))

option_list <- list( 
make_option(c("--in.file"), default='', help = 'flowClust result file'),
make_option(c("--out.file"), default='', help = ""),
make_option(c("--plot.file"), default='', help = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

in.file <- 'flowclust-4-pstat5-response.csv'
in.file <- opt$in.file
head(pstat5 <- read.csv(in.file,stringsAsFactors=FALSE))

head(snps <- read.csv('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/cb-IL2RA.csv',stringsAsFactors=F))
dim(snps.pstat5 <- merge(snps, pstat5))


pdf(opt$plot.file, width=20, height=15)
pstat5.pheno <- do.call( 'rbind', lapply( sort(unique(snps.pstat5$cluster)),
    function(cluster2) {
        l <- list()
        print(l$cluster <- cluster2)
        #mclust mean response
        d <- subset(snps.pstat5,cluster==cluster2)
        print(length(dup.individuals <- names(which(table(d$individual)==2))))
        #dim(d <- merge(snps, merge(pstat5, pstat5norm)))
        length(unique(d$individual)) 
        #dim(d2 <- d[-which(duplicated(d$individual)),])
        d[d=='N-N'] <- NA
        snps <- grep('^rs', colnames(d), value=T)
        snps <- names(which(apply(d[,snps], 2, function(x)!all(is.na(x)))))
        snps <- names(which(apply(d[,snps], 2, function(x) length(table(x))>1)))
        snps <- names(sort(colSums(apply(apply(d[,snps], 2, is.na),2,as.numeric))))
        #apply(d[,snps], 2, table) 
        #table(d$snp)
        d[,snps] <- lapply(d[,snps], factor)
        d[,snps] <- lapply(d[,snps], as.numeric)
        cbind(data.frame(cluster=cluster2),
              do.call('rbind', lapply(c(grep('^pheno', colnames(d), value=TRUE),'CD25'),
               function(p) {
                   l <- list()
                   l$pheno <- p
                    #association with SNPs
                    par(mfrow=c(4,5)) #title(main=cluster2, outer=true)
                    for (snp in snps) {
                        fm <- as.formula(paste(p, snp, sep=' ~ '))
                        #m <- lme(fm,random=~ 1|as.factor(individual), data=d)
                        m <- lme(fm,random=~ 1|as.factor(individual), data=na.omit(d[,c(p,'individual', snp)]))
                        p.value <- summary(m)$tTable[,'p-value'][[2]]
                        l[[paste('p.value',snp,sep='.')]] <- p.value
                        print(main <- tryCatch(sprintf('pvalue=%.3f',p.value), warning=function(e) {''}, finally=''))
                        if (!is.nan(p.value))
                        beeswarm(fm,data=d, main=main, ylim=range(d[,p]),col=cluster2, pch=19, bg=cluster2)
                    }
                    #association with T1D
                    fm <- as.formula(paste(p,'t1d',sep=' ~ '))
                    m <- lme(fm,random=~ 1|as.factor(individual), data=na.omit(d[,c(p,'individual', 't1d')]))
                    p.value <- summary(m)$tTable[,'p-value'][[2]]
                    l[['p.value.t1d']] <- p.value
                    main <- tryCatch(sprintf('pvalue=%.3f', p.value), warning=function(e) {''}, finally='')
                    if (!is.nan(p.value))
                    beeswarm(fm,data=d, main=main, ylim=range(d[,p]),col=cluster2, pch=19, bg=cluster2)
                    #association with sex
                    fm <- as.formula(paste(p,'sex',sep=' ~ '))
                    m <- lme(fm,random=~ 1|as.factor(individual), data=na.omit(d[,c(p,'individual', 'sex')]))
                    p.value <- summary(m)$tTable[,'p-value'][[2]]
                    l[['p.value.sex']] <- p.value
                    main <- tryCatch(sprintf('pvalue=%.3f', p.value), warning=function(e) {''}, finally='')
                    if (!is.nan(p.value))
                    beeswarm(fm,data=d, main=main, ylim=range(d[,p]),col=cluster2, pch=19, bg=cluster2)
                    #repeatability
                    x <- d[ d$individual %in% dup.individuals, p]
                    print( x1 <- x[c(TRUE,FALSE)] )
                    print( x2 <- x[c(FALSE,TRUE)] )
                    l$cor <- r.squared <- cor(x1,x2)**2
                    plot(x1,x2,xlab='day1',ylab='day2',xlim=range(x),ylim=range(x),main=sprintf('cor=%.3f',r.squared),col=cluster2)
                    print(d[9:10,])
                    points(t(d[9:10,p]),pch='x')
                    abline(b=1,a=0)
                    return(as.data.frame(l))
        })))
    } ) )
dev.off()


write.csv(pstat5.pheno, quote=FALSE, row.names=FALSE, file=opt$out.file)



