library(reshape)

CELL.TYPES <- c('Memory Eff','Memory Treg','Naive Eff','Naive Treg')
DOSES <- c('0U','01U','10U','1000U')


base.dir <- '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR'
setwd(base.dir)

counts <- data.frame()
for (f in list.files(base.dir,pattern='.*_.*U_.*.RData')) {
    print( individual.dose.day <- unlist(strsplit(gsub('.RData','',f),'_')) )
    load(f)
    counts <- rbind(counts,
                data.frame( individual=individual.dose.day[1],
                dose=individual.dose.day[2],
                day=individual.dose.day[3],
                t(sapply(CELL.TYPES, function(cell.type) sum(CLR[,cell.type]))) ) 
    )
}

counts <- na.omit(counts)

boxplot.matrix(
cbind(
'Naive Eff'=apply(cast(counts, individual + day ~ dose, value='Naive.Eff')[,DOSES],1,function(x)mean(x)/var(x)),
'Naive Treg'=apply(cast(counts, individual + day ~ dose, value='Naive.Treg')[,DOSES],1,function(x)mean(x)/var(x)),
'Memory Eff'=apply(cast(counts, individual + day ~ dose, value='Memory.Eff')[,DOSES],1,function(x)mean(x)/var(x)),
'Memory Treg'=apply(cast(counts, individual + day ~ dose, value='Memory.Treg')[,DOSES],1,function(x)mean(x)/var(x))
)
)

boxplot.matrix(
cbind(
'Naive Eff'=apply(cast(counts, individual + day ~ dose, value='Naive.Eff')[,DOSES],1,function(x)log(mean(x)/var(x))),
'Naive Treg'=apply(cast(counts, individual + day ~ dose, value='Naive.Treg')[,DOSES],1,function(x)log(mean(x)/var(x))),
'Memory Eff'=apply(cast(counts, individual + day ~ dose, value='Memory.Eff')[,DOSES],1,function(x)log((mean(x)/var(x)))),
'Memory Treg'=apply(cast(counts, individual + day ~ dose, value='Memory.Treg')[,DOSES],1,function(x)log(mean(x)/var(x)))
)
)


