#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

option_list <- list( 
make_option(c("-f","--fcs"), help = "fcs file")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowCore, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(cluster, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(mixtools, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))

read.FCS(opt$fcs)[,-9] -> fcs.data
c(1,2,3,5,7,8) -> p
fcs.data[rowSums(fcs.data@exprs[,p]>1)==dim(fcs.data@exprs[,p])[[2]],]->fcs.data

exprs <- function(fcs.data) {
    fcs.data@exprs -> x
    log10(x[,3:8]) -> x[,3:8]
    return(x)
}

### auto gating strategy
exprs(fcs.data)->x

suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(modeest, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))

# pe < 2
# cd4 > 2
# cd25 < 2
# cd45ra < 2
# naive: cd45ra > 1
# memory: cd45ra < 1

#scatter gate
x[x[,2] < 200 & x[,1] < 2000,] -> x

x[,6] -> pe
#x[pe < mfv(pe[pe>0])-var(pe[pe>0]),]->x

#CD4 gate
x[,7] -> cd4
#x[mfv(cd4)-var(cd4) < cd4 & cd4 < mfv(cd4)+var(cd4),]->x
#fixed gate
x[cd4>2 & pe<2,]->x

#non tregs gate
x[,5] -> apc
x[,3] -> cd127
x[,8] -> cd45ra

#x[apc < mfv(apc)+2*var(apc) & cd127 > mfv(cd127)-3*var(cd127),]->x
#x[apc < 2 & cd127 > mfv(cd127)-3*var(cd127),]->x

#fixed gates
x[apc < 2 & cd45ra < 2,]->x


pam(x[,8], k=2)->res

sort(res$medoids)->medoids
print(medoids)


### parametric mixture model
normalmixEM(x[,8], mu=medoids) -> mm
#save(mm, file=sprintf('mm.%s.obj',opt$fcs))

x[mm$posterior[,1] > .95,]  -> memory
x[mm$posterior[,2] > .95,] -> naive

mean(memory[,5]) -> cd25pos.threshold
memory[memory[,5]>=cd25pos.threshold,] -> memory.cd25pos
naive[naive[,5]>=cd25pos.threshold,] -> naive.cd25pos

pdf('rplots.pdf')

plot(x[,c(5,8)], pch=".")
abline(v=cd25pos.threshold)

dev.off()


print(dim(x))
print(100*dim(naive.cd25pos)[[1]]/dim(naive)[[1]])
print(100*dim(memory.cd25pos)[[1]]/dim(x)[[1]])
print(100*dim(naive)[[1]]/dim(x)[[1]])
print(100*dim(memory)[[1]]/dim(x)[[1]])
print(mean(10**memory[,5]))
print(median(10**memory[,5]))
print(mean(10**memory.cd25pos[,5]))

dim(x)[[1]] -> nontregs.count
dim(memory)[[1]] -> memory.count
cat('>>fcsFile,APC.mfi,memory.freqpar,memory.count\n')
cat('>>')
cat(opt$fcs, mean(10**memory[,5]), 100*memory.count/nontregs.count, memory.count, sep=",")
cat('\n')

