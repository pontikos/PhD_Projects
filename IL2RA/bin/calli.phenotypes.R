read.csv("~nikolas/CellPhenotypes/Calli.Gating/nonTregs.stats", header=T) -> nonTregs
nonTregs$experiment.date->nonTregs$date
nonTregs$file.name->nonTregs$fcsFile
read.csv("~nikolas/CellPhenotypes/Calli.Gating/memory.stats", header=T) -> memory
memory$experiment.date->memory$date
memory$file.name->memory$fcsFile
read.csv("~nikolas/CellPhenotypes/Calli.Gating/naive.stats", header=T) -> naive
naive$experiment.date->naive$date
naive$file.name->naive$fcsFile
read.csv("~nikolas/CellPhenotypes/Calli.Gating/naive.cd25neg.stats", header=T) -> naive.cd25neg
naive.cd25neg$experiment.date->naive.cd25neg$date
naive.cd25neg$file.name->naive.cd25neg$fcsFile
read.csv("~nikolas/CellPhenotypes/Calli.Gating/naive.cd25pos.stats", header=T) -> naive.cd25pos
naive.cd25pos$experiment.date->naive.cd25pos$date
naive.cd25pos$file.name->naive.cd25pos$fcsFile


(merge(nonTregs, memory, by=c('individual', 'date', 'fcsFile'), suffixes=c('.nonTregs','.memory'))) -> memory
memory$APC.length.memory->memory$memory.count
memory$APC.length.nonTregs->memory$nonTregs.count
100*memory$memory.count/memory$nonTregs.count -> memory$memory.freqpar
memory$APC.mean.memory->memory$APC.mfi
memory[,c("date","individual","fcsFile", "APC.mfi", "memory.freqpar", "memory.count")]->memory


(merge(nonTregs, naive, by=c('individual', 'date', 'fcsFile'), suffixes=c('.nonTregs','.naive'))) -> naive
naive$APC.length.naive->naive$naive.count
naive$APC.length.nonTregs->naive$nonTregs.count
100*naive$naive.count/naive$nonTregs.count -> naive$naive.freqpar
naive[,c("date","individual","fcsFile", "naive.freqpar", "naive.count")]->naive



(merge(naive, naive.cd25neg, by=c('individual', 'date', 'fcsFile'), suffixes=c('.naive','.naive.cd25.neg'))) -> naive.cd25neg
naive.cd25neg$APC.length->naive.cd25neg$naive.cd25neg.count
100*naive.cd25neg$naive.cd25neg.count/naive.cd25neg$naive.count -> naive.cd25neg$naive.cd25neg.freqpar
naive.cd25neg[,c("date","individual","fcsFile", "naive.cd25neg.freqpar", "naive.cd25neg.count")]->naive.cd25neg


(merge(naive, naive.cd25pos, by=c('individual', 'date', 'fcsFile'), suffixes=c('.naive','.naive.cd25.pos'))) -> naive.cd25pos
naive.cd25pos$APC.length->naive.cd25pos$naive.cd25pos.count

### naive.cd25pos.count + naive.cd25neg.count > naive.count !?
#100*naive.cd25pos$naive.cd25pos.count/naive.cd25pos$naive.count -> naive.cd25pos$naive.cd25pos.freqpar
100*naive.cd25pos$naive.cd25pos.count/(naive.cd25pos$naive.cd25pos.count+naive.cd25neg$naive.cd25neg.count) -> naive.cd25pos$naive.cd25pos.freqpar

naive.cd25pos[,c("date","individual","fcsFile", "naive.cd25pos.freqpar", "naive.cd25pos.count")]->naive.cd25pos


merge(memory, naive, by=c('individual', 'date', 'fcsFile'))->m
merge(m, naive.cd25neg, by=c('individual', 'date', 'fcsFile'))->m
merge(m, naive.cd25pos, by=c('individual', 'date', 'fcsFile'))->m
m[order(m$individual,m$date),]->calli

### Calculate mef based on coefficients for that day from regression on bead stats.
### Two types of regression: ours and vincent's
compute.mef <- function(memory.apc.mfi, beads.stats, mef=c(4100, 10300, 25500, 67300, 139100), base=exp(1)) {

    merge(memory.apc.mfi, beads.stats, by="date") -> memory.apc.mfi
    #normalisation: log(MEF, base) = alpha x log(MFI, base) + beta
    if (is.numeric(base)) {
        log(mef, b=base)->mef
        fun <- function (d) {
            log(as.numeric((d[paste("mfi",2:6,sep=".")])), b=base)->mfi
            lm(mef ~ mfi)->m
            return(c(alpha=m$coefficients[2], beta=m$coefficients[1], rse=summary(m)$sigma))
        }
        t(apply(memory.apc.mfi, 1, fun))->d
        colnames(d)<-c("alpha","beta","rse")
        cbind(memory.apc.mfi, d)->memory.apc.mfi
        base^(memory.apc.mfi$beta)*(memory.apc.mfi$APC.mfi**memory.apc.mfi$alpha) -> memory.apc.mfi$APC.mef
        return(memory.apc.mfi[,c("date", "individual", "APC.mfi", "APC.mef")])
    #vincent normalisation: MEF = alpha x MFI
    } else {
        fun <- function (d) {
            as.numeric((d[paste("mfi",2:6,sep=".")]))->mfi
            lm(mef ~ mfi - 1)->m
            return(c(alpha=m$coefficients[1], rse=summary(m)$sigma))
        }
        t(apply(memory.apc.mfi, 1, fun))->d
        colnames(d)<-c("alpha","rse")
        cbind(memory.apc.mfi, d)->memory.apc.mfi
        memory.apc.mfi$APC.mfi*memory.apc.mfi$alpha -> memory.apc.mfi$APC.mef
        return(memory.apc.mfi[,c("date", "individual", "APC.mfi", "APC.mef")])
    }
}

#bead data necessary for normalisation to MEF
read.csv("~nikolas/dunwich/Beads.Stats/beads.pam.fcs2.stats") -> beads.stats
#unlist(strsplit(toupper(beads.stats$file.name),".FCS"))->beads.stats$file.name
tolower(beads.stats$file.name) -> beads.stats$file.name
gsub('(cad.*_treg).*', '\\1' ,beads.stats$file.name)->beads.stats$experiment
as.Date(beads.stats$date) -> beads.stats$date

gsub('(cad.*_treg).*', '\\1' ,calli$fcsFile)->calli$experiment
#leave out dates
merge(beads.stats[,-2], calli, by=c('experiment'))->calli

#compute MEF
exp(calli$mef.beta)*(calli$APC.mfi**calli$mef.alpha) -> calli$APC.mef

#merge(m, mef, by=c('individual', 'date', 'APC.mfi'))->m

#use vincent normalisation
#read.csv("/home/vplagnol/export/FACS/rawData/Calli_CD25bright_CBR200.csv")->d
#gsub(".fcs$", "", tolower(d$fcsFile))->d$fcsFile
#d[,c("fcsFile", "individual", "APC.slope")]->d
#merge(m, d, by=c('individual', 'fcsFile'))->m
#m$APC.slope*m$APC.mfi->m$APC.mef

calli[,c("date","individual","fcsFile", "APC.mfi", "APC.mef", "memory.freqpar", "memory.count", "naive.count", "naive.cd25pos.count", "naive.cd25neg.count", "naive.cd25pos.freqpar")]   -> calli.phenotypes

write.csv(calli.phenotypes, quote=F, row.names=F, file="~nikolas/CellPhenotypes/calli.cell.phenotypes")
