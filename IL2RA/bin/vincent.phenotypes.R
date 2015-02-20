
read.csv("/home/vplagnol/export/FACS/rawData/Calli_CD25bright_CBR200.csv")->d

#memory
colnames(d)[grep('Lymphocytes.CD4.CD127hi.127hi.CD45RAneg.ct',colnames(d))] <- 'memory.count'
colnames(d)[grep('Lymphocytes.CD4.CD127hi.127hi.CD45RAneg.freqPar',colnames(d))] <- 'memory.freqpar'
colnames(d)[grep('Lymphocytes.CD4.CD127hi.127hi.CD45RAneg.MeanAPC', colnames(d))] <- 'APC.mfi'
colnames(d)[grep('Lymphocytes.CD4.CD127hiCD25pos.CD45RAneg.MeanAPC_norm', colnames(d))] <- 'APC.mef'
#naive count
colnames(d)[grep('Lymphocytes.CD4.CD127hi.127hi.CD45RApos.ct', colnames(d))] <- 'naive.count'
#naive cd25 pos count
colnames(d)[grep( 'Lymphocytes.CD4.CD127hiCD25pos.CD45RApos.ct', colnames(d))] <- 'naive.cd25pos.count'
#naive cd25 neg count
colnames(d)[grep( 'Lymphocytes.CD4.CD127hiCD25neg.CD45RApos.ct', colnames(d))] <- 'naive.cd25neg.count'
#naive cd25 pos freq
colnames(d)[grep('CD4.127hi.CD45RApos.CD25pos.freqPar', colnames(d))] <- 'naive.cd25pos.freqpar'

d[order(d$individual,d$date),]->d
#d[!duplicated(d$individual),]->d

#bead data necessary for normalisation to MEF
read.csv("~nikolas/dunwich/Beads.Stats/beads.kmeans.fcs2.stats")->beads.kmeans.fcs2.stats
read.csv("~nikolas/dunwich/Beads.Stats/beads.pam.fcs2.stats")->beads.pam.fcs2.stats
#upper case all filenames and discard .FCS extensions
unlist(strsplit(toupper(beads.pam.fcs2.stats$file.name),".FCS"))->beads.pam.fcs2.stats$file.name
read.csv("~nikolas/dunwich/Beads.Stats/beads.pam.fcs3.stats")->beads.pam.fcs3.stats
#upper case all filenames and discard .FCS extensions
unlist(strsplit(toupper(beads.pam.fcs3.stats$file.name),".FCS"))->beads.pam.fcs3.stats$file.name
read.csv("~nikolas/dunwich/Beads.Stats/beads.manual.fcs2.stats")->beads.manual.fcs2.stats
read.csv("~nikolas/dunwich/Beads.Stats/manual.fcs2.stats")->manual.fcs2.stats




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

#compute.mef(d, beads.manual.fcs2.stats) -> d.mef
#merge(d[,-grep("APC.mef", names(d))], d.mef, by=c('individual', 'date', 'APC.mfi'))->d

d$APC.slope*d$APC.mfi->d$APC.mef


d[,c("date","individual","fcsFile", "APC.mfi", "APC.mef", "memory.freqpar", "memory.count", "naive.count", "naive.cd25pos.count", "naive.cd25neg.count", "naive.cd25pos.freqpar")]->vincent.phenotypes

na.omit(vincent.phenotypes) -> vincent.phenotypes
gsub(".fcs$", "", tolower(vincent.phenotypes$fcsFile))->vincent.phenotypes$fcsFile

write.csv(vincent.phenotypes, quote=F, row.names=F, file="~nikolas/CellPhenotypes/vincent.cell.phenotypes")

