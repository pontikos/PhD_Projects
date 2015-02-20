### LOAD DATA

### get individual info (age, sex, genotype) from vincent file
read.csv("/home/vplagnol/export/FACS/rawData/Calli_CD25bright_CBR200.csv")->vincent.data
#read.csv("~/thor/Calli_CD25bright_CBR200.csv")->vincent.data
gsub(".fcs$", "", tolower(vincent.data$fcsFile))->vincent.data$fcsFile
colnames(vincent.data)[grep('rs12722495', colnames(vincent.data))] <- 'rs12722495'
colnames(vincent.data)[grep('rs2104286', colnames(vincent.data))] <- 'rs2104286'
colnames(vincent.data)[grep('rs11594656', colnames(vincent.data))] <- 'rs11594656'
vincent.data[,c('fcsFile', 'date', 'individual', 'sex', 'age', 'rs12722495', 'rs2104286', 'rs11594656')]->INDIVIDUAL.INFO
colnames(INDIVIDUAL.INFO)<-c('fcsFile', 'date', "individual", "sex", "age", "snp95", "snp86", "snp56")
na.omit(INDIVIDUAL.INFO)->INDIVIDUAL.INFO
#drop individual CB00503W
INDIVIDUAL.INFO[as.character(INDIVIDUAL.INFO$individual)!='CB00503W',]->INDIVIDUAL.INFO
#no bead data for the 2008-05-12
#INDIVIDUAL.INFO[INDIVIDUAL.INFO$date!='2008-05-12',]->INDIVIDUAL.INFO

#remove duplicates
INDIVIDUAL.INFO[!duplicated(as.character(INDIVIDUAL.INFO$individual)),]->UNIQUE.INDIVIDUAL.INFO

#duplicated individuals contains mapping of fcs to individual and date
read.csv("~nikolas/CellPhenotypes/recalled.individuals.pch")  -> RECALLED.INDIVIDUALS.PCH
#read.csv("~/thor/CellPhenotypes/recalled.individuals.pch")  -> RECALLED.INDIVIDUALS.PCH

### Only keep recalled individuals for phenotype repeatability
only.keep.recalled.individuals <- function(all.individuals, recalled.individuals.pch, drop.individuals=c('CB00503W')) {
    #only keep recalled individuals
    merge(all.individuals, recalled.individuals.pch)->recalled.individuals
    recalled.individuals[which(! recalled.individuals$individual %in% drop.individuals ),]->recalled.individuals
    #memory.apc.mfi[which(memory.apc.mfi$individual != 'CB00496N'),]->memory.apc.mfi
    return(recalled.individuals)
}

### Build the recalled table so that for each individual we have phenotype on day 1
### and phenotype on day 2 (.day1 and .day2 suffixes)
recalled.phenotypes.table <- function(phenotype) {
    as.Date(phenotype$date)->phenotype$date
    phenotype->m2
    merge(phenotype, m2, by=c("individual","pch"), suffixes=c(".day1", ".day2"))->recalled.phenotype
    #only keep the rows were day1 earlier than day2 otherwise we have duplicates
    recalled.phenotype[recalled.phenotype$date.day1<recalled.phenotype$date.day2,]->recalled.phenotype
    recalled.phenotype[order(recalled.phenotype$individual),]->recalled.phenotype
    #pch char for plotting
    recalled.phenotype$pch <- as.character(recalled.phenotype$pch)  
    return(recalled.phenotype)
} 

### Make recalled table from cell phenotypes, only need the fcsFile name to do the mapping
make.recalled <- function(cell.phenotypes) {
    return(recalled.phenotypes.table(only.keep.recalled.individuals(cell.phenotypes, recalled.individuals.pch=RECALLED.INDIVIDUALS.PCH)))
}

###
load.individual.cell.phenotypes <- function(phenotype.file) {
    read.csv(phenotype.file)->cell.phenotypes
    #dates in file do not always agree with dates in vincent.data which ignores rows on merge
    #so we ignore date column to avoid complications
    cell.phenotypes[,names(cell.phenotypes)!="date"]->cell.phenotypes
    #merge with information about individual (age, sex, genotype)
    na.omit(merge(INDIVIDUAL.INFO, cell.phenotypes))->individual.cell.phenotypes
    #merge with information about unique individual (age, sex, genotype)
    na.omit(merge(UNIQUE.INDIVIDUAL.INFO, cell.phenotypes))->unique.individual.cell.phenotypes
    #make recalled table: phenotype.day1 phenotype.day2
    recalled.phenotypes.table(only.keep.recalled.individuals(cell.phenotypes, recalled.individuals.pch=RECALLED.INDIVIDUALS.PCH))->recalled.individual.cell.phenotypes
    return(list(all=individual.cell.phenotypes, unique=unique.individual.cell.phenotypes, recalled=recalled.individual.cell.phenotypes))
}



#load.individual.cell.phenotypes("~nikolas/CellPhenotypes/calli.cell.phenotypes")->calli.cell.phenotypes
#load.individual.cell.phenotypes("~nikolas/CellPhenotypes/vincent.cell.phenotypes")->vincent.phenotypes
#load.individual.cell.phenotypes("~nikolas/CellPhenotypes/auto2.cell.phenotypes")->auto.phenotypes
#load.individual.cell.phenotypes("~nikolas/CellPhenotypes/sp.mm.cell.phenotypes")->sp.mm.cell.phenotypes
#load.individual.cell.phenotypes("~nikolas/CellPhenotypes/mm.cell.phenotypes")->mm.cell.phenotypes
#load.individual.cell.phenotypes("~nikolas/CellPhenotypes/manual.cell.phenotypes")->manual.cell.phenotypes
#load.individual.cell.phenotypes("~nikolas/CellPhenotypes/gates.2.cell.phenotypes")->gates.2.cell.phenotypes



