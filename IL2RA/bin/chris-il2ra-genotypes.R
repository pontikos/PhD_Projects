library(snpStats)

### FIRST STEP: get genotypes of ichip.individuals in PLINK files

#Um, why not use the PLINK files developed for that purpose then?  I thought these were already in use in the stats group?
#Niko: these are in:
#/ipswich/data/Immunochip/PLINK/distribution/
#and you can either use plink itself to carve a subset of chromosome 10 out of the files of interest, or read them into R and proceed from there.
#There are 3 datasets in which you will find these IP samples:
#[a] cbr100-preqc.* - these are the IP samples (IDs beginning CB) done by David van Heel, and labelled "cbr100" if you grep for your list in
#/ipswich/data/Immunochip/support/cbr/cbr-immunochip-2012-06-07.tab 
#[b] cbrmicro-preqc.* - these are the ones labelled cbrmicro
#[c] cbr-preqc.* - these are the ones labelled vasculitis.
#Note [a] and [c] are parts of someone else's experiments, hence the odd naming.
# individuals in this study which have been ichipped

head(ichip.individuals <- read.table('~nikolas/Projects/IL2RA/individuals-ichip.txt'))
# all cbr indviduals which have been ichipped
#read.table('~nikolas/Projects/IL2RA/cbr-immunochip-2012-06-07.tab')
print(il2ra <- read.table('~nikolas/Projects/IL2RA/snps-for-niko-jan-2014.csv', header=T))
print(il2ra <- read.table('~nikolas/Projects/IL2RA/summw2.csv', header=T))
print(il2ra$snp <- gsub('\\.', '-', as.character(il2ra$snp)))
all.ids <- unique(as.character(unlist(ichip.individuals[,2:4])))
## READ PLINK FILES
basedir <- '/ipswich/data/Immunochip/PLINK/distribution'
#[a] cbr100-preqc.* - these are the IP samples (IDs beginning CB) done by David van Heel, and labelled "cbr100" if you grep for your list in
fam <- file.path(basedir, 'cbr100-preqc.fam')
bim <- file.path(basedir, 'cbr100-preqc.bim')
bed <- file.path(basedir, 'cbr100-preqc.bed')
cbr100 <- read.plink(bed, bim, fam)
dim(cbr100.genotypes <- cbr100$genotypes[na.omit(pmatch(all.ids, rownames(cbr100$genotypes))),intersect(il2ra$snp,colnames(cbr100$genotypes))])
#[b] cbrmicro-preqc.* - these are the ones labelled cbrmicro
fam <- file.path(basedir, 'cbrmicro-preqc.fam')
bim <- file.path(basedir, 'cbrmicro-preqc.bim')
bed <- file.path(basedir, 'cbrmicro-preqc.bed')
cbrmicro <- read.plink(bed, bim, fam)
dim(cbrmicro.genotypes <- cbrmicro$genotypes[na.omit(pmatch(all.ids, rownames(cbrmicro$genotypes))),intersect(il2ra$snp,colnames(cbrmicro$genotypes))])
#[c] cbr-preqc.* - these are the ones labelled vasculitis.
fam <- file.path(basedir, 'cbr-preqc.fam')
bim <- file.path(basedir, 'cbr-preqc.bim')
bed <- file.path(basedir, 'cbr-preqc.bed')
cbr <- read.plink(bed, bim, fam)
dim(cbr.genotypes <- cbr$genotypes[na.omit(pmatch(all.ids, rownames(cbr$genotypes))),intersect(il2ra$snp,colnames(cbr$genotypes))])
# rbind [a] [b] [c]
cbr.il2ra.geno <- as.data.frame(rbind( cbr.genotypes, cbrmicro.genotypes, cbr100.genotypes ))
# there are individuals which were analysed in both VanHeel and (UVA xor Sanger)
cbr.il2ra.geno$individual <- ichip.individuals[pmatch(rownames(cbr.il2ra.geno), ichip.individuals[,2]),4]
dim(cbr.il2ra.geno)
dim(cbr.il2ra.geno <- cbr.il2ra.geno[!duplicated(cbr.il2ra.geno),])
# SNPs we have found
print(snps <- intersect(il2ra$snp,colnames(cbr.il2ra.geno)))
# convert genotypes to numeric
cbr.il2ra.geno[,snps] <- apply(cbr.il2ra.geno[,snps],2,as.numeric)
write.csv(cbr.il2ra.geno, file='~nikolas/Projects/IL2RA/cbr-il2ra-geno.csv', row.names=row.names(cbr.il2ra.geno))


### SECOND STEP: CD25-related phenotype to IL2RA genotype correlation
dim(cbr.il2ra.geno <- read.csv('~nikolas/Projects/IL2RA/cbr-il2ra-geno.csv', row.names=1))
print(il2ra <- read.table('~nikolas/Projects/IL2RA/snps-for-niko-jan-2014.csv', header=T))
# SNPs we have found
print(snps <- intersect(il2ra$snp,colnames(cbr.il2ra.geno)))
# missing SNPs
setdiff(il2ra$snp,colnames(cbr.il2ra.geno))

# disagreement of genotype between UVA and vanHeel for SNP imm_10_6148346
# turns out the MAF is the other way round since this is a 50:50 SNP
apply(cbr.il2ra.geno[grep('CB00310L', cbr.il2ra.geno$individual),],2,unique)

# Calliope Dendrou's cell phenotypes
calli<-read.csv('~nikolas/Projects/IL2RA/Calli_CD25bright_CBR200.csv')
cat('number of individuals overlap\n')
length(unique(intersect(cbr.il2ra.geno$individual, calli$individual)))
cat('number of individuals which are not in Calli_CD25bright_CBR200.csv\n')
length(unique(setdiff(cbr.il2ra.geno$individual, calli$individual)))
dim(d<-merge(cbr.il2ra.geno, calli, by='individual'))
# APC is the fluorochrome used for CD25
print(cd25.phenotypes <- c(grep('MeanAPC_norm',colnames(d),value=T),grep('CD25pos.*freqPar',colnames(d),value=T)))
#ignore cd25 phenotypes which have NA
print(cd25.phenotypes <- cd25.phenotypes[apply(d[,cd25.phenotypes],2,function(x) sum(is.na(x))==0)])


# this is to convert 1, 2, 3 into TT, TC, CC
num.geno <- t(apply(snp.map[gsub('\\.','-',snps),c('allele.1','allele.2')], 1, function(x) c(paste(x[1],x[1],sep=''),paste(x[1],x[2],sep=''),paste(x[2],x[2],sep=''))))
dim(d[,snps] <- sapply(snps, function(snp) num.geno[gsub('\\.','-',snp),d[,snp]]))
write.csv(d[,c('individual', 'date', 'age', 'sex', snps, cd25.phenotypes)], file='~nikolas/Projects/IL2RA/cd25phen-gen.csv', quote=F, row.names=F)


### Association with SNPs
# snp map from PLINK (bim file)
snp.map <- cbr100$map
assoc <- function(cd25.pheno, snp) {
        pheno <- d[,cd25.pheno]
        snp <- gsub('\\.','-',snp)
        s <- d[,snp]
        m.params <- summary(lm( pheno ~ s ))$coefficients[2,]
        names(m.params) <- c('beta', 'se.beta', 't.value', 'p.value')
        data.frame(
        snp=snp,
        position=snp.map[snp,'position'],
        allele1=snp.map[snp,'allele.1'],
        allele2=snp.map[snp,'allele.2'],
        phenotype=cd25.pheno,
        p.value=m.params[['p.value']],
        beta=m.params[['beta']],
        se.beta=m.params[['se.beta']]) 
}
# for all phenotype x SNP combinations (15 phenotypes x 29 snps)
phen.gen.assoc <- do.call('rbind', apply(expand.grid(phenotype=as.character(cd25.phenotypes), snp=as.character(snps), stringsAsFactors=F), 1, function(x) assoc(x[1], x[2])))
write.csv(phen.gen.assoc, file='~nikolas/Projects/IL2RA/cd25phen-gen-assoc.csv', row.names=F)

dim(phen.gen.assoc <- read.csv('~nikolas/Projects/IL2RA/cd25phen-gen-assoc.csv'))
phen.gen.assoc$phenotype <- gsub('Lymphocytes.CD4.(.*).MeanAPC_norm', '\\1', phen.gen.assoc$phenotype)
# Manhattan plot showing SNP assoc with the different phenotypes
require(ggplot2)
ggplot(phen.gen.assoc, aes(x=position, y=-log10(p.value), labels=snp, group=phenotype, col=phenotype)) + geom_point() + ylab('-log10 p-value') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_x_continuous(breaks=unique(phen.gen.assoc$position), labels=unique(as.character(phen.gen.assoc$snp))) + geom_line(size=.25)
#+ theme(legend.position="none") 


# correlation between CD25 MFI phenotypes
pairs(d[,cd25.phenotypes])
cor(d[,cd25.phenotypes])

# correlation plot of phenotypes
require(GGally)
X<-d[,cd25.phenotypes]
colnames(X) <- gsub('Lymphocytes.CD4.(.*).MeanAPC_norm', '\\1', colnames(X))
ggcorr(X, size=2)




