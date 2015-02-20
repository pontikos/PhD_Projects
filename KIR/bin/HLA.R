setwd('~nikolas/stats-archive/Papers/KIR/Data') 

print(load('qPCR/MI-10-datasets.RData'))

# HLA data
head(hla.data <- read.table('Data/HLA/cc-hla-2013-06-21.tab',header=T))
dim( hla.data )
hla.data$t1d <- hla.data$t1d-1
hla.data$dil_subjectid <- hla.data$dil_subject
dim( hla.data <- subset(hla.data, t1d %in% 0:1) )
hla.data$HLAA <- with(hla.data,paste(HLAA_Bw4_aa80_1, HLAA_Bw4_aa80_2, sep='-'))
hla.data$HLAB <- with(hla.data,paste(HLAB_Bw4_Bw6_aa80_1, HLAB_Bw4_Bw6_aa80_2, sep='-'))
hla.data$HLA_Bw <- ifelse(grepl('6N', with(hla.data, paste(HLAA,HLAB,sep='-'))), '6N', '0')
hla.data$HLA_Bw[grepl('4T', with(hla.data, paste(HLAA,HLAB,sep='-')))] <- '4T'
hla.data$HLA_Bw[grepl('4I', with(hla.data, paste(HLAA,HLAB,sep='-')))] <- '4I'
#drop column collection
head(hla.data <- hla.data[,-which(colnames(hla.data)=='collection')])
# 20,445 = 9,174 : 11,271  


case.only.chisq.test <- function(b) {
    dim( cases <- subset(cases, !is.na(cases[,paste(b,'kir',sep='.')])) )
    # kir presence/absence
    table( kir3ds1 <- ifelse(as.numeric(gsub('(.)-.','\\1',cases[,paste(b,'kir3ds1',sep='.')])>0), 'KIR3DS1+', 'KIR3DS1-') )
    table( kir3dl1 <- ifelse(as.numeric(gsub('.-(.)','\\1',cases[,paste(b,'kir3dl1',sep='.')])>0), 'KIR3DL1+', 'KIR3DL1-') )
    table( kir <- paste( kir3ds1, kir3dl1, sep='/' ) )
    # hlabw4 presence/absence
    table( hlabw4 <- c('HLA-Bw4-','HLA-Bw4+')[1+as.numeric(grepl('4',cases$HLA_Bw))] )
    table( hlabw4I <- c('HLA-Bw4-80I-','HLA-Bw4-80I+')[1+as.numeric(grepl('4I',cases$HLA_Bw))] )
    # kir3dl1
    print(xtable(table(kir3dl1, hlabw4)))
    print(chisq.test(table(kir3dl1, hlabw4)))
    # kir3ds1
    print(xtable(table(kir3ds1, hlabw4I)))
    print(chisq.test(table(kir3ds1, hlabw4I)))
    # kir
    print(xtable(table(kir, hlabw4)))
    print(chisq.test(table(kir, hlabw4)))
    #xtable(table(kir, hlabw4I))
    #chisq.test(table(kir, hlabw4I))
} 


mi <- lapply(imputations, function(i) {
       i <- subset(i, t1d==1)
       i <- merge(i, hla.data)
       # kir presence/absence
       i$kir3ds1 <- ifelse(as.numeric(gsub('(.)-.','\\1',i$geno)>0), 'KIR3DS1+', 'KIR3DS1-')
       i$kir3dl1 <- ifelse(as.numeric(gsub('.-(.)','\\1',i$geno)>0), 'KIR3DL1+', 'KIR3DL1-')
       # hlabw4 presence/absence
       i$hlabw4 <- c('HLA-Bw4-','HLA-Bw4+')[1+as.numeric(grepl('4',i$HLA_Bw))] 
       i$hlabw4I <- c('HLA-Bw4-80I-','HLA-Bw4-80I+')[1+as.numeric(grepl('4I',i$HLA_Bw))] 
       return(i)
})


