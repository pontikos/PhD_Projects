options(stringsAsFactors = FALSE)

x <- read.csv('~nikolas/Projects/IL2RA/Calli_CD25bright_CBR200.csv')
x$fcsFile <- gsub('.fcs','', tolower(x$fcsFile))

v <- read.csv('~nikolas/Projects/IL2RA/CellPhenotypes/vincent.cell.phenotypes')
v$fcsFile <- tolower(v$fcsFile)

print(dim(x <- merge(x,v,by=c('fcsFile','individual'))))
#CB00503W is not included
print(dim(x <- x[which(x$individual != 'CB00503W'),]))

x$rs12722495 <- x$DIL9620.CD25.rs12722495
x$rs2104286 <- x$DIL8103.CD25.rs2104286
x$rs11594656 <- x$DIL10847.CD25.rs11594656

#individual CB00556D has no genotype information
x[which(is.na(x$rs12722495)),'individual']
x <- x[-which(x$individual=='CB00556D'),]




cbind(x[,c('rs12722495','rs2104286','rs11594656')])

table(x[,c('rs12722495','rs2104286','rs11594656')])

#rs12722495
table(x[,c('rs12722495','sex')])
#tapply(x$age, list(x$rs12722495,x$sex), median)
tapply(x$age, x$rs12722495, mean)

#rs2104286
table(x[,c('rs2104286','sex')])
tapply(x$age, x$rs2104286, mean)

#rs11594656
cbind(table(x[,c('rs11594656','sex')]),tapply(x$age, x$rs11594656, mean))



cbind(x[,c('rs12722495','rs2104286','rs11594656')])
table(paste(x$rs12722495, x$rs2104286, x$rs11594656, sep='-'))

unphased.genotype <- paste(x$rs12722495, x$rs2104286, x$rs11594656, sep='-')
write.csv( cbind( table(unphased.genotype,x$sex) , mean.age=round(tapply(x$age, unphased.genotype, mean),1) ), quote=FALSE )

