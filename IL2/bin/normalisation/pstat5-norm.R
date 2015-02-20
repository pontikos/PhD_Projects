
d <- read.csv('pstat5-peaks.csv', stringsAsFactors=F)

mapply(
       function(f, a, b) {
           X <- read.csv(f)
           X[,'PSTAT5'] <- b*X[,'PSTAT5']+a
           write.csv( X, file=gsub('.fcs.csv','_pstat5-norm.csv',f), quote=F, row.names=F )
       },
       file.path('Lymphocytes5/RData', sprintf('%s.csv', basename(d$fcsFile))),
       d$pSTAT5.loc.a,
       d$pSTAT5.loc.b
)

