source('~nikolas/bin/FCS/fcs.R')

head(d <- read.table('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/all-FCS-panel-date-ids-dose.csv', sep=';', header=T, stringsAsFactors=FALSE)) 
channels <- c('FSCW','SSCW','FSCH','SSCH','FSCA','SSCA','CD4','CD25','CD45RA','FOXP3','PSTAT5')

out.dir <- '/chiswick/data/store/facs/Tony-RData/PSTAT5/CD25/CD45RA/CD4/FOXP3'

for (i in 1:nrow(d)) {
    fcs <- d$fcsFile[i]
    dose <- d$dose[i]
    date <- d$date[i]
    print(id <- grep('^[C|K]', unlist(strsplit(d$ids[i],',')), value=TRUE)[[1]])
    print(fcs)
    print(out.file <- file.path(out.dir, sprintf('%s.RData',paste(id, dose, date, sep='_'))))
    if (!file.exists(out.file)) {
        print(dim(fcs.data <- getChannels(read.FCS(fcs,channels=channels,TRANS=NULL),channels)))
        save(fcs.data,file=out.file)
    } else {
        print('FILE EXISTS!')
        print(load(out.file))
        print(dim(fcs.data))
    }
}


channels <- c('FSCW','SSCW','FSCH','SSCH','FSCA','SSCA', 'CD25','CD3','CD4','CD45RA','CD56','CD8','FOXP3','PSTAT5')
for (fcs in list.files('~/dunwich/Projects/IL2/Tony-FCS/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/',pattern='.*.fcs',full.names=TRUE)) {
    print(dim(fcs.data <- getChannels(read.FCS(fcs,channels=channels,TRANS=NULL),channels)))
    out.file <- file.path('~/dunwich/Projects/IL2/Tony-RData/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/',gsub('.fcs','.RData',basename(fcs)))
    save(fcs.data,file=out.file)
}


