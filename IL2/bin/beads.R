
beads <- read.table('beads-tony-2014-01-07.csv', sep=';', header=T)

source('~nikolas/bin/FCS/fcs.R')

basedir <- '~nikolas/facs'
b.single <- read.FCS(file.path(basedir,'/Tony/Tony/19_June_2102_pSTAT5_TT_NKCD4CD8_singles/beads_pSTAT5_AF488.fcs'),'PSTAT5')
b.comp <- read.FCS(file.path(basedir,'/Tony/Tony/03_sept_2102_compbeads/beads_pSTAT5_AF488.fcs'),'PSTAT5')
b.cytocal <- read.FCS(file.path(basedir,'/Tony/Tony/02_August_2102_pSTAT5_D-GAP/cytocal_020812_beads.fcs'), 'PSTAT5')

b <- read.FCS(file.path(basedir,'/Tony/Tony/130213_NSG_Pro_tissues_Treg/Beads_pSTAT5.fcs'),'HCD45')

par(mfrow=c(3,1))
plot(density(getChannels(b.single, 'PSTAT5'),bw=.1),main='single',xlab='pstat5')
plot(density(getChannels(b.comp, 'PSTAT5'),bw=.1),main='comp',xlab='pstat5')
plot(density(getChannels(b.cytocal, 'PSTAT5'),bw=.1),main='cytocal',xlab='pstat5')

plot(density(getChannels(b, 'HCD45'),bw=.1),main='single',xlab='hcd45')

