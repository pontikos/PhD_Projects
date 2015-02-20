suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowCore, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
#suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowViz, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
library("geneplotter")

as.character(read.csv('~nikolas/bad.APC.mfi.fcs')[,1]) -> fcs.file.names
#as.character(read.csv('~nikolas/good.APC.mfi.fcs')[,1]) -> fcs.file.names

wsp.obj <- function(x)
    sprintf('~nikolas/dunwich/WSP.R.obj/%s_%s_Treg.wsp.Rdata', toupper(x[[1]]), format(as.Date(x[[2]], format='%Y%b%d'), '%d%m%y'))

#wsp objects
sapply( strsplit(fcs.file.names, '_'), wsp.obj ) -> wsp.files

#fcs files
sapply(fcs.file.names, function(x) sprintf('~nikolas/dunwich/FCS.Gated/Calli/NON_T_REGS/%s.fcs', x)) -> fcs.files


data.frame(base.name=fcs.file.names, fcs=as.character(fcs.files), wsp=as.character(wsp.files))->files

#read.csv('~nikolas/dunwich/cd45ra.thresholds.stats') -> thresholds
#thresholds$base.name <- gsub('.fcs', '', thresholds$fcsFile)
#merge(thresholds, files, by=c('base.name'))->files

#pdf('good.gates.pdf')

par(mfrow=c(1,1), mar=c(2, 2, .75, .75))
pdf('~/Plots/Learning/bad.gates.pdf')
for (r in 1:nrow(files)) {

    files[r,] -> d
    as.character(files[r,'wsp'])->wsp
    as.character(files[r,'fcs'])->fcs
    as.character(files[r,'base.name'])->base.name

    cat(base.name, '\n')

    mm.env <- new.env()
    load(sprintf('~nikolas/dunwich/FCS.Gated/Calli/NON_T_REGS/mm.%s.fcs.obj', base.name), env=mm.env)
    load(sprintf('~nikolas/dunwich/FCS.Gated/Calli/NON_T_REGS/np.mm.%s.fcs.obj', base.name), env=mm.env)
    mm.env$mm -> mm
    mm.env$np.mm  -> np.mm

    cbind(mm$x, mm$posterior) -> mm$xy
    mm$xy[order(mm$xy[,1]),]->mm$xy

    cbind(np.mm$data, np.mm$posteriors) -> np.mm$xy
    np.mm$xy[order(np.mm$xy[,1]),]->np.mm$xy

    #load flowlist
    print(wsp)
    load(wsp)
    #print(flowlist@FlowJoWorkspace)
    flowlist@flowJoData[[1]]->f1
    sprintf('%s.fcs:lymphocytes:cd4:cd127hi:127hi:cd45raneg$', base.name)->p
    #print(p)
    f1@filter[grep(p, f1@filterName, ignore.case=T)][[1]]@filters[[2]]->cd45raneg
    print(cd45raneg@min)
    print(cd45raneg@max)
    sprintf('%s.fcs:lymphocytes:cd4:cd127hi:127hi:cd45rapos$', base.name)->p
    f1@filter[grep(p, f1@filterName, ignore.case=T)][[1]]@filters[[2]]->cd45rapos
    print(cd45rapos@min)
    print(cd45rapos@max)

    sprintf('%s.fcs:lymphocytes:cd4:cd127hicd25pos:cd45rapos$', base.name)->p
    sprintf('%s.fcs:lymphocytes:cd4:cd127hicd25pos$', base.name)->p
    min((f1@filter[grep(p, f1@filterName, ignore.case=T)][[1]]@filters[[2]])@boundaries[,1])->cd25pos
    

    read.FCS(fcs)->fcs.data

    #~nikolas/dunwich/Beads.Stats/beads.auto.fcs3.stats
    read.csv('~nikolas/Beads.Stats/beads.manual.fcs2.stats') -> beads.stats
    #unlist(strsplit(toupper(beads.stats$file.name),".FCS"))->beads.stats$file.name
    toupper(beads.stats$file.name)->beads.stats$file.name
    as.Date(beads.stats$date) -> beads.stats$date

    as.Date(fcs.data@description[["$DATE"]],format="%d-%b-%Y")->experiment.date
    fcs.data@description[["EXPERIMENT NAME"]]->experiment.name

    #find the bead data based on experiment name
    grep(experiment.name, ignore.case=T, beads.stats$file.name)->i

    blank.beads.mfi <- beads.stats[i,'mfi.1']
    dim.beads.mfi <- beads.stats[i,'mfi.2']

    #pacific cd45ra
    fcs.data@exprs[,8]->x
    length(x) -> n
    round(100*length(x[x<cd45raneg@max])/n,2)->naive.prop
    round(100*length(x[x>cd45rapos@min])/n,2)->memory.prop
    round(100*length(x[cd45raneg@max < x & x < cd45rapos@min])/n,2)->neither.prop
    log10(x[x>1]) -> y

    plot(density(y, bw=.1), ylim=c(0,1), xlim=c(0,2), main=base.name)
    #abline(v=log10(cd45raneg@max), col="green")
    #text(log10(cd45raneg@max), 0.4, naive.prop)
    #abline(v=log10(cd45rapos@min), col="green")
    #text(log10(cd45rapos@min), 0.2, memory.prop)
    #text((log10(cd45raneg@max)+log10(cd45rapos@min))/2, 0, neither.prop)
    #abline(v=files[r, 'medoids.mean'], col="black", lty=2)
    #abline(v=files[r, 'mean'], col="red")
    #abline(v=files[r, 'valley'], col="blue")
    read.csv('~/CellPhenotypes/Gates/mm.csv')->mm.gates
    abline(v=mm.gates[base.name==mm.gates$fcsFile, 'cd45ra.neg'], col="blue")
    abline(v=mm.gates[base.name==mm.gates$fcsFile, 'cd45ra.pos'], col="blue")
    #mm.gates[base.name==mm.gates$fcsFile, 'cd45ra.pos']
    #abline(v=files[r, 'mm.1'], col="blue")
    #abline(v=files[r, 'mm.2'], col="blue")
    lines(sort(y), mm$lambda[[1]]*dnorm(sort(y), mm$mu[[1]], mm$sigma[[1]])+mm$lambda[[2]]*dnorm(sort(y), mm$mu[[2]], mm$sigma[[2]]), lty=2)
    #
    #lines(rep(d$medoids.1, 2), c(0, dnorm(d$medoids.1, mm$mu, mm$sigma)), lty=2)
   #lines(rep(d$medoids.2, 2), c(0, dnorm(d$medoids.2, mm$mu, mm$sigma)), lty=2)
    #
    #lines(np.mm$xy[,c(1,2)], col='blue')
    #lines(np.mm$xy[,c(1,3)], col='blue')
    #
    #lines(mm$xy[,1], 1-mm$xy[,2], col='red')
    #lines(mm$xy[,1], dy, col='green')
#plot(np.mm, lty=2, newplot=FALSE, addlegend=TRUE)
#plot(ecdf(np.mm$data))
#plot(ecdf(np.mm$posteriors[,2]))
#cdf.mm<-compCDF(cbind(mm$xy[,1]), mm$xy[,2:3], lwd=2)
#lines(mm$xy[,1], mm$xy[,3], col='red')
    #abline(v=1.09, col="red")
    break
}
dev.off()

# 2d plot: cd25 and cd45ra
def.par <- par(no.readonly = TRUE)
par("mar" = c(1,1,1,1))
layout(matrix(c(2,1,1,2,1,1,4,3,3), 3, 3, byrow = TRUE))

#CD127
#log10(fcs.data@exprs[,c(3)])->x
#APC
log10(fcs.data@exprs[,c(5)])->x
#CD45RA
log10(fcs.data@exprs[,c(8)])->y
x>0 & y>0->i
x[i]->x
y[i]->y
#density(x, bw=.01) -> d.x
#density(y, bw=.01) -> d.y
plot(x, y, pch=".", col="gray")
points(mean(x), mean(y), pch="x")
mixtools::ellipse(mu=c(mean(x),mean(y)), sigma=cov(cbind(x,y)), alpha=.5)
#cd45ra
abline(h=log10(cd45raneg@max), col="green")
abline(h=log10(cd45rapos@min), col="green")
abline(h=files[r, 'valley'], col="blue")
abline(h=d$medoids.1, col="red")
abline(h=d$medoids.2, col="red")
#cd25
abline(v=log10(cd25pos), col="green")
abline(v=log10(blank.beads.mfi), col="blue")
#abline(v=log10(dim.beads.mfi), col="blue")
#text(log10(cd45raneg@max), 0.4, naive.prop)
breaks <- 100
yh <- hist(y, breaks = breaks, plot = FALSE)
#ignore first count
#yh$counts[[1]]<-0
bar.colours <- rep('gray', length(yh$counts))
bar.colours[[which.min(abs(yh$mids-log10(cd45raneg@max)))]] <- 'green'
bar.colours[[which.min(abs(yh$mids-log10(cd45rapos@min)))]] <- 'green'
#bar.colours[[which.min(abs(yh$mids-files[r, 'valley']))]] <- 'blue'
bar.colours[[which.min(abs(yh$mids-d$medoids.1))]] <- 'red'
bar.colours[[which.min(abs(yh$mids-d$medoids.2))]] <- 'red'
barplot(-(yh$counts),space=0,horiz=T, axes = FALSE, col=bar.colours)
xh <- hist(x, breaks = breaks, plot = FALSE)
#ignore first count
#xh$counts[[1]]<-0
bar.colours <- rep('gray', length(xh$counts))
bar.colours[[which.min(abs(xh$mids-log10(cd25pos)))]] <- 'green'
bar.colours[[which.min(abs(xh$mids-log10(blank.beads.mfi)))]] <- 'blue'
barplot(-(xh$counts),space=0,horiz=F, axes = FALSE, col=bar.colours)

#plot(d.x$x, d.x$y, xlim=range(x), type='l')
#plot(d.y$y, d.y$x, ylim=range(y), xlim=rev(range(d.y$y)), type='l')

par(def.par)

#plot without log10 transform
    #plot((fcs.data@exprs[,c(5,8)]), pch=".")
    #abline(h=(cd45raneg@max), col="green")
    #abline(h=(cd45rapos@min), col="green")
    #abline(v=cd25pos, col="green")
    #abline(h=10**(files[r, 'valley']), col="blue")
    #abline(v=blank.beads.mfi, col="blue")
}

dev.off()




