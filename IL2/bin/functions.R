library(ggplot2)
library(gtable)
library(scales)
library(reshape2)
library(flowBeads)
library(spade)
library(fastcluster)


dens.fun <- function(x,...) {
    d <- density(x,...)
    return(splinefun(d$x,d$y))
}





load.FCS <- function(ext='.fcs', subdir='', doses=c('0U', '1U', '10U', '1000U'), verbose=T) {
    rep.fcs <- list()
    for (individual in rep.individuals) {
        rep.fcs[[individual]] <- list()
        for (day in c('day1', 'day2')) {
            rep.fcs[[individual]][[day]] <- list()
            base.dir <- file.path(dir, individual, day, subdir)
            for (dose in doses) {
                file.name <- file.path(base.dir, paste(dose, ext, sep=''))
                if (verbose) print(file.name)
                rep.fcs[[individual]][[day]][[dose]] <- read.FCS(file.name)
            }
        }
    }
    rep.fcs
}


update.FCS <- function(ext='.fcs', subdir='', doses=c('0U', '1U', '10U', '1000U'), verbose=T) {
    for (individual in rep.individuals) {
        for (day in c('day1', 'day2')) {
            base.dir <- file.path(dir, individual, day, subdir)
            for (dose in doses) {
                file.name <- file.path(base.dir, paste(dose, ext, sep=''))
                if (verbose) print(file.name)
                update.FCS.description(file.name, keywords=list('INDIVIDUAL'=individual, 'DOSE'=dose))
            }
        }
    }
}


multi.density <- function(fcs.data, channels=c('FSC-A', 'SSC-A', 'CD25', 'CD4$', 'CD45RA'), kernel_mult=5.0, apprx_mult=1.5, med_samples=2000) {
    SPADE.density(getChannels(fcs.data, channels), kernel_mult, apprx_mult, med_samples)
}


downsample.indexes <- function(fcs.data, density, desired_samples=20000, verbose=TRUE) {
    fcs.data <- fcs.data@exprs
    fcs.data <- cbind(fcs.data, density)
    # exclude_pctile
    boundary1 <- quantile(density,0.01,names=FALSE)
    out_data <- subset(fcs.data, density > boundary1) # Exclusion    
    if (desired_samples >= nrow(out_data)) return( which(density > boundary1) )
    # Need to find target density such there are approximately desired_samples
    # remaining after downsampling. To do so we solve for the density such that
    # the sum of samples below that density plus the expected value of
    # samples retained above that density equals approximately the desired
    # number of samples
    density_s <- sort( out_data[,'density'] )
    cdf       <- rev(cumsum(1.0/rev(density_s)))
    # Default solution if target density smaller than any present
    boundary2 <- desired_samples/cdf[1] 
    # Boundary actually falls amongst densities present
    if (boundary2 > density_s[1]) {
      targets <- (desired_samples-1:length(density_s)) / cdf 
      boundary2 <- targets[which.min(targets-density_s > 0)]
    }
    if (verbose) print(c(boundary1, boundary2))
    return( which( (boundary1 < density) & (boundary2/density > runif(length(density))) ) )
}


downsample.FCS <- function(fcs.data, density, desired_samples=20000, verbose=FALSE) {
    indexes <- downsample.indexes(fcs.data, density, desired_samples, verbose)
    return(flowFrame(fcs.data@exprs[indexes,],parameters(fcs.data),description=description(fcs.data)))
}



# Cluster observations into ~k clusters
SPADE.cluster <- function(fcs.data, k=200, channels=c('FSC-A', 'SSC-A', 'CD25', 'CD4$', 'CD45RA')) {
    i <- sapply(channels, function(x) grep( x, parameters(fcs.data)@data$desc, ignore.case=T ))
    fcs.data <- fcs.data@exprs[,i]
	if (nrow(fcs.data) > 60000) warning("Potentially too many observations for the clustering step",immediate=TRUE);
  # Transpose table before call into row major order
	clust <- .Call("SPADE_cluster",t(fcs.data),as.integer(k))
	# Invalid clusters have assgn == 0
	centers = c()
    is.na(clust$assgn) <- which(clust$assgn == 0)
	for (i in c(1:max(clust$assgn, na.rm=TRUE))) {  
		obs <- which(clust$assgn == i)
		if (length(obs) > 1) {
			centers <- rbind(centers,colMeans(fcs.data[obs,,drop=FALSE]))
			clust$assgn[obs] <- nrow(centers)
		} else {
			is.na(clust$assgn) <- obs
		}
    }
    return(list(centers=centers,assign=clust$assgn))
}



flatten <- function(rep.fcs) {
    flat.fcs <- unlist(rep.fcs)
    names(flat.fcs) <- sapply(flat.fcs, function(x) paste(getDate(x), getIndividual(x), getDose(x), sep='.'))
    flat.fcs
}


extractDate <- function(x) unique(simplify2array(strsplit(names(x), '\\.'))[1,])
extractIndividual <- function(x) unique(simplify2array(strsplit(names(x), '\\.'))[2,])
extractDose <- function(x) unique(simplify2array(strsplit(names(x), '\\.'))[3,])

unflatten.day <- function(flat, doses, f=identity) {
    day1 <- sapply(doses, function(dose) flat[grep(sprintf('\\.%s',dose), names(flat))][c(TRUE,FALSE)], simplify='array', USE.NAMES=FALSE)
    day2 <- sapply(doses, function(dose) flat[grep(sprintf('\\.%s',dose), names(flat))][c(FALSE,TRUE)], simplify='array', USE.NAMES=FALSE)
    return(list(day1=f(day1), day2=f(day2)))
}

unflatten.dose <- function(flat, resting='0U', stimulated='1000U', do.unlist=TRUE) {
    resting.dose <- flat[grep(sprintf('\\.%s',resting), names(flat))]
    stimulated.dose <- flat[grep(sprintf('\\.%s',stimulated), names(flat))]
    if (do.unlist)
        return(list(resting=unlist(resting.dose), stimulated=unlist(stimulated.dose)))
    else
        return(list(resting=(resting.dose), stimulated=(stimulated.dose)))
}




# upsampling
SPADE.assignToCluster <- function(tbl, cluster_data, cluster_assign)
    .Call("SPADE_assign",t(tbl),t(cluster_data),as.integer(cluster_assign))




### beads
plot.beads <- function(beads.fcs, dir='~nikolas/IL2/Plots/beads/') {
    dir.create(dir, recursive=T, showWarnings=F)
    for (p in getParams(beads.fcs[[1]])) {
        pdf(file.path(dir, sprintf('%s.pdf',p)))
        for (b in beads.fcs) {
            plot(b, p)
        }
        dev.off()
    }
}

channel.timeseries <- function(flow.data, beads.fcs, channels=c('CD25', 'CD4', 'pSTAT5', 'CD45RA'), dir='~nikolas/IL2/Plots/time-series/') {
    dir.create(dir, recursive=T, showWarnings=F)
    panel <- data.frame(bead.channel=gsub('-', '.', gsub('-A$', '', gsub(' ', '.', toupper(parameters(flow.data[[1]])@data[,'name'])))), channel=parameters(flow.data[[1]])@data[,'desc'])
    all.d <- data.frame()
    for (channel in channels) {
        bead.channel <- panel[which(channel==panel$channel),'bead.channel']
        print(bead.channel)
        lymph.d <- data.frame(date=as.Date(sapply(flow.data, getDate)), dose=factor(sapply(flow.data, getDose), levels=c('0U','1U','10U','1000U'), ordered=T), mfi=as.numeric(sapply(flow.data, function(x) getMFI(x,channel))), row.names=NULL)
        beads.d <- data.frame(date=as.Date(sapply(beads.fcs, getDate)), t(as.matrix(sapply(beads.fcs, function(x) x@clustering.stats['mean.fi',bead.channel,]))) , row.names=NULL, stringsAsFactors=FALSE)
        colnames(beads.d) <- gsub('X', 'beads', colnames(beads.d))
        beads.d <- melt(beads.d,'date')
        beads.d$value <- lgcl(beads.d$value)
        colnames(lymph.d) <- colnames(beads.d)
        d <- rbind(lymph.d, beads.d)
        d <- subset(d, variable %in% c('beads1', 'beads2','beads3', '0U', '1000U'))
        d$channel <- channel
        all.d <- rbind(all.d, d)
        #print(all.d)
    }
    g <- ggplot(all.d, aes(x=date, y=value, colour=variable, linetype=as.factor(grepl('beads',variable)), group=variable)) + geom_line() + xlab('') + ylab('MFI') + theme_bw() 
    #+ scale_x_date(labels=date_format('%b %y')) 
    #grid
    g <- g + facet_wrap( ~ channel, ncol=2)
    #horizontal
    #g <- g + facet_grid(. ~ channel)
    #vertical
    #g <- g + facet_grid(channel ~ .)
    g <- g + scale_colour_manual(values=c('black', 'green', 'blue', 'blue', 'blue'), name="", breaks=c('0U', '1000U', "beads1", 'beads2', 'beads3'), labels=c("resting", "1000U", 'beads1', 'beads2', 'beads3'))
    g <- g+ scale_linetype_discrete(guide=FALSE)
    #g <- g + scale_colour_discrete(name="", breaks=c('0U', '1000U', "beads1", 'beads2', 'beads3'), labels=c("resting", "1000U", 'beads1', 'beads2', 'beads3'))
    pdf(file.path(dir, 'all.pdf'), width=10, height=5)
    print(g)
    dev.off()
}

### repeatability per channel at two different days or same day different doses
### default is for different days
channel.repeatability <- function(rep.fcs, dir='~nikolas/IL2/Plots/day-effect-mean', day=c('day1', 'day2'), dose=rep('0U',2),
                                  channels=c('FSC-A', 'SSC-A', 'CD25', 'CD4$', 'pSTAT5', 'CD45RA'), fun=mean, verbose=TRUE, normalise=NULL, span=10, ...) {
    dir.create(dir, recursive=T, showWarnings=F)
    #R.squared <- function(x, y) return( 1-(mean((x-y)**2))/var(c(x,y)) )
    R.squared <- function(x, y) cor(x, y)**2
    d <- data.frame()
    class.f <- class(fun(1:10))[[1]]
    print(class.f)
    cor.channels <- data.frame(r.squared=numeric(length(channels)), row.names=channels, stringsAsFactors=FALSE)
    for (channel in channels) {
        if (verbose) print(channel)
        if (class.f %in% c('density','ecdf')) {
            pdf(file.path(dir, sprintf('%s.pdf',channel)), height=5, width=10) 
            par(mfrow=c(2,5), mar=c(1,1,1,1), oma=c(1,1,1,1))
        } else {
            pdf(file.path(dir, sprintf('%s.pdf',channel)))
            par(mfrow=c(1,1))
            x <- c()
            y <- c()
        }
        m.cor <- data.frame(row.names=rep.individuals, d1=numeric(length(rep.individuals)), d2=numeric(length(rep.individuals)), stringsAsFactors=FALSE)
        for (individual in rep.individuals) {
            if (verbose) print(individual)
            i <- grep( channel, parameters(rep.fcs[[individual]][[day[[1]]]][[dose[[1]]]])@data$desc, ignore.case=T )
            x1 <- rep.fcs[[individual]][[day[[1]]]][[dose[[1]]]]@exprs[,i]
            x2 <- rep.fcs[[individual]][[day[[2]]]][[dose[[2]]]]@exprs[,i]
            if (!is.null(normalise)) x2 <- normalise(x1, x2, ...)
            d1 <- fun(x1)
            d2 <- fun(x2)
            m.cor[individual,] <- cbind(median(x1), median(x2))
            if (class.f %in% c('density')) {
              x <- c(d1$x, d2$x)
              y <- c(d1$y, d2$y)
              #try and identify peaks
              plot(d1, xlab='', ylab='', main='', xaxt='n', yaxt='n', xlim=range(x), ylim=range(y), type='l')
              lines(d2, col='red')
              legend('topright', pch[pch$individual==individual,'pch'], cex=2, bty='n')
              abline(v=extract.landmarks(x1)$lms, col='black')
              abline(v=median(x1), lwd=2, lty=2, col='black')
              abline(v=extract.landmarks(x2)$lms, col='red')
              abline(v=median(x2), lwd=2, lty=2, col='red')
            } else if (class.f %in% c('ecdf')) {
              plot(d1, xlab='', ylab='', main='', xaxt='n', yaxt='n')
              lines(d2, col='red')
              legend('topright', pch[pch$individual==individual,'pch'], cex=2, bty='n')
            } else {
              x <- c(x, d1)
              y <- c(y, d2)
            }
        }
        cor.channels[channel,'r.squared'] <- R.squared(m.cor[,1],m.cor[,2])
        if (!(class.f %in% c('density','ecdf'))) {
            r <- range(c(x, y))
            R2 <- round(R.squared(x, y), digits=3)
            cat(channel, round(R2, 3), '\n')
            plot(x, y, xlab=paste(day[[1]], dose[[1]]), ylab=paste(day[[2]], dose[[2]]), xlim=r, ylim=r, pch=pch[,'pch'], cex=2, main=as.expression(bquote( r^2 == .(R2))), cex.main=2)
            #legend('topright', as.expression(bquote( r^2 == .(R2))), cex=2, bty='n')
            abline(b=1,a=0)
        }
        dev.off()
    } 
    return(cor.channels)
}




### boxplot distribution by channel across all people and order by date
channel.boxplot <- function(flat.fcs, dir='~nikolas/IL2/Plots/boxplot/', dose='.*') {
    dir.create(dir, showWarnings=F, recursive=T)
    flat.fcs <- flat.fcs[sort(names(flat.fcs))]
    #names(flat.fcs) <- t(as.data.frame(strsplit(names(flat.fcs), '\\.')))[,1]
    flat.fcs <- flat.fcs[grepl(dose, t(as.data.frame(strsplit(names(flat.fcs), '\\.')))[,3])]
    for (channel in  c('FSC-A', 'SSC-A', 'CD25', 'CD4$', 'pSTAT5', 'CD45RA')) {
    #for (channel in  c('FSC-A')) {
        pdf(file.path(dir, sprintf('%s.pdf',channel)), height=5, width=10) 
        par(mar=c(7,5,1,1))
        #par(mar=c(1,1,1,1), oma=c(1,1,1,1))
        boxplot(sapply(flat.fcs, function(x) getChannel(x, channel)), col=sapply(flat.fcs, function(x) pch[which(pch$name==getIndividual(x)),'col']), outline=F, las=2, main=channel,
                cex.axis=.5,
                cex.names=.5,
                cex.lab=.5,
                names=names(flat.fcs))
        #sapply(rep.fcs, getIndividual)
        #sapply(rep.fcs, getDose)
        #ggplot(data, aes(x=factor(Purpose), y=Rate)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5))
        dev.off()
    }
}


###
channel.downsample <- function(rep.fcs, rep.down.fcs, dir='~nikolas/IL2/Plots/downsampling-density', channels=c('FSC-A', 'SSC-A', 'CD25', 'CD4$', 'pSTAT5', 'CD45RA'), fun=mean, verbose=TRUE, dose='0U', day='day1') {
    dir.create(dir, recursive=T, showWarnings=F)
    #R.squared <- function(x, y) return( 1-(mean((x-y)**2))/var(c(x,y)) )
    R.squared <- function(x, y) cor(x, y)**2
    d <- data.frame()
    class.f <- class(fun(1:10))[[1]]
    print(class.f)
    for (channel in channels) {
        if (verbose) print(channel)
        if (class.f %in% c('density','ecdf')) {
            pdf(file.path(dir, sprintf('%s.pdf',channel)), height=5, width=10) 
            par(mfrow=c(2,5), mar=c(1,1,1,1), oma=c(1,1,1,1))
        } else {
            pdf(file.path(dir, sprintf('%s.pdf',channel)))
            par(mfrow=c(1,1))
            x <- c()
            y <- c()
        }
        for (individual in rep.individuals) {
            if (verbose) print(individual)
            i <- grep( channel, parameters(rep.fcs[[individual]][[day]][[dose]])@data$desc, ignore.case=T )
            d1 <- fun(rep.fcs[[individual]][[day]][[dose]]@exprs[,i])
            d2 <- fun(rep.down.fcs[[individual]][[day]][[dose]]@exprs[,i])
            if (class.f %in% c('density')) {
              x <- c(d1$x, d2$x)
              y <- c(d1$y, d2$y)
              #try and identify peaks
              plot(d1, xlab='', ylab='', main='', xaxt='n', yaxt='n', xlim=range(x), ylim=range(y), type='l')
              lines(d2, col='red')
              legend('topright', pch[pch$name==individual,'pch'], cex=2, bty='n')
            } else if (class.f %in% c('ecdf')) {
              plot(d1, xlab='', ylab='', main='', xaxt='n', yaxt='n')
              lines(d2, col='red')
              legend('topright', pch[pch$name==individual,'pch'], cex=2, bty='n')
            } else {
              x <- c(x, d1)
              y <- c(y, d2)
            }
        }
        if (!(class.f %in% c('density','ecdf'))) {
            r <- range(c(x, y))
            R2 <- round(R.squared(x, y), digits=3)
            plot(x, y, xlab='all', ylab='downsampled', xlim=r, ylim=r, pch=pch[,'pch'], cex=2, main=as.expression(bquote( r^2 == .(R2))), cex.main=2)
            #legend('topright', as.expression(bquote( r^2 == .(R2))), cex=2, bty='n')
            abline(b=1,a=0)
        }
        dev.off()
    } 
}

# Main scatterplot
marginal.2D.density <- function(fcs.data, fcs.data2, dir='~nikolas/IL2/Plots/marginal-density/', channels=c('FSC-A', 'SSC-A', 'CD25', 'CD4$', 'CD45RA', 'pSTAT5'), zoom=1) {
    dir.create(dir, recursive=T, showWarnings=F)
    comb <- t(combn(channels, 2))
    for (i in 1:nrow(comb)) {
        channels <- comb[i,]
        print(channels)
        d <- data.frame(getChannels(fcs.data, channels))
        d$col <- ifelse('col' %in% colnames(fcs.data@exprs), as.factor(fcs.data@exprs[,'col']), 'black') 
        #d$pct <- , fcs.data@exprs[,'pct'], .5) 
        if('pct' %in% colnames(fcs.data@exprs)) d$pct <- fcs.data@exprs[,'pct']
        else d$pct <- .5
        colnames(d) <- c('x', 'y', 'col', 'pct') 
        d2 <- data.frame(getChannels(fcs.data2, channels))
        d2$col <- ifelse('col' %in% colnames(fcs.data2@exprs), as.factor(fcs.data2@exprs[,'col']), 'green') 
        if('pct' %in% colnames(fcs.data2@exprs)) d2$pct <- fcs.data2@exprs[,'pct']
        else d2$pct <- .5
        colnames(d2) <- c('x', 'y', 'col', 'pct') 
        p1 <- ggplot(d, aes(x, y)) + 
          xlab(channels[1]) +
          ylab(channels[2]) +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_continuous(expand = c(0, 0)) +
          expand_limits(y = c(min(d$x) - .1*diff(range(d$x)), max(d$x) + .1*diff(range(d$x)))) +
          expand_limits(x = c(min(d$y) - .1*diff(range(d$y)), max(d$y) + .1*diff(range(d$y)))) +
          theme(plot.margin= unit(c(.2, .2, .5, .5), "lines"), plot.background=element_rect(fill='white')) +
          geom_point(data=d, colour='black', size=zoom*d$pct) +
          geom_density2d(data=d, colour='black') +
          #geom_point(data=d2, colour='green', size=d2$pct)
          geom_point(data=d2, colour='green', size=zoom*d2$pct) +
          geom_density2d(data=d2, colour='green')
          #geom_point(data=d3, colour='blue', size=.5) +
          #geom_point(data=d4, colour='green', size=d4$z) 
        # Horizontal marginal density plot - to appear at the top of the chart
        p2 <- ggplot(d, aes(x = x)) + 
          geom_density() +
          geom_density(data=d2, mapping=aes(x=x), colour='green') +
          #geom_density(data=d3, mapping=aes(x=x), colour='blue') +
          #geom_density(data=d4, mapping=aes(x=x), colour='green') +
          scale_x_continuous(expand = c(0, 0)) +
          expand_limits(x = c(min(d$y) - .1*diff(range(d$y)), max(d$y) + .1*diff(range(d$y)))) +
          theme(axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                plot.margin= unit(c(1, .2, -.5, 0.5), "lines"),
                plot.background=element_rect(fill='white')) 
        # Vertical marginal density plot - to appear at the right of the chart
        p3 <- ggplot(d, aes(x = y)) + 
          geom_density() +
          geom_density(data=d2, mapping=aes(x=y), colour='green') +
          #geom_density(data=d3, mapping=aes(x=y), colour='blue') +
          #geom_density(data=d4, mapping=aes(x=y), colour='green') +
          scale_x_continuous(expand = c(0, 0)) +
          expand_limits(x = c(min(d$x) - .1*diff(range(d$x)), max(d$x) + .1*diff(range(d$x)))) +
          coord_flip() +
          theme(axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                plot.margin= unit(c(.2, 1, 0.5, -.5), "lines"),
                plot.background=element_rect(fill='white')) 
        # Get the gtables
        gt1 <- ggplot_gtable(ggplot_build(p1))
        gt2 <- ggplot_gtable(ggplot_build(p2))
        gt3 <- ggplot_gtable(ggplot_build(p3)) 
        # Get maximum widths and heights for x-axis and y-axis title and text
        maxWidth = unit.pmax(gt1$widths[2:3], gt2$widths[2:3])
        maxHeight = unit.pmax(gt1$heights[4:5], gt3$heights[4:5]) 
        # Set the maximums in the gtables for gt1, gt2 and gt3
        gt1$widths[2:3] <- as.list(maxWidth)
        gt2$widths[2:3] <- as.list(maxWidth) 
        gt1$heights[4:5] <- as.list(maxHeight)
        gt3$heights[4:5] <- as.list(maxHeight) 
        # Combine the scatterplot with the two marginal boxplots
        # Create a new gtable
        gt <- gtable(widths = unit(c(7, 2), "null"), height = unit(c(2, 7), "null")) 
        # Instert gt1, gt2 and gt3 into the new gtable
        gt <- gtable_add_grob(gt, gt1, 2, 1)
        gt <- gtable_add_grob(gt, gt2, 1, 1)
        gt <- gtable_add_grob(gt, gt3, 2, 2) 
        # And render the plot
        pdf(file.path(dir, sprintf('%s.pdf',paste(channels, collapse='-'))))
        print(grid.newpage())
        print(grid.draw(gt))
        print(grid.rect(x = 0.5, y = 0.5, height = .995, width = .995, default.units = "npc", gp = gpar(col="black", fill = NA, lwd = 1)) ) 
        dev.off()
    }
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

colMedians <- function(x) apply(x, 2, median)

compute.cluster.stats <- function(fcs.data, cluster.assignment, channels=c('FSC-A', 'SSC-A', 'CD25', 'CD4$', 'CD45RA', 'pstat5')) {
    fcs.data <- cbind( getChannels(fcs.data, channels), cluster.assignment )
    total.count <- nrow(fcs.data)
    head(x <- as.matrix(do.call('rbind', by( fcs.data, cluster.assignment, function(x) data.frame(t(colMedians(x)),count=nrow(x),pct=100*nrow(x)/total.count)))))
    colnames(x) <- c(gsub('\\$', '', channels), 'cluster.assignment', 'count', 'pct')
    build.flowFrame(x)
}


# Assumption bottom 1% and top 1% do not react to IL2
# Points which are invariant across doses, should be invariant across days.
quantile.normalise <- function(resting, stimulated, q0=.01, q1=.99, channel='pSTAT5') {
    conc0 <- getChannels(resting, channel)
    conc1 <- getChannels(stimulated, channel)
    q.conc0 <- quantile(conc0, probs=c(q0,q1))
    q.conc1 <- quantile(conc1, probs=c(q0,q1))
    m <- lm(q.conc0 ~ q.conc1)
    return( function(x) cbind(1,x)%*%coefficients(m) )
}

quantile.normalize <- function(f0, f1, quantiles=c(.25,.5,.75), channel=NULL) {
    if (is.null(channel)) x0 <- f0
    else x0 <- getChannels(f0, channel)
    if (is.null(channel)) x1 <- f1
    else x1 <- getChannels(f1, channel)
    q.x0 <- quantile(x0, probs=quantiles)
    q.x1 <- quantile(x1, probs=quantiles)
    m <- lm(q.x0 ~ q.x1)
    #f1 <- setChannels(f1, channel, x1.trans)
    #return(f1)
    return( cbind(1,x1)%*%coefficients(m) )
}

peak.normalize <- function(f0, f1, k=2, channel=NULL) {
    if (is.null(channel)) x0 <- f0
    else x0 <- getChannels(f0, channel)
    if (is.null(channel)) x1 <- f1
    else x1 <- getChannels(f1, channel)
    q.x0 <- pam(x0[sample(x0,1000)], k=k)$medoids
    q.x1 <- pam(x1[sample(x1,1000)], k=k)$medoids
    m <- lm(q.x0 ~ q.x1)
    #f1 <- setChannels(f1, channel, x1.trans)
    #return(f1)
    return( cbind(1,x1)%*%coefficients(m) )
}

my.normalize <- function(f0, f1, channel=NULL) {
    if (is.null(channel)) x0 <- f0
    else x0 <- getChannels(f0, channel)
    if (is.null(channel)) x1 <- f1
    else x1 <- getChannels(f1, channel)
    q.x0 <- c(mean(x0), Mode(x0), median(x0))
    q.x1 <- c(mean(x1), Mode(x1), median(x1))
    m <- lm(q.x0 ~ q.x1)
    #f1 <- setChannels(f1, channel, x1.trans)
    #return(f1)
    return( cbind(1,x1)%*%coefficients(m) )
}

quantile.transforms <- function(flow.data) {
    transforms <- list()
    for ( f in flow.data[grep('\\.0U', names(flow.data))] ) {
        individual <- getIndividual(f)
        day <- getDate(f)
        #print(paste(day, individual, sep='.'))
        resting <- f
        stimulated <- flow.data[[ grep( paste(day, individual, '1000U', sep='.'), names(flow.data) )]]
        transforms[[paste(day, individual, sep='.')]] <- quantile.normalise(resting, stimulated)
    }
    transforms
}


plot.dose.response <- function(flow.data, q.transforms, q0=.01, q1=.99, dose='1000U', dir='~nikolas/dunwich/Projects/IL2/Plots/ungated/dose-response/') { # e0=baseline.conc.ecdf; e1=high.conc.ecdf
    abc <- data.frame()
    dir.create(dir, recursive=T, showWarnings=F)
    for (individual in unique(simplify2array(strsplit(names(flow.data), '\\.'))[2,])) {
          #png(file.path(dir, sprintf('%s.png',individual)), width=400, height=400)
          pdf(file.path(dir, sprintf('%s.pdf',individual)))
          par(mfrow=c(2,2), mar=c(2,3,2,2), mgp=c(2,1,0), las=1)
       f <- function(normalised) {
           abc <- data.frame(individual, normalised)
            for (day in sapply(flow.data[grep(paste(individual, '0U', sep='.'), names(flow.data))] , getDate)) {
                cat(individual, day, normalised, '\n')
                x <- seq(-.5, 4, length.out=20000)
                e0 <- ecdf(getChannels(flow.data[[grep(paste(day, individual, '0U', sep='.'), names(flow.data))]], 'pstat5'))(x)
                if (normalised) 
                    e1 <- ecdf(q.transforms[[paste(day,individual,sep='.')]](getChannels(flow.data[[grep(paste(day, individual, dose, sep='.'), names(flow.data))]], 'pstat5')))(x)
                else
                    e1 <- ecdf(getChannels(flow.data[[grep(paste(day, individual, dose, sep='.'), names(flow.data))]], 'pstat5'))(x)
                #plot(x, e0, xlab='pStat5 signal', col='white', main=sprintf("Day %s",day),ylab="", xlim=c(-.5,3))
                plot(x, e0, col='white', xlim=c(-.5,3), xlab='', ylab='')
                #axis(side=1)
                #axis(side=2)
                lines(x, e0, col=ifelse(normalised, 'pink', 'grey'), lwd=2)
                lines(x, e1, col=ifelse(normalised, 'red', 'black'), lwd=2)
                ecdf.diff <- e0-e1
                lines(x, ecdf.diff, lty=2, col='black',lwd=2)
                addpoints(x, e0,q0,q1,0)
                addpoints(x, e1,q0,q1,1)
                #text <- c("0U", "1000U")
                #legend("right",lwd=rep(2,length(text)),col=ifelse(normalised, c('green','darkgreen'), c('pink','red')),legend=text)
                ecdf.diff.area <- abs(sum(ecdf.diff)*diff(x)[1])
                abc.pct <- round(ecdf.diff.area,digits=2)
                abc <- cbind(abc, abc.pct)
                #title(
                text(0,1,sprintf('%.2f', abc.pct)) 
                }
            return(abc)
       }
       abc <- rbind(abc, f(FALSE))
       #abc <- rbind(abc, f(TRUE))
      dev.off()
    }
    colnames(abc) <- c('individual', 'normalised', 'abc.pct.day1', 'abc.pct.day2')
    return(abc)
}

addpoints <- function(x, y,q0=0.01,q1=0.99,wh=0) {
  iq <- which.min(abs(y-q0))
  iQ <- which.min(abs(y-q1))
  #cols <- c('green', 'darkgreen')
  fills <- c("#ff000099","#0000ff99")
  points(x[c(iq,iQ)],y[c(iq,iQ)],pch=25-wh,col='black',bg=fills[wh+1],cex=1.5)
}

plot.dose.response.agreement <- function(d,dir='~nikolas/IL2/Plots/dose-response/', main='') {
    #d$abc.pct.day1 <- d$abc.pct.day1/100
    #d$abc.pct.day2 <- d$abc.pct.day2/100
    dir.create(dir, recursive=T, showWarnings=F)
    pdf(file.path(dir, 'agreement.pdf'))
    plot(d$abc.pct.day1, d$abc.pct.day2, pch=d$pch,
         col=as.numeric(d$normalised)+1,
         xlim=range(as.numeric(c(d$abc.pct.day1, d$abc.pct.day2))),
         ylim=range(as.numeric(c(d$abc.pct.day1,d$abc.pct.day2))),
         xlab='Area Between Curves Day 1',
         ylab='Area Between Curves Day 2',
         main=main)
    abline(b=1, a=0)
    r <- round(cor(d[!d$normalised, 'abc.pct.day1'], d[!d$normalised, 'abc.pct.day2'])**2, digits=3)
    r.norm <- round(cor(d[d$normalised, 'abc.pct.day1'], d[d$normalised, 'abc.pct.day2'])**2, digits=3)
    #legend('topleft', legend=c(as.expression(bquote(r^2 == .(r))), as.expression(bquote(r^2 == .(r.norm)))), text.col=c('black', 'red'))
    legend('topleft', legend=as.expression(bquote(r^2 == .(r))))
    dev.off()
}

