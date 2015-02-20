library(spade)
source('~nikolas/bin/FCS/fcs.R')

downsample.FCS <- function(fcs, comp=TRUE, verbose=FALSE, cols=c('FSC-A', 'SSC-A'), ...) {
   #if (class(fcs) != 'flowFrame' && file.exists(fcs)) 
   fcs <- read.FCS(fcs, comp, verbose, ...)
    #assume it's already a flow frame

   in_data <- fcs@exprs
   density <- SPADE.density(in_data[,cols], kernel_mult=5.0, apprx_mult=1.5, med_samples=2000)

   desired_samples <- NULL

   # exclude_pctile,target_pctile
   boundary <- quantile(density,c(0.01, 0.05),names=FALSE)
    
   in_data <- cbind(in_data, density)
   out_data <- subset(in_data, density > boundary[1]) # Exclusion    
   desired_samples <- 20000
   density <- out_data[,'density']
    if (desired_samples < nrow(out_data)) {
		# Need to find target density such there are approximately desired_samples
		# remaining after downsampling. To do so we solve for the density such that
		# the sum of samples below that density plus the expected value of
		# samples retained above that density equals approximately the desired
		# number of samples
		density_s <- sort(density)
		cdf       <- rev(cumsum(1.0/rev(density_s)))
		
		# Default solution if target density smaller than any present
		boundary <- desired_samples/cdf[1] 
		if (boundary > density_s[1]) {  # Boundary actually falls amongst densities present
			targets <- (desired_samples-1:length(density_s)) / cdf 
            #looks for first FALSE
			boundary <- targets[which.min(targets-density_s > 0)]
		}
		out_data  <- subset(out_data,boundary/density > runif(length(density)))
    }
    out_data <- out_data[,-which(colnames(out_data)=='density')]

    flowFrame(out_data,parameters(fcs),description=description(fcs))
}


#§downsample


