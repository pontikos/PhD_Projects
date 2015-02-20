source('~nikolas/bin/FCS/fcs.R')
suppressPackageStartupMessages(library(mixtools))
suppressPackageStartupMessages(library(iterators))
suppressPackageStartupMessages(library(RANN))

split.response <- function(d, i) {
    d <- as.data.frame(d)
    #diff.PSTAT5.1
    if (i < 2) return(list(d=d))
    if (nrow(d)<10) return(list(d=d))
    print(diff.pstat5 <- paste('diff','PSTAT5',i,sep='.'))
    x <- d[,diff.pstat5]
    res <- kmeans(x, centers=c(0,1))
    mu <- as.numeric(res$centers)
    lambda <- as.numeric(prop.table(table(res$cluster)))
    sigsqrd <- as.numeric(tapply(x,res$cluster,var))
    m <- mixtools::normalmixEM2comp(x,mu=mu,lambda=lambda,sigsqrd=sigsqrd)
    i1 <- which( (d[,diff.pstat5] > 0) & (m$posterior[,1]>.99))
    d1 <- d[i1,]
    i2 <- which( (d[,diff.pstat5] > 0) & (m$posterior[,2]>.99))
    d2 <- d[i2,]
    return(list(i1=i1, d1=split.response(d1, i-1), i2=i2, d2=split.response(d2, i-1), d=d, m=m))
}

plot.tree <- function(X,fcs.data,CLR,cells,cols=c('purple','pink','lightblue')) {
    figure.labels <- iter(paste(letters,')',sep=''))
    layout(matrix(1:15,nrow=3,ncol=5,byrow=T),widths=rep(4,5),heights=rep(1,3))
    par(mai=c(1,.25,.1,.12),cex.lab=1)
    plot.new()
    plot.new()
    #arrows(1,.5,.5,0)
    my.plot.mixEM(X$m,which=2,main2='',xlab2='pSTAT5 response at 1000U')
    title(nextElem(figure.labels), adj=0)
    for (i in 1:length(cells)) lines(density(fcs.data[as.logical(CLR[,cells[i]]),'diff.PSTAT5.4']),col=cols[i],lwd=2)
    plot.new()
    #arrows(0,.5,.5,0)
    plot.new()
    #
    plot.new()
    my.plot.mixEM(X$d1$m,which=2,main2='low',xlab2='pSTAT5 response at 10U')
    title(nextElem(figure.labels), adj=0)
    graphics::box(lty='solid',col='red')
    plot.new()
    legend('center','10 units',cex=1,bty='n')
    my.plot.mixEM(X$d2$m,which=2,main2='high',xlab2='pSTAT5 response at 10U')
    title(nextElem(figure.labels), adj=0)
    for (i in 1:length(cells)) lines(density(fcs.data[as.logical(CLR[,cells[i]]),'diff.PSTAT5.3']),col=cols[i],lwd=2)
    graphics::box(lty='solid',col='green')
    plot.new()
    legend('center',ncol=1,fill=cols[1:length(cells)],legend=cells,bty='n',cex=1)
    # d)
    tryCatch({
      my.plot.mixEM(X$d1$d1$m,which=2,main2='low',xlab2='pSTAT5 response at 0.1U')
    }, warning = function(w) {
        #plot.new()
      plot(normalised.density(X$d1$d[,'diff.PSTAT5.2']),main='low',xlab='pSTAT5 response at 0.1U')
    }, error = function(e) {
      #plot.new()
      plot(normalised.density(X$d1$d[,'diff.PSTAT5.2']),main='low',xlab='pSTAT5 response at 0.1U')
    }, finally = {
      graphics::box(lty='solid',col='red')
      title(nextElem(figure.labels), adj=0)
    })
    # e)
    tryCatch({
    my.plot.mixEM(X$d1$d2$m,which=2,main2='high',xlab2='pSTAT5 response at 0.1U')   
    }, warning = function(w) {
        plot.new()
    }, error = function(e) {
        plot.new()
    }, finally = {
      graphics::box(lty='solid',col='green')
      title(nextElem(figure.labels), adj=0)
    })
    plot.new()
    legend('center','0.1 units', cex=1, bty='n')
    # f)
    tryCatch({
      my.plot.mixEM(X$d2$d1$m,which=2,main2='low',xlab2='pSTAT5 response at 0.1U')
    }, warning = function(w) {
#        plot.new()
      plot(normalised.density(X$d2$d[,'diff.PSTAT5.2']),main='low',xlab='pSTAT5 response at 0.1U')
    }, error = function(e) {
        #plot.new()
      plot(normalised.density(X$d2$d[,'diff.PSTAT5.2']),main='low',xlab='pSTAT5 response at 0.1U')
    }, finally = {
      graphics::box(lty='solid',col='red')
      title(nextElem(figure.labels), adj=0)
    })
    # g)
    tryCatch({
      my.plot.mixEM(X$d2$d2$m,which=2,main2='high',xlab2='pSTAT5 response at 0.1U')
      for (i in 1:length(cells)) lines(density(fcs.data[as.logical(CLR[,cells[i]]),'diff.PSTAT5.2']),col=cols[i],lwd=2)
    }, warning = function(w) {
        plot.new()
    }, error = function(e) {
        plot.new()  
    }, finally = {
      graphics::box(lty='solid',col='green')
      title(nextElem(figure.labels), adj=0)
    })
}



my.plot.mixEM <- function (x, whichplots = 1, loglik = 1 %in% whichplots, density = 2 %in% 
            whichplots, xlab1 = "Iteration", ylab1 = "Log-Likelihood", 
          main1 = "Observed Data Log-Likelihood", col1 = 1, lwd1 = 2, 
          xlab2 = NULL, ylab2 = NULL, main2 = NULL, col2 = NULL, lwd2 = 2, 
          alpha = 0.05, marginal = FALSE, ...) 
{
  def.par <- par(ask = (loglik + density > 1), "mar")
  mix.object <- x
  if (!inherits(mix.object, "mixEM")) 
    stop("Use only with \"mixEM\" objects!")
  if (loglik) {
    plot(mix.object$all.loglik, xlab = xlab1, ylab = ylab1, 
         main = main1, type = "l", lwd = lwd1, col = col1, 
         ...)
  }
  if (density) {
    if (mix.object$ft == "logisregmixEM") {
      if (ncol(mix.object$x) != 2) {
        stop("The predictors must have 2 columns!")
      }
      if (sum((mix.object$y == 1) + (mix.object$y == 0)) != 
            length(mix.object$y)) {
        stop("The response must be binary!")
      }
      k = ncol(mix.object$beta)
      x = mix.object$x[, 2]
      if (is.null(main2)) {
        main2 <- "Most Probable Component Membership"
      }
      if (is.null(xlab2)) {
        xlab2 <- "Predictor"
      }
      if (is.null(ylab2)) {
        ylab2 <- "Response"
      }
      if (is.null(col2)) {
        col2 <- 2:(k + 1)
      }
      plot(x, mix.object$y, main = main2, xlab = xlab2, 
           ylab = ylab2, col = col2[apply(mix.object$posterior, 
                                          1, which.max)], ...)
      a = cbind(x, mix.object$y)
      a = a[order(a[, 1]), ]
      for (i in 1:k) {
        lines(a[, 1], plogis(mix.object$beta[1, i] + 
                               mix.object$beta[2, i] * a[, 1]), col = col2[i])
      }
    }
    if (mix.object$ft == "normalmixEM") {
      k <- ncol(mix.object$posterior)
      x <- sort(mix.object$x)
      a <- hist(x, plot = FALSE)
      maxy <- max(max(a$density), 0.3989 * mix.object$lambda/mix.object$sigma)
      if (is.null(main2)) {
        main2 <- "Density Curves"
      }
      if (is.null(xlab2)) {
        xlab2 <- "Data"
      }
      if (is.null(col2)) {
        col2 <- 2:(k + 1)
      }
      hist(x, prob = TRUE, main = main2, xlab = xlab2,  ylim = c(0, maxy), border='white', ...)
      if (length(mix.object$mu) == 1) {
        arbvar <- TRUE
        mix.object$sigma <- mix.object$scale * mix.object$sigma
        arbmean <- FALSE
      }
      if (length(mix.object$mu) == k && length(mix.object$sigma) ==  1) {
        arbmean <- TRUE
        arbvar <- FALSE
      }
      if (length(mix.object$sigma) == k && length(mix.object$mu) ==  k) {
        arbmean <- TRUE
        arbvar <- TRUE
      }
      for (i in 1:k) {
        lines(x, mix.object$lambda[i] * dnorm(x,
                          mean = mix.object$mu[i *  arbmean + (1 - arbmean)],
                          sd = mix.object$sigma[i *  arbvar + (1 - arbvar)]), col = col2[i], lwd = lwd2)
      }
    }
    if (mix.object$ft == "repnormmixEM") {
      x <- as.vector(as.matrix(x))
      k <- ncol(mix.object$posterior)
      x <- sort(mix.object$x)
      a <- hist(x, plot = FALSE)
      maxy <- max(max(a$density), 0.3989 * mix.object$lambda/mix.object$sigma)
      if (is.null(main2)) {
        main2 <- "Density Curves"
      }
      if (is.null(xlab2)) {
        xlab2 <- "Data"
      }
      if (is.null(col2)) {
        col2 <- 2:(k + 1)
      }
      hist(x, prob = TRUE, main = main2, xlab = xlab2, 
           ylim = c(0, maxy), ...)
      if (length(mix.object$mu) == 1) {
        arbvar <- TRUE
        mix.object$sigma = mix.object$scale * mix.object$sigma
        arbmean <- FALSE
      }
      if (length(mix.object$mu) == k && length(mix.object$sigma) == 
            1) {
        arbmean <- TRUE
        arbvar <- FALSE
      }
      if (length(mix.object$sigma) == k && length(mix.object$mu) == 
            k) {
        arbmean <- TRUE
        arbvar <- TRUE
      }
      for (i in 1:k) {
        lines(x, mix.object$lambda[i] * dnorm(x, mean = mix.object$mu[i * 
                                                                        arbmean + (1 - arbmean)], sd = mix.object$sigma[i * 
                                                                                                                          arbvar + (1 - arbvar)]), col = col2[i], lwd = lwd2)
      }
    }
    if (mix.object$ft == "regmixEM.mixed") {
      x.1 = mix.object$x
      n = sum(sapply(x.1, nrow))
      x.1.sum = sum(sapply(1:length(x.1), function(i) length(x.1[[i]][, 
                                                                      1])))
      if (x.1.sum == n) {
        x = lapply(1:length(x.1), function(i) matrix(x.1[[i]][, 
                                                              -1], ncol = 1))
      }
      else {
        x = x.1
      }
      post.beta(x = x, y = mix.object$y, p.beta = mix.object$posterior.beta, 
                p.z = mix.object$posterior.z)
    }
    if (mix.object$ft == "mvnormalmixEM") {
      x = mix.object$x
      if (ncol(x) != 2) {
        stop("The data must have 2 columns!")
      }
      post = apply(mix.object$posterior, 1, which.max)
      k <- ncol(mix.object$posterior)
      if (is.list(mix.object$sigma)) {
        sigma = mix.object$sigma
      }
      else {
        sigma = lapply(1:k, function(i) mix.object$sigma)
      }
      if (is.list(mix.object$mu)) {
        mu = mix.object$mu
      }
      else {
        mu = lapply(1:k, function(i) mix.object$mu)
      }
      if (is.null(xlab2)) {
        xlab2 <- "X.1"
      }
      if (is.null(ylab2)) {
        ylab2 <- "X.2"
      }
      if (is.null(col2)) {
        col2 <- 2:(k + 1)
      }
      if (marginal == FALSE) {
        if (is.null(main2)) {
          main2 <- "Density Curves"
        }
        plot(x, col = col2[post], xlab = xlab2, ylab = ylab2, 
             main = main2, ...)
        lapply(1:k, function(i) points(mu[[i]][1], mu[[i]][2], 
                                       pch = 19))
        for (i in 1:k) {
          for (j in 1:length(alpha)) {
            ellipse(mu = mu[[i]], sigma = sigma[[i]], 
                    alpha = alpha[j], col = col2[i])
          }
        }
      }
      else {
        if (is.null(main2)) {
          main2 <- ""
        }
        x <- mix.object$x[, 1]
        y <- mix.object$x[, 2]
        xhist <- hist(x, plot = FALSE)
        yhist <- hist(y, plot = FALSE)
        top <- max(c(xhist$counts, yhist$counts))
        xrange <- range(x)
        yrange <- range(y)
        nf <- layout(matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE), 
                     c(4, 1), c(1, 4), TRUE)
        layout.show(nf)
        par(mar = c(3, 3, 1, 1))
        plot(mix.object$x[, 1], mix.object$x[, 2], col = col2[post], 
             xlab = xlab2, ylab = ylab2, main = main2, ...)
        lapply(1:k, function(i) points(mu[[i]][1], mu[[i]][2], 
                                       pch = 19))
        for (i in 1:k) {
          for (j in 1:length(alpha)) {
            ellipse(mu = mu[[i]], sigma = sigma[[i]], 
                    alpha = alpha[j], col = (i + 1))
          }
        }
        par(mar = c(0, 3, 1, 1))
        barplot(xhist$counts, axes = FALSE, ylim = c(0, 
                                                     top), space = 0, ...)
        par(mar = c(3, 0, 1, 1))
        barplot(yhist$counts, axes = FALSE, xlim = c(0, 
                                                     top), space = 0, horiz = TRUE, ...)
      }
    }
    if (mix.object$ft == "regmixEM") {
      if (ncol(mix.object$x) != 2) {
        stop("The predictors must have 2 columns!")
      }
      post <- apply(mix.object$posterior, 1, which.max)
      k <- ncol(mix.object$posterior)
      x <- mix.object$x[, 2]
      y <- mix.object$y
      n <- length(y)
      if (is.null(main2)) {
        main2 <- "Most Probable Component Membership"
      }
      if (is.null(xlab2)) {
        xlab2 <- "Predictor"
      }
      if (is.null(ylab2)) {
        ylab2 <- "Response"
      }
      if (is.null(col2)) {
        col2 <- 2:(k + 1)
      }
      plot(x, y, main = main2, xlab = xlab2, ylab = ylab2, 
           type = "n", ...)
      a = cbind(mix.object$x[, 2], mix.object$y, post)
      for (i in 1:k) {
        xy = subset(cbind(a, mix.object$posterior[, i]), 
                    a[, 3] == i)[, -3]
        xy = matrix(xy, ncol = 3)
        points(xy[, 1], xy[, 2], col = col2[i])
        if (is.matrix(mix.object$beta) == FALSE) {
          abline(coef = mix.object$beta)
          beta = matrix(mix.object$beta, ncol = k, nrow = 2)
        }
        else {
          abline(coef = mix.object$beta[, i], col = col2[i])
          beta = mix.object$beta
        }
        out = lm(y ~ x, weights = mix.object$posterior[, 
                                                       i])
        fit = beta[1, i] + beta[2, i] * x
        out.aov = anova(out)
        MSE = out.aov$Mean[2]
        xy.f = cbind(x, y, fit)
        xy.sort = xy.f[order(xy.f[, 1]), ]
        x.new = seq(from = min(x), to = max(x), length.out = 100)
        y.new = beta[1, i] + beta[2, i] * x.new
        s.h <- sqrt(MSE * (1/n + (x.new - mean(xy.sort[, 
                                                       1]))^2/var(xy.sort[, 1])/(n - 1)))
        for (j in 1:length(alpha)) {
          W = sqrt(qf(1 - alpha[j], 2, n - 2))
          upper = y.new + W * s.h
          lower = y.new - W * s.h
          lines(x.new, upper, col = (i + 1))
          lines(x.new, lower, col = (i + 1))
        }
      }
    }
  }
  par(def.par)
}


