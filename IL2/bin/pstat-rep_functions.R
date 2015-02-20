daycol <- function(day) {
  #if(day==1) return(c("blue","darkblue"))
  return(c("green","darkgreen"))
}


f <- function(conc0, conc1, q0, q1) {
    q.conc0 <- quantile(conc0, probs=c(q0,q1))
    q.conc1 <- quantile(conc1, probs=c(q0,q1))
    m <- lm(q.conc0 ~ q.conc1)
    cbind(1,conc1)%*%coefficients(m)
}

myp <- function(d0,d1, day=1, q0=0.01, q1=0.99, adjust=FALSE, new=FALSE) { # e0=baseline.conc.ecdf; e1=high.conc.ecdf

  conc0 <- trans(d0@exprs[,7])
  conc1 <- trans(d1@exprs[,7])
  if(adjust) {
    conc1 <- f(conc0, conc1, q0, q1)
  }
  e0 <- ecdf(conc0)
  e1 <- ecdf(conc1)

  if(new) {
    par(mar=c(3,2,1,0),mgp=c(2,1,0))
    plot(e0, xlab='pStat5 signal', col='white', main=sprintf("Day %s",day),ylab="")
  }
  cols <- daycol(day)
  
  ## x <- seq(0,1,.01)
  ## inv.baseline.conc <- function(x) sapply(x, function(x) inverse(e0, -1, 5)(x)$root)
  ## inv.high.conc <- function(x) sapply(x, function(x) inverse(e1, -1, 5)(x)$root)
  ## inv.ecdf.diff <- function(x) return(inv.high.conc(x)-inv.baseline.conc(x))
  ## lines(inv.ecdf.diff(x), x, col=cols[1], lty=2, lwd=2)
  ## 
  lines(e0, col=cols[1], lwd=2)
  lines(e1, col=cols[2], lwd=2)

  ## inv.ecdf.diff.area <- abs(sum(inv.ecdf.diff(x)))
  ## fi.auc <- round(inv.ecdf.diff.area,digits=2)
  ## if(day==1) {
  ##   text(0, 1, fi.auc, col=cols[1], adj=c(0.5,0))
  ## }
  #  x <- seq(min(baseline.conc), max(baseline.conc), .1)
  x <- knots(e0)
  ecdf.diff <- function(x) return(-e1(x)+e0(x))
  lines(x, (ecdf.diff(x)), lty=2, col=cols[1],lwd=2)
  ecdf.diff.area <- abs(sum(ecdf.diff(x)))
  pct.auc <- round(ecdf.diff.area,digits=2)
  ## if(day==1) {
    ## text(4, 0, pct.auc, col=cols[1], adj=c(0,0))
  ## } else {
  ##   text(4, 0, pct.auc, col=cols[1], adj=c(0,1))
  ## }
  addpoints(e0,q0,q1,day,0)
  addpoints(e1,q0,q1,day,1)
}

addpoints <- function(e,q0=0.01,q1=0.99,day=1,wh=0) {
  x <- knots(e)
  y <- e(x)
  iq <- which.min(abs(y-q0))
  iQ <- which.min(abs(y-q1))
  cols <- daycol(day)
  fills <- c("#ff000099","#0000ff99")
  points(x[c(iq,iQ)],y[c(iq,iQ)],pch=25-wh,col=cols[wh+1],bg=fills[wh+1],cex=1.5)
}

addlegend <- function(days=c(1,2)) {
  cols <- unlist(lapply(days, daycol))
  if(length(days)==1) {
    text <- c("0U", "1000U")
  } else {
    text <- c("Day 1, 0U", "Day 1, 1000U", "Day 2, 0U", "Day 2, 1000U")
  }
  legend("right",lwd=rep(2,length(text)),col=cols,legend=text)
}


make.plots <- function(individual, cell.type=NULL)  {
  if(is.null(cell.type)) cell.type <- ''
  else cell.type <- paste(cell.type,'.',sep='')
  print(cell.type)
  day1.1 <- read.FCS(file.path(dir, individual, 'day1', sprintf('0U.%sfcs',cell.type)))
  day1.2 <- read.FCS(file.path(dir, individual, 'day1', sprintf('1000U.%sfcs',cell.type))) 
  day2.1 <- read.FCS(file.path(dir, individual, 'day2', sprintf('0U.%sfcs',cell.type)))
  day2.2 <- read.FCS(file.path(dir, individual, 'day2', sprintf('1000U.%sfcs',cell.type)))
  png(sprintf("%s/raw.%spng",individual,cell.type), width = 800, height = 400)
  par(mfrow=c(1,2))
  myp(day1.1, day1.2, day=1, new=TRUE)
  #addlegend(1)
  myp(day2.1, day2.2, day=2,new=TRUE)
  #addlegend(2)
  dev.off()
  png(sprintf("%s/norm.%spng",individual, cell.type), width = 800, height = 400)
  par(mfrow=c(1,2))
  myp(day1.1, day1.2, day=1, new=TRUE, adjust=TRUE)
  #addlegend(1)
  myp(day2.1, day2.2, day=2,new=TRUE,  adjust=TRUE)
  #addlegend(2)
  dev.off()

}
