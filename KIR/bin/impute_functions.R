## data read
read.snp <- function(f, cols, train=TRUE) {
  cat("reading",f,"\t\t\t\r")
  data <- read.table(f,sep="\t", col.names=cols, row.names=1)
  data <- cbind(data, geno[rownames(data), c("collection", "KIR.3DS1", "KIR.3DL1","KIR")])
  data$f <- sub("lrr-data/","",f)
  if(train) {
    return(subset(data,!is.na(KIR.3DS1)))
  } else {
    return(data)
  }
}
read.snp.mi <- function(f, cols) {
  cat("reading",f,"\t\t\t\r")
  data <- read.table(f,sep="\t", col.names=cols, row.names=1)
  data$f <- sub("lrr-data/","",f)
  return(data)
}



## data prep
make.wide <- function(data, cols, mi=FALSE) {
  wide <- reshape(data[,c("sample","SNP_Name", cols)],
                  v.names=cols,
                  timevar="SNP_Name",
                  idvar="sample",
                  direction="wide")
  if(mi) { return(wide) }
  cbind(wide,
        KIR=as.factor(geno[as.character(wide$sample),"KIR"]),
        KIR.3DS1=as.factor(geno[as.character(wide$sample),"KIR.3DS1"]),
        KIR.3DL1=as.factor(geno[as.character(wide$sample),"KIR.3DL1"]))
}

## loss functions
L.KIR <- function(true,pred) {
#  print(table(true,pred))
  sum(true!=pred)/length(true)
}
L.DS1 <- function(true,pred) {
  true <- substr(true,3,3)
  pred <- substr(pred,3,3)
 # print(table(true,pred))
  sum(true!=pred)/length(true)
}
L.DL1 <- function(true,pred) {
  true <- substr(true,4,4)
  pred <- substr(pred,4,4)
 # print(table(true,pred))
  sum(true!=pred)/length(true)
}

loss <- function(y,pred,samples, loss.function) {
  ncv <- max(samples)
  nit <- ncol(pred)
  errors <- matrix(NA,ncv,ncol(pred))
  for(icv in 1:ncv) {
    isample <- which(samples==icv) 
    errors[icv, ] <- sapply(1:nit, function(i) {
      loss.function(y[isample], pred[isample,i]) })
  }
  return(errors)
}

varycv.knn <- function(wdata, train=0.5, KIR="KIR", k=1:20) {
  kcols <- grep("KIR",colnames(wdata))
  wdata <- wdata[complete.cases(wdata[,kcols]),]
  if(length(train)==1) {
    n <- nrow(wdata)
    i <- sample(1:nrow(wdata),round(train*n))
  }
  x.train <- wdata[i,-c(1,kcols)]
  x.test <- wdata[-i,-c(1,kcols)]
  y.train <- wdata[i,KIR]
  y.test <- wdata[-i,KIR]
  e <- numeric(length(k))
  for(i in 1:length(k)) {
   m <-  knn(train=x.train, test=x.test, cl=y.train, k=k[i])
   tt <- table(m,y.test)
   e[i] <- (sum(tt) - sum(diag(tt)))/sum(tt)    
  }
  return(data.frame(k=k,train=train,e=e))
}

mycv.knn <- function(wdata, KIR="KIR", k=1:20,trim=TRUE,trim.pattern="rs592645|60034052|60056721",
                     do.plot=TRUE) {
  kcols <- grep("KIR",colnames(wdata))
  x <- wdata[,-c(1,kcols)]
  if(trim)
    x <- x[,grep(trim.pattern,colnames(x))]
  message("predicting using ",ncol(x)," columns:\n",paste(colnames(x),collapse=" "))
  y <- wdata[,KIR]
  e <- numeric(length(k))
  for(i in 1:length(k)) {
    m <- knn.cv(train=x, cl=y, k=k[i], prob=TRUE)
    tt <- table(m,y)
    e[i] <- (sum(tt) - sum(diag(tt)))/sum(tt)    
  }
  if(do.plot)
    plot(k,e,type="b")
  return(data.frame(k=k,e=e))
}

my.knn <- function(wdata, KIR="KIR", k=5, trim.pattern="rs592645|60034052|60056721") {
  kcols <- grep("KIR",colnames(wdata))

  ## no missing allowed
  use <- complete.cases(wdata[,-kcols])
  wdata <- wdata[use,]

  samples <- !is.na(wdata[,KIR])
  x <- wdata[,-c(1,kcols)]
  x <- x[,grep(trim.pattern,colnames(x))]
  x.train <- x[samples,]
  x.test <- x ## all samples
  y <- wdata[samples,KIR]

  m <-  knn(train=x.train, test=x.test, cl=y, k=k, prob=TRUE)
  return(data.frame(sample=wdata$sample,
                    KIR.niko=as.character(wdata[,KIR]),
                    KIR.pred=as.character(m),
                    prob=attributes(m)$prob,
                    stringsAsFactors=FALSE))
}

## 10-fold CV nnet
mycv <- function(wdata,ncv=10,nit=10,KIR="KIR") {
  kcols <- grep("KIR",colnames(wdata))
  x <- wdata[,-c(1,kcols)]
  y <- wdata[,KIR]
  samples <- sample(rep(1:ncv,length=nrow(wdata)))
  mi.data <- matrix(as.character(NA),nrow(x),nit)
  rownames(mi.data) <- as.character(wdata[,1])
  for(icv in 1:ncv) {
    isample <- which(samples==icv)
    mnet <- multinom(y[-isample] ~ ., data=x[-isample,], maxit=10000)
    mpred <- t(apply(predict(mnet, newdata=x[isample,], type="probs"),1,cumsum))
    for(i in 1:nit) {
      cat(".")
      U <- runif(nrow(mpred))
      mi.data[isample,i] <- colnames(mpred)[ apply(U < mpred, 1, function(x) which(x)[1]) ]
    }
  }
  return(list(y=y,
              mi.data=mi.data,
              errors.KIR=loss(y, mi.data, samples, L.KIR),
              errors.3DS1= loss(y, mi.data, samples, L.DS1),
              errors.3DL1=loss(y, mi.data, samples, L.DL1)))
}

mynet <- function(wdata,nit=10,KIR="KIR") {
  kcols <- grep("KIR",colnames(wdata))

  ## no missing allowed
  use <- complete.cases(wdata[,-kcols])
  wdata <- wdata[use,]

  samples <- !is.na(wdata[,KIR])
  x.train <- wdata[samples,-c(1,kcols)]
  x.test <- wdata[,-c(1,kcols)] ## all samples
  y <- wdata[samples,KIR]
  mi.data <- matrix(as.character(NA),nrow(x.test),nit)
  mnet <- multinom(y ~ ., data=x.train, maxit=10000)
  mpred <- t(apply(predict(mnet, newdata=x.test, type="probs"),1,cumsum))
  for(i in 1:nit) {
    cat(".")
    U <- runif(nrow(mpred))
    mi.data[,i] <- colnames(mpred)[ apply(U < mpred, 1, function(x) which(x)[1]) ]
  }
  rownames(mpred) <- rownames(mi.data) <- wdata[,"sample"]
  return(list(probs=mpred,mi.data=mi.data))
}


## summarize errors
serrors <- function(errors, stat="Mean") {
  sapply(errors, function(e) {
    sapply(e[grep("errors",names(e))], function(x) summary(rowMeans(x)))[stat,] })
}

## tabulate example errors - where do we go wrong?
terrors <- function(data,errors,mi=1,samples=NULL) {
  if(!is.null(samples))  {
    pred <- errors$mi.data[samples,mi]
    true <- data[samples,"KIR"]
  } else {
    pred <- errors$mi.data[,mi]
    true <- data[,"KIR"]
  }
  tt <- table(true=true,pred=pred)
  print(1-sum(diag(tt))/sum(tt))
  invisible(tt)
}

## nearest neighbour
knn.impute <- function(wdata,trim.pattern="seq.rs592645") {
  kcols <- grep("KIR",colnames(wdata))
  tcols <- grep(trim.pattern,colnames(wdata))
  message("using ",length(tcols)," variables to predict ",kcols,"\n",paste(colnames(wdata)[tcols],collapse=" "))
  use <- complete.cases(wdata[,tcols])

  ## no missing allowed
  print(table(use))
  ## use
  ## FALSE  TRUE 
  ##     2 12142 
  wdata <- wdata[use,]

  train <- !is.na(wdata[,"KIR"])
  x <- wdata[,tcols]
  x.train <- x[train,]
  x.test <- x ## all samples
  y <- wdata[train,"KIR"]
#  print(summary(x.train))
  m<-knn(train=x.train, test=x.test, cl=y, k=3, prob=TRUE)
  data.frame(sample=rownames(x.test),
             knn.call=m,
             knn.prob=attributes(m)$prob)
}

## nearest neighbour
knn.loo <- function(wdata, k=3) {
  kcols <- grep("KIR",colnames(wdata))

  ## no missing allowed
  use <- complete.cases(wdata[,-kcols])
  print(table(use))
  ## use
  ## FALSE  TRUE 
  ##     2 12142 
  wdata <- wdata[use,]

  train <- !is.na(wdata[,"KIR"])
  x <- wdata[train,-c(1,kcols)]
  y <- wdata[train,"KIR"]
  y.nn <- character(length(y))
  prob <- numeric(length(y))
  for(i in 1:nrow(x)) {
    x.train <- x[-i,]
    x.test <- x[i,,drop=FALSE]
    y.train <- y[-i]
    m<-knn(train=x.train, test=x.test, cl=y.train, k=k, prob=TRUE)
    y.nn[i] <- as.character(m)
    prob[i] <- attributes(m)$prob
  }
  data.frame(sample=rownames(x.test),
             y=y,
             knn.call=y.nn,
             knn.prob=prob)
}
