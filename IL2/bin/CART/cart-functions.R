fun <- var
fun <- function(x) sum((x-mean(x))**2)

# The numbering of the nodes in the binary tree is breadthwise.
# There is 2**branch.length nodes in the tree.
# The data.frame splits contains the splits in the non-leaf nodes.

# median absolute difference
mad <- function(x) median(abs(x-median(x)))

## leaf.min=1000
## branch.length=2
## splitting stops when either of these conditions is false
create.splits <- function(X,markers,fmarker,fun=var,branch='1',leaf.min=100, branch.length=2, splits=c()) {
    #needs to be a data.frame rather than a matrix because rownames need to be defined
    X <- as.data.frame(X)
    branch.index <- strtoi(branch,base=2)
    #leaf node
    if (nrow(X)<leaf.min | nchar(branch)>branch.length) {
        X <- cbind(X,cluster=branch.index)
        X <- X[order(as.numeric(rownames(X))),]
        return(list(splits=NULL,X=X))
    }
    splits <- do.call('rbind', lapply(markers, function(marker) {
        x <- X[,marker]
        sp <- quantile(x,probs=seq(0.1,.9,length=100))
        fun.sum <- function(p) fun(X[which(x<p),fmarker])+fun(X[which(x>p),fmarker])
        sp.hat <- sp[which.min(sapply(sp, fun.sum))]
        return(data.frame(marker=marker,sp=sp.hat,cost=fun.sum(sp.hat),stringsAsFactors=FALSE))
        }
    ))
    s <- splits[which.min(splits$cost),]
    splits <- data.frame(index=branch.index,param=s$marker,split=s$sp,cost=s$cost,n=nrow(X),stringsAsFactors=FALSE)
    X1 <- X[which(X[,s$marker] < splits$split),]
    X2 <- X[which(X[,s$marker] >= splits$split),]
    r1 <- create.splits(X1,markers,fmarker,fun,paste(branch,'0',sep=''),leaf.min,branch.length,splits)
    r2 <- create.splits(X2,markers,fmarker,fun,paste(branch,'1',sep=''),leaf.min,branch.length,splits)
    X <- rbind(r1$X,r2$X)
    X <- X[order(as.numeric(rownames(X))),]
    splits <- rbind(splits,r1$splits,r2$splits)
    return(list(splits=splits,X=X))
}

## use.split=TRUE
## Can either use split points defined in splits or just the parameter and redefine the
## split point.
## Currently this algorithm will not work if you have empty bins.
apply.splits <- function(X, splits, branch='1', use.split=TRUE) {
    #needs to be a data.frame rather than a matrix because rownames need to be defined
    X <- as.data.frame(X)
    branch.index <- strtoi(branch,base=2)
    if (!branch.index %in% splits$index) {
        X <- cbind(X,cluster=branch.index)
        X <- X[order(as.numeric(rownames(X))),]
        return(list(splits=NULL,X=X))
    }
    s <- splits[which(branch.index==splits$index),]
    if (!use.split) s$split <- median(X[,s$param])
    X1 <- X[which(X[,s$param] < s$split),]
    X2 <- X[which(X[,s$param] >= s$split),]
    r1 <- apply.splits(X1,splits,paste(branch,'0',sep=''))
    r2 <- apply.splits(X2,splits,paste(branch,'1',sep=''))
    X <- rbind(r1$X,r2$X)
    X <- X[order(as.numeric(rownames(X))),]
    return(list(splits=splits,X=X))
}

##
cluster.stats <- function(fcs.data,cs) {
    cs2 <- apply.splits(fcs.data, cs$splits)
    K <- cs2$X$cluster
    cluster <- as.data.frame(do.call('rbind',by(fcs.data,K,colMedians)))
    cluster$n <- as.numeric(table(K))
    cluster$prop <- prop.table(as.numeric(table(K)))
    cluster$cluster <- as.numeric(rownames(cluster))
    return(cluster)
} 


