library(spade)
source('~nikolas/bin/FCS/fcs.R')

# The numbering of the nodes in the binary tree is breadthwise.
# There is 2**branch.length nodes in the tree.
# The data.frame splits contains the splits in the non-leaf nodes.

# median absolute difference
mad <- function(x) median(abs(x-median(x)))

## leaf.min=1000
## branch.length=2
## splitting stops when either of these conditions is false
create.splits <- function(X,fun=var,branch='1',leaf.min=100, branch.length=2, splits=c()) {
    #needs to be a data.frame rather than a matrix because rownames need to be defined
    X <- as.data.frame(X)
    branch.index <- strtoi(branch,base=2)
    #leaf node
    if (nrow(X)<leaf.min | nchar(branch)>branch.length) {
        X <- cbind(X,cluster=branch.index)
        X <- X[order(as.numeric(rownames(X))),]
        return(list(splits=NULL,X=X))
    }
    param <- as.character(colnames(X)[[which.max(apply(X,2,fun))]])
    med <- median(X[,param])
    splits <- data.frame(index=branch.index,param=param,split=med,stringsAsFactors=FALSE)
    X1 <- X[which(X[,param] < med),]
    X2 <- X[which(X[,param] >= med),]
    r1 <- create.splits(X1,fun,paste(branch,'0',sep=''),leaf.min,branch.length,splits)
    r2 <- create.splits(X2,fun,paste(branch,'1',sep=''),leaf.min,branch.length,splits)
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

##
compute.mst <- function(X) {
    ### make MST from cluster centers
    adjacency  <- as.matrix(dist(X, method='manhattan'))
    full_graph <- graph.adjacency(adjacency,mode="undirected",weighted=TRUE)
    mst_graph  <- minimum.spanning.tree(full_graph)
    return(mst_graph)
}


##
compute.mst.layout <- function(X) {
    mst_graph <- compute.mst(X)
    layout.table <- SPADE.layout.arch(mst_graph)
    return(layout.table)
}

##
plot.mst <- function(layout.table,color,vsize=1) {
    plot(layout.table, cex=vsize, col=color, pch=20, yaxt='n', xaxt='n', ann=FALSE,frame.plot=FALSE)
}


