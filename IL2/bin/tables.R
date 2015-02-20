

#'Alexa Fluor 700-A'    'CD4'
#'Pacific Blue-A' 'CD45RA'
#'APC-A'   'CD25'
#'PE YG-A'  'FoxP3'
#'Alexa Fluor 488-A' 'pSTAT5'




tables <- list()
for ( f in
     c(
        'T1D_TC14_pSTAT5_Treg_I022267C_CB00010K',
        'T1D_pSTAT5_TC6_Treg_CB00055J_CB00055J',
        'T1D_pSTAT5_TC2_Treg_I021196N_CB00086S',
        'T1D_TC17_pSTAT5_Treg_CB00165D_I022596K_CB00165D',
        'T1D_pSTAT5_TC5_Treg_I021505Z_CB00357M',
        'T1D_TC13_pSTAT5_Treg_I022215W_CB00366X',
        'T1D_TC15_pSTAT5_Treg_I022371Q_CB00370B',
        'T1D_TC12_pSTAT5_Treg_I022106C_CB00394C',
        'T1D_pSTAT5_TC3_Treg_I021336Q_CB00396E',
        'T1D_TC9_pSTAT5_Treg__I021909N_CB00406Q',
        'T1D_TC9_rep_pSTAT5_Treg_CB00406Q_I023293S',
        'T1D_TC9_rep_pSTAT5_Treg_new_pro_CB00406Q_I023293S',
        'T1D_TC10_pSTAT5_Treg_I021958R_CB00782Z',
        'T1D_pSTAT5_TC1_Treg_I021129Q_CB01145T',
        'T1D_TC7_pSTAT5_Treg_I021697H_CB01179F',
        'T1D_TC14_pSTAT5_Treg_I022288A_CB01338D',
        'T1D_TC16_pSTAT5_Treg_I022420T_CB01340F',
        'T1D_TC11_pSTAT5_Treg_I022056Y_CB01422V',
        'T1D_TC7_pSTAT5_Treg_I021755W_CB01438M',
        'T1D_pSTAT5_TC4_Treg_I021420G_CB01477E',
        'T1D_pSTAT5_TC1_Treg_I021141D_CB01480H',
        'T1D_pSTAT5_TC2_Treg_I021212F_CB01481J',
        'T1D_pSTAT5_TC2_Treg_I021228Y_CB01482K',
        'T1D_TC14_pSTAT5_Treg_I022309Y_CB01483L',
        'T1D_pSTAT5_TC3_Treg_I021347C_CB01484M',
        'T1D_pSTAT5_TC4_Treg_I021445J_CB01486P',
        'T1D_TC10_pSTAT5_Treg_I021976L_CB01487Q',
        'T1D_pSTAT5_TC5_Treg_I021532D_CB01488R',
        'T1D_pSTAT5_TC5_Treg_I021586M_CB01489S',
        'T1D_pSTAT5_TC5_Treg_I021559H_CB01490T',
        'T1D_pSTAT5_TC6_Treg_CB01491V_CB01491V',
        'T1D_TC7_pSTAT5_Treg_I021726P_CB01492W',
        'T1D_TC8_pSTAT5_Treg_I021767J_CB01494Y',
        #'T1D_TC8_pSTAT5_Treg_I021767J_CB01494Y_0U',
        'T1D_TC8_rep_pSTAT5_Treg_CB01494Y_I023360Q',
        #'T1D_TC8_rep_pSTAT5_Treg_CB01494Y_I023360Q_0U',
        'T1D_TC8_pSTAT5_Treg_I021794N_CB01495Z',
        'T1D_TC8_rep_pSTAT5_Treg_CB01495Z_I023343X',
        'T1D_TC15_pSTAT5_Treg_I022395R_CB01496A',
        'T1D_TC9_pSTAT5_Treg__I021927H_CB01498C',
        'T1D_TC9_rep_pSTAT5_Treg_CB01498C_I023276Z',
        'T1D_TC9_rep_pSTAT5_Treg_new_pro_CB01498C_I023276Z',
        'T1D_TC10_pSTAT5_Treg_I022003Q_CB01499D',
        'T1D_TC11_pSTAT5_Treg_I022035A_CB01500E',
        'T1D_TC12_pSTAT5_Treg_I022124X_CB01502G',
        'T1D_TC12_pSTAT5_Treg_I022153D_CB01503H',
        'T1D_TC13_pSTAT5_Treg_I022188R_CB01504J',
        'T1D_TC16_pSTAT5_Treg_I022464R_CB01509P',
        'T1D_TC16_pSTAT5_Treg_I022438N_CB01510Q',
        'T1D_TC17_pSTAT5_Treg_CB01522D_I022623P_CB01522D'
        ) ) {

tables[[f]] <- read.table(sprintf('%s/clusters.table',f), header=T)

}

names(tables) <- gsub('T1D_.*_Treg_', '', names(tables))
names(tables) <- gsub('.*(CB[^_]*).*', '\\1', names(tables))


tables.princomp <- lapply(tables, princomp)

xy <- lapply(tables.princomp, function(x) x$scores[,1:2])

dist.xy <- lapply(xy, dist)



sum(dist(do.call('rbind', tables[1:2]))) - do.call('sum', dist.xy[1:2])


d <- list()
n1 <- 'I021767J_CB01494Y'
for (n in names(tables)) {
    d[[paste(n1, n)]] <- sum(dist(do.call('rbind', dist.xy[c(n1,n)]))) - do.call('sum', dist.xy[c(n1,n)])
}
cbind( names(d)[order(as.numeric(d))] , sort(as.numeric(d)) )



dist.tables <- lapply(tables, dist)
d <- list()
n1 <- 'I021767J_CB01494Y'
for (n in names(tables)) {
    d[[paste(n1, n)]] <- sum(dist(do.call('rbind', tables[c(n1,n)]))) - do.call('sum', tables[c(n1,n)])
}
cbind( names(d)[order(as.numeric(d))] , sort(as.numeric(d)) )



t1 <- tables[[1]]

library(igraph)
g <- graph.adjacency(as.matrix(dist(t1)), mode='undirected', weighted=TRUE)
v.g <- V(g)
g2 <- simplify(g)

set.seed(3952)
layout1 <- layout.fruchterman.reingold(g2)
plot(g2, layout=layout1)

mst <- minimum.spanning.tree(g)
layout.g <- layout.kamada.kawai(mst)
plot(mst, layout=layout.g) 

graphs <- lapply(tables, function(x) graph.adjacency(as.matrix(dist(x)), mode='undirected', weighted=TRUE))

graphs.mst <- lapply(graphs, minimum.spanning.tree)
graphs.diameter <- lapply(graphs, diameter)


library(ggplot2)
#ggplot(df, aes(x=individual, y=diameter, colour=factor(individual)))+geom_bar(stat='identity', position='dodge')+coord_flip()

plot.graphs.fun <- function(graphs, fun, ylab) {
    graphs.fun <- lapply(graphs, fun)
    df <- data.frame(individual=names(graphs.fun), ylab=as.numeric(graphs.fun))
    df <- df[order(df$individual),]
    print(ggplot(df, aes(x=individual, y=ylab))+geom_point(stat='identity')+coord_flip()+theme(legend.position='none')+ylab(ylab))
    return(graphs.fun)
}


tables.all <- do.call('rbind', tables)
tables.all.pc <- princomp(tables.all)
indexes <- cumsum(sapply(tables, function(x) dim(x)[[1]]))
plot(tables.all.pc$scores[,1:2], col='gray', pch=20)
#rownames(
i <- which(names(indexes)=='CB01482K')
points(tables.all.pc$scores[indexes[i]:(indexes[i+1]-1),1:2], col='red', pch=20)
i <- which(names(indexes)=='CB01481J')
points(tables.all.pc$scores[indexes[i]:(indexes[i+1]-1),1:2], col='green', pch=20)
i <- which(names(indexes)=='CB00086S')
points(tables.all.pc$scores[indexes[i]:(indexes[i+1]-1),1:2], col='blue', pch=20)
        
#diameter
graphs.diameter <- plot.graphs.fun(graphs=graphs, fun=igraph::diameter, ylab='diameter')
#mean graph strength
graphs.mean_graph_strength <- plot.graphs.fun(graphs=graphs, fun=function(x) mean(igraph::graph.strength(x)), ylab='mean graph strength')

f <- function(x) mean(graph.knn(x)[[1]])
f <- function(x) determinant(graph.laplacian(x))$modulus
barplot(sapply(graphs, f) , horiz=T, las=2, cex.names=.5, col=as.factor(names(graphs)))

f <- function(x) sum(E(x)$weight)
f <- function(x) mean(E(x)$weight)
f <- function(x) diameter(x)
f <- function(x) transitivity(x, type='undirected')
barplot(sapply(graphs.mst, f), horiz=T, las=2, cex.names=.5, col=as.factor(names(graphs)))

lapply(graphs.mst, cliques)

plot(ecdf(E(graphs.mst[[1]])$weight))
for (i in 1:length(graphs.mst)) lines(ecdf(E(graphs.mst[[i]])$weight))

plot(sort(E(graphs.mst[[1]])$weight), cumsum(E(graphs.mst[[1]])$weight)/sum(E(graphs.mst[[1]])$weight), type='l')
for (i in 1:length(graphs.mst)) lines(sort(E(graphs.mst[[i]])$weight), cumsum(E(graphs.mst[[i]])$weight)/sum(E(graphs.mst[[i]])$weight))
for (i in c("CB00086S", "CB01481J", "CB01482K")) lines(sort(E(graphs.mst[[i]])$weight), cumsum(E(graphs.mst[[i]])$weight)/sum(E(graphs.mst[[i]])$weight), col='red')


plot(sort(E(graphs[[1]])$weight), cumsum(E(graphs[[1]])$weight)/sum(E(graphs[[1]])$weight), type='l')
for (i in 1:length(graphs)) lines(sort(E(graphs[[i]])$weight), cumsum(E(graphs[[i]])$weight)/sum(E(graphs[[i]])$weight))


