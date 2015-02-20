
create.graph <- function(cluster.stats, channels) {
    channels <- toupper(channels)
    graph <- list()
        for (dose in doses) {
            graph1 <- list()
            graph2 <- list()
            ### try out various network properties
            for (individual in rep.individuals) {
                l <- lapply(cluster.stats[grep(paste(individual, dose, sep='.'), names(cluster.stats))], exprs)
                l <- lapply(l, function(x) { colnames(x) <- toupper(colnames(x)); return(x) })
                #day 1
                g <- graph.adjacency(as.matrix(dist(l[[1]][,channels])), mode='undirected', weighted=TRUE)
                g <- set.vertex.attribute(g, name='pct', value=l[[1]][,'PCT'])
                g2 <- graph.adjacency(as.matrix(dist(l[[2]][,channels])), mode='undirected', weighted=TRUE)
                #day 2
                g2 <- set.vertex.attribute(g2, name='pct', value=l[[2]][,'PCT'])
                graph1[[individual]] <- g
                graph2[[individual]] <- g2
            }
            graph[[dose]] <- list(graph1, graph2)
        }
    return(graph)
}




graph.f.density <- function(graph=graph.core, dir='~nikolas/IL2/Plots/graph-core-distances/', f=function(x) E(x)$weight) {
  dir.create(dir, recursive=T, showWarnings=F)
  for (individual in rep.individuals) {
    pdf(file.path(dir, sprintf('%s.pdf',individual)))
    day1.resting <- f(graph[['0U']][[1]][[individual]])
    day2.resting <- f(graph[['0U']][[2]][[individual]])
    day1.stimulated <- f(graph[['1000U']][[1]][[individual]])
    plot(density(day1.resting), main=individual)
    #dose effect
    #lines(density(E(graph[['0U']][[1]][[individual]])$weight), col='black', lty=2)
    abline(v=median(day1.resting), col='black')
    lines(density(day1.stimulated), col='red')
    abline(v=median(day1.stimulated), col='red')
    #day effect resting
    lines(density(day2.resting), col='black', lty=2)
    abline(v=median(day2.resting), col='black', lty=2)
    dev.off()
  }
}


degree.norm <- function(graph, ...) {
  deg <- degree(g, ...)
  n <- vcount(graph)
  deg/(n-1)
}


betweenness.norm <- function(graph, ...) {
  bet <- betweenness(graph, ...)
  n <- vcount(graph)
  (2 * bet) / (n*n - 3*n + 2)
}



degree.centralisation <- function(graph, ...) {
  deg <- degree(g, ...)
  n <- vcount(graph)
  (max(deg)*n - sum(deg)) / (n*n - 3*n +2)
}


betweenness.centralisation <- function(graph, ...) {
  bet <- betweenness(graph, ...)
  n <- vcount(graph)
  2 * (max(bet)*n - sum(bet)) / (n*n*n - 4*n*n + 5*n - 2)
}


closeness.centralisation <- function(graph, ...) {
  if (!is.connected(graph)) warning("This graph is disconnected; calculation\
  of closeness centrality on a disconnected graph is only possible because\
  igraph assumes a non-infinite distance (that is equal to the number of\
  vertices in the network) between two disconnected vertices.")
  clo <- closeness(g, ...)
  n <- vcount(graph)
  (max(clo)*n - sum(clo))*(2*n - 3) / (n*n - 3*n + 2)
}


graph.summary <- function(graph, f, name, file.name, dir='~nikolas/IL2/Plots/') {
    dir.create(dir, recursive=T, showWarnings=F)
    pdf(file.path(dir, sprintf('%s.pdf',file.name)), width=10, height=5)
    #par(oma=c(0,0,2,0))
    par(mfrow=c(1,2))
    x <- sapply(graph[['0U']][[1]], f)
    y <- sapply(graph[['1000U']][[1]], f)
    plot(x, y, xlab='Dose 0U', ylab='Dose 1000U', pch=letters[1:10], xlim=range(c(x,y)), ylim=range(c(x,y)))
    abline(b=1,a=0)
    r <- round(cor(x,y)**2, digits=3)
    legend('bottomright', legend=as.expression(bquote(r^2 == .(r))), text.col='black')
    x <- sapply(graph[['0U']][[1]], f)
    y <- sapply(graph[['0U']][[2]], f)
    plot(x, y, xlab='Day 1', ylab='Day 2', pch=letters[1:10], xlim=range(c(x,y)), ylim=range(c(x,y)))
    r <- round(cor(x,y)**2, digits=3)
    legend('bottomright', legend=as.expression(bquote(r^2 == .(r))), text.col='black')
    abline(b=1,a=0)
    #title(name, outer=TRUE)
    dev.off()
}


plot.agreement <- function(dose0, dose1, day0, day1, file.name, dir='~nikolas/IL2/Plots/') {
    dir.create(dir, recursive=T, showWarnings=F)
    pdf(file.path(dir, sprintf('%s.pdf',file.name)), width=10, height=5)
    #par(oma=c(0,0,2,0))
    par(mfrow=c(1,2))
    x <- dose0
    y <- dose1
    plot(x, y, xlab='Dose 0U', ylab='Dose 1000U', pch=letters[1:10], xlim=range(c(x,y)), ylim=range(c(x,y)))
    abline(b=1,a=0)
    r <- round(cor(x,y)**2, digits=3)
    legend('bottomright', legend=as.expression(bquote(r^2 == .(r))), text.col='black')
    x <- day0
    y <- day1
    plot(x, y, xlab='Day 1', ylab='Day 2', pch=letters[1:10], xlim=range(c(x,y)), ylim=range(c(x,y)))
    r <- round(cor(x,y)**2, digits=3)
    legend('bottomright', legend=as.expression(bquote(r^2 == .(r))), text.col='black')
    abline(b=1,a=0)
    #title(name, outer=TRUE)
    dev.off()
}


plot.doses.per.day <- function(d, file.name, dir) {
    dir.create(dir, recursive=T, showWarnings=F)
    pdf(file.path(dir, sprintf('%s.pdf',file.name)), width=10, height=5)
    par(mfrow=c(1,2), las=1)
    day1 <- d$day1
    day2 <- d$day2
    ylim <- range(c(day1,day2))
    plot(NULL, axes=FALSE, xlim=c(1,4), ylim=ylim, xlab='', ylab='', main='Day 1')
    axis(1, 1:4, labels=doses)
    axis(2)
    apply(day1, 1, lines)
    plot(NULL, axes=FALSE, xlim=c(1,4), ylim=ylim, xlab='', ylab='', main='Day 2')
    axis(1, 1:4, labels=doses)
    axis(2)
    apply(day2, 1, lines) 
    dev.off()
}




plot.doses.per.day.best.fit <- function(d, file.name, dir, labels=c('0U', '0.1U', '10U', '1000U')) {
    dir.create(dir, recursive=T, showWarnings=F)
    pdf(file.path(dir, sprintf('%s.pdf',file.name)), width=10, height=5)
    par(mfrow=c(1,2), las=1)
    day1 <- d$day1
    day2 <- d$day2
    ylim <- range(c(day1,day2))
    plot(NULL, axes=FALSE, xlim=c(1,4), ylim=ylim, xlab='', ylab='', main='Day 1')
    axis(1, 1:4, labels=labels)
    axis(2)
    apply(day1, 1, points)
    m.day1 <- glm(value ~ Var2, data=melt(day1))
    print(summary(m.day1))
    abline(coef(m.day1))
    #apply(day1, 1, lines) 
    #apply(day1, 1, function(x) abline(coef(line(1:4, x))))
    plot(NULL, axes=FALSE, xlim=c(1,4), ylim=ylim, xlab='', ylab='', main='Day 2')
    axis(1, 1:4, labels=labels)
    axis(2)
    apply(day2, 1, points) 
    #apply(day2, 1, lines) 
    #apply(day2, 1, function(x) abline(coef(line(1:4, x))))
    m.day2 <- glm(value ~ Var2, data=melt(day2))
    print(summary(m.day2))
    abline(coef(m.day2))
    dev.off()
}


