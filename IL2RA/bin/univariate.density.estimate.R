


np.density <- function(mm) {
x <- np.mm$data
u <- x 
bw <- mm$bandwidth

Km <- exp(outer(x/bw, u/bw, function(a, b) -(a - b)^2/2))

normw <- matrix(mm$posteriors[,1]/sum(mm$posteriors[,1]), nrow = 1)
as.vector(normw %*% Km)/(bw * sqrt(2 * pi)) -> result

normw <- matrix(mm$posteriors[,2]/sum(mm$posteriors[,2]), nrow = 1)
as.vector(normw %*% Km)/(bw * sqrt(2 * pi)) -> result2

return(mm$lambdahat[1]*result+mm$lambdahat[2]*result2)
}



pdf('Rplots.pdf')
plot(density(x,bw=.01), col="gray")
lines(density(x,bw=bw.nrd0(x)), col="darkgreen")
lines(density(x,bw=bw.nrd(x)), col="green")
lines(density(x,bw=bw.ucv(x)), col="blue")
lines(density(x,bw=bw.bcv(x)), col="purple")
lines(density(x,bw=bw.SJ(x)), col="darkred")
points(x, np.density(np.mm), col="red", pch="." )
points(x, np.density(sp.mm), col="red", pch="x" )
sort(x)->x.sort
lines(x.sort, mm$lambda[1]*dnorm(x.sort, mm$mu[1], mm$sigma[1])+mm$lambda[2]*dnorm(x.sort, mm$mu[2], mm$sigma[2]))
dev.off()


> wkde
function (x, u = x, w = rep(1, length(x)), bw = bw.nrd0(as.vector(x)), 
    sym = FALSE) 
{
    if (sym) {
        return((wkde(x, u, w, bw) + wkde(x, -u, w, bw))/2)
    }
    Km <- exp(outer(x/bw, u/bw, function(a, b) -(a - b)^2/2))
    normw <- matrix(w/sum(w), nrow = 1)
    as.vector(normw %*% Km)/(bw * sqrt(2 * pi))
}
<environment: namespace:mixtools>
