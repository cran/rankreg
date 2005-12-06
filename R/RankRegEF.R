RankRegEF <- function(y, d, x, beta, type="Gehan") {
n <- length(y)                   ## dimension of x must be n x q.
x <- as.matrix(x)                ## x must NOT including intercept.
xdim <- dim(x)
if( xdim[1] != n ) stop("check dim of x")
if( length(beta) != xdim[2] ) stop("check dim of beta and x")
if(any((d!=0)&(d!=1)))
     stop("d must be 0(right-censored) or 1(uncensored)")
if(length(d) != n) stop("check the length of d")

e <- y - as.vector( x %*% beta )
ordere <- order(e, -d)
esort <- e[ordere]
dsort <- d[ordere]
xsort <- as.matrix(x[ordere,])
dsort[length(dsort)] <- 1       #last one as uncensored always

##xbar <- rev(cumsum(rev(xsort)))/(n:1)  #for 1 dim

xbar <- xsort
####for(j in 1:(n-1)) xbar[j,] <- colMeans(xsort[j:n,])
for(j in 1:xdim[2]) xbar[,j] <- rev(cumsum(rev(xsort[,j])))/(n:1)

if(type == "Gehan") {A <- (n:1) * (xsort - xbar)}
 else {if(type == "Logrank") A <- (xsort - xbar)
        else stop("type must be either Gehan or Logrank") }
A <- dsort*A

rootnV <- sqrt(n) * colMeans(A)
return(rootnV)
}
