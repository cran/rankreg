\name{RankRegV}
\alias{RankRegV}
\title{Test and Variance estimator of rank estimating function 
in censored AFT model}
\usage{
RankRegV(y, d, x, beta, betanull = beta, type="Gehan")  
}
\arguments{
    \item{y}{a vector containing the censored responses in the AFT model.}
    \item{d}{a vector of 1's and 0's. censoring indicator. 1(uncensor), 
             0(censored). Both y and d should be of length n.}
    \item{x}{the design matrix, of size n by q. 
                        Should not have a column of 1s.}
    \item{beta}{a vector of dimension q. Usually it should be the solution
               of the estimation equation, from \code{rankreg()}. }
    \item{betanull}{the null value (under H0) of beta to be tested.}
    \item{type}{either equal to "Gehan" (default) or "Logrank".}
}
\description{
Compute the score type test,  with 
the chi square value for testing Ho: betanull.

It also resturns the variance-covariance matrix of the 
rank estimating function (score function)
at beta for the censored AFT model.
}
\details{
It returns the chi square value and P-value which are based on
\code{EF(betanull)*Vtest(beta)^(-1)* EF(betanull)} . 

The input \code{beta} usually should be the solution
of the estimating function. 
This solution can be obtained by using \code{rankaft}.
But other values of beta may also make sense sometime.

When \code{betanull} equals \code{beta} from \code{rankaft},
then the chi square
statistic should be zero (since \code{EF(betanull)=0})
and the P-value should be 1.
}
\value{
A list with \code{VEF}, which is the variance estimator of the
estimating function at \code{beta}; 

\code{chisqvalue} which
is the quadratic form \code{EF(betanull)*VEF^(-1)*EF(betanull)};

and \code{Pval} that is the P-value from the chi square quantile.
}
\references{
Kalbfleisch, J. and Prentice, R. (2002)
{\em The Statistical Analysis of Failure Time Data.} 2nd Ed.
Wiley, New York.  (Chapter 7) 
}
\author{ Mai Zhou. }
\examples{
data(myeloma)
RankRegV(y=myeloma[,1],d=myeloma[,2],x=cbind(myeloma[,3],myeloma[,4]),
         beta=c(-15.01117, 1.317596), betanull=c(-15.01117, 1.317596))
### you should get a chisquare = 0.001124574 , Pval = 0.9994379
GP1 <- c(143,164,188,188,190,192,206,209,213,216,220,227,230,234,246,265,304,
       216,244)
GP2 <- c(142,156,163,198,205,232,232,233,233,233,233,239,240,261,280,280,296,
       296,323,204,344)
status1 <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0)
status2 <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0)
RankRegV(y=c(GP1,GP2),d=c(status1,status2),x=c(rep(0,19),rep(1,21)),
         beta=0,betanull=0)
}
\keyword{nonparametric}
\keyword{htest}
