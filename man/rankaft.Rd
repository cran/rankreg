\name{rankaft}
\alias{rankaft}
\title{Compute Rank estimator in right censored data AFT model}
\usage{
rankaft(x, y, delta)  
}
\arguments{
    \item{x}{the design matrix, of size n by q.}
    \item{y}{a vector containing the censored responses in the AFT model.}
    \item{delta}{a vector of 1's and 0's. censoring indicator. 1(uncensor), 
             0(censored). Both y and delta should be of length n.}
}
\description{
   Compute the Gehan and Logrank type rank regression estimators in 
   the censored AFT model, using linear programming. 
   This function is similar to \code{aft.fun()} except we strip away
   the re-sampling part in order to speedup things.
}
\details{
This program is memory hungry. Caution: at least
1G of RAM needed for sample size 1000; 
at least 512MB RAM for sample size 400.

We cut the re-sampling part (from \code{ aft.fun ()} )
to save computing time/memory, and concentrate on the bottleneck.

For statistical infernce, there are three options: (1)
re-sampling method to estimate the var-cov matrix
(available in aft.fun)
(2) score type test available from function \code{ RankRegV() }
and (3) by empirical likelihood (see the reference).
}
\value{
   A list with \code{beta} which is the Gehan (betag)
   and Logrank type (betal) estimate 
   rbinded together; and \code{residuals}. 
}
\references{
Kalbfleisch, J. and Prentice, R. (2002)
{\em The Statistical Analysis of Failure Time Data}. 2nd Ed.
Wiley, New York.  (Chapter 7) 

Jin, Z., Lin, D.Y.,
Wei, L. J. and Ying, Z. (2003).
Rank-based inference for the accelerated failure time model.
{\em Biometrika}, {\bf 90}, 341-353.

Zhou, M. (2005).
Empirical likelihood analysis of the rank estimator for the 
censored AFT model.
{\em Biometrika}, {\bf 92}, 492-498.
}
\examples{
data(myeloma)
rankaft(x=cbind(myeloma[,3],myeloma[,4]),y=myeloma[,1],delta=myeloma[,2])
}
\author{ Original Splus code by Jin Z. Adapted to R by Mai Zhou.}
\keyword{nonparametric}
\keyword{htest}
