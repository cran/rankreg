\name{aft.fun}
\alias{aft.fun}
\title{Rank estimator in censored AFT model}
\usage{
aft.fun(x, y, delta, randomseed=10, weight="logrank", nstep=3, mcsize=100) 
}
\arguments{
    \item{x}{the design matrix, of size n by q.}
    \item{y}{a vector containing the censored responses in the AFT model.}
    \item{delta}{a vector of 1's and 0's. censoring indicator. 1(uncensor), 
             0(censored). Both y and d should be of length n.}
    \item{randomseed} {the random seed used to generate re-samples.}
    \item{weight} {}
    \item{nstep}{an integer. The number of iterations used to compute the
        logrank type estimator starting from the Gehan estimator.  }
    \item{mcsize}{number of resamples used to compute the variance
                  estimator. }
}
\description{
   Use linear programming to solve the Gehan rank estimation 
   equation for censored AFT model. Use iterated Gehan type solutions
   to solve the Logrank estimation equation.
   Finally, it computes the variance-covariance estimator for both 
   rank regression estimators in the censored AFT model by re-sampling.
}
\details{
For data sets with more than 400 observations, this function is slow.
The reason is that it needs to solve linear programming problems
of size n square. So 400 becomes 160000.
}
\value{
   A list with the following components.
   \code{beta}: first column is the Gehan estimator, the rest are
    logrank type estimators.
   \code{betaw}: the estimates from re-sampling.
   \code{covw}: variance-covariance estimator of \code{beta} from resampling.
}
\references{
Jin, Z., Lin, D.Y., Wei, L. J. and Ying, Z. (2003).
Rank-based inference for the accelerated failure time model.
{\em Biometrika}, {\bf 90}, 341-353.

Kalbfleisch, J. and Prentice, R. (2002)
The Statistical Analysis of Failure Time Data. 2nd Ed.
Wiley, New York.  (In particular Chapter 7)
}
\author{ Original Splus code by Z. Jin. Adapted to R by Mai Zhou. }
\keyword{nonparametric}
\keyword{htest}
