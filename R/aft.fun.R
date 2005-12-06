aft.fun<-function(x,y,delta,randomseed=10,weight="logrank",nstep=3, mcsize=100)
  {
# Required Input: x,y and delta
#                 x: covariates (in matrix or data.frame form)
#                 y: response (y should be tranformed to logarithm scale 
#                         in advance, no transformation is done in program) 
#                 delta: indicator of censoring, 1 means uncensored, 
#                        0 means censored (this program can not handle 
#                        1 and 2 censoring indicator)
# Optional Input:  randomseed, nstep, mcsize.
#                 randomseed: random seed used in random number generation
#                 nstep: number iteration to compute logrank estimate
#                 mcsize: Monte Carlo run numbers used in computiong the 
#                         covariance matrix of coefficient estimates.
#  for factor covariates, the program uses the dummy variables
# output: beta, betaw, covw
#    beta: first colume of it is Gehan estimate, and rest columns are
#          logrank like estimates
#    betaw: estimates from resampling 
#    covw: covariance of beta, the 1st matrix in the array is the
#          covariance matrix of the Gehan estimate, and rest of them
#          are corresponding covriance matrix of the logran-like estimates
# Reference: Jin, Ying, Wei (2002?) Biometrika.

if(any((delta!=0)&(delta!=1)))
     stop("delta must be 0(right-censored) or 1(uncensored)")
 
  set.seed(randomseed)   #set seed for random number generator
  ynew<-1000*(length(y))^2
  data1<-data.frame(y,x)
  options(contrasts=c("contr.treatment","contr.poly"))
  tempfit<-lm(y~.,x=TRUE,y=TRUE,data=data1)
  x<-as.matrix(tempfit$x[,-1])
  xn<-dimnames(x)[[2]]
  y<-tempfit$y
  dimnum<-dim(x) # dimnum[1]=number of observation
                    # dimnum[2]=number of covariates, excluding intercept
  n1<-dimnum[1]
  n2<-dimnum[2]
  betagw<-betagc<-matrix(0, nrow=mcsize,ncol=n2)
  if(weight=="logrank"){
    nst<-nstep
  }
  else {nst<-0}
 
  betagm<-betalw<-array(0, dim=c(n2,mcsize))
  betagc<-betalc<-array(0, dim=c(n2,,mcsize))
  covw<-array(0, dim=c(n2,n2, nst+1))
  residuals<-matrix(0, nrow=n1,ncol=nst+1)
  
  yy0<-rep(y,rep(n1,n1))
  delta1<-rep(delta,rep(n1,n1))
  yy1<-rep(y,n1)
  yy<-delta1*(yy0-yy1)
  xx0<-matrix(rep(as.vector(x),rep(n1,n1*n2)),nrow=n1*n1)
  xx1<-t(matrix(rep(as.vector(t(x)),n1),nrow=n2))
  xx<-xx0-xx1   #difference of xx0-xx1
  xxdif<-xx*delta1 #difference of xx0-xx1 multiplying censoring indicator
  xnew<-apply(xxdif,2,sum)
  xnew<-rbind(xxdif,-xnew)
  yynew<-c(yy,ynew)
  dimnames(xnew)<-list(NULL,xn)

  fit<-rq.fit(xnew,yynew, tau=0.5)
  betag<-fit$coef
  residn<-fit$resid
  residn<-(!(residn>0))
  residn<-residn[-(length(residn))]
  betal<-betag
  if (weight=="logrank"){
     for (i in 1:nstep){
      fitted<-x%*%betal
      eb<-y-fitted
      ss0b<-(n1+1-rank(eb))/n1
      ss0b1<-rep(ss0b,rep(n1,n1))
      xxdifl<-xxdif/ss0b1
      xnewl<-apply(xxdifl,2,sum)
      xnewl<-rbind(xxdifl,-xnewl)
      yyl<-c(yy/ss0b1,ynew)
      fitl<-rq.fit(xnewl,yyl, tau=0.5)
      betal<-fitl$coef
    }
  }
  
  zi<-matrix(rexp(n1*mcsize),nrow=mcsize,ncol=n1)
  for(i in 1:mcsize) {
      zzi<-rep(as.vector(zi[i,]),rep(n1,n1))
      xm<-xxdif*zzi
      xmnew<-apply(xm,2,sum)
      xmnew<-rbind(xm,-xmnew)
      ymnew<-c(yy*zzi,ynew)
      fitm<-rq.fit(xmnew,ymnew,tau=0.5)
      betagm[,i]<-fitm$coef
      betalw[,i]<-betagm[,i]
      betagc[,i]<-betagm[,i]-betag
      if (weight=="logrank"){
        for (j in 1:nstep){
          fitted<-x%*%betalw[,i]
          eb<-y-fitted
          ss0b<-(n1+1-rank(eb))/n1
          ss0b1<-rep(ss0b,rep(n1,n1))
          xxdifl<-xm/ss0b1
          xnewl<-apply(xxdifl,2,sum)
          xnewl<-rbind(xxdifl,-xnewl)
          yyl<-c(zzi*yy/ss0b1,ynew)
          fitml<-rq.fit(xnewl,yyl,tau=0.5)
          betalw[,i]<-fitml$coef
        }  
      betalc[,i]<-betalw[,i]-betal
      }
    }
   predmatrix<-x-t(matrix(rep(apply(x,2,mean),n1),ncol=n1))
   covw[,,1]<-(betagc)%*%t(betagc)/mcsize
   residuals[,1]<-y-predmatrix%*%as.matrix(betag)
   covw[,,2]<-(betalc)%*%t(betalc)/mcsize
   residuals[,2]<-y-predmatrix%*%as.matrix(betal)

  object<-list(beta=rbind(betag,betal),betacov=covw,residuals=residuals,
               betagm=betagm,betalw=betalw,
               mcsize=mcsize,message=fit$message,warning=fit$warning,
               weight=weight)
  class(object)<-"AFT"
  object
}

##!date
##
##fit1aft<-aft.fun(x=cbind(pbc02.dat[,3],log(pbc02.dat[,4]),log(pbc02.dat[,5]),
## pbc02.dat[,6],log(pbc02.dat[,7])),y=log(pbc02.dat[,1]),delta=pbc02.dat[,2],
## nstep=5, mcsize=100)
##fit1aft1<-aft.fun(x=cbind(pbc02.dat[,3],log(pbc02.dat[,4]),log(pbc02.dat[,5]),
## pbc02.dat[,6],log(pbc02.dat[,7])),y=log(pbc02.dat[,1]),delta=pbc02.dat[,2],
### nstep=5, mcsize=100)
###fit1aft1$beta
### (From Jin Zhezhen, 2004, May for Splus)
###  R version ported by Mai Zhou
