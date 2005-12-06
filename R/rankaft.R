rankaft <- function(x,y,delta) {
#  Input:   x: covariates (in matrix or data.frame form)
#           y: response (y should be tranformed to logarithm scale 
#              in advance, no transformation is done in this program) 
#       delta: indicator of censoring, 1 means uncensored, 
#              0 means censored (this program can not handle 
#                      1 and 2 censoring indicator)
#  output:    (betag and betal).
#        beta: first colume of it is Gehan estimator, and 
#              rest columns are logrank like estimator (via iteration).
#  Notice: this function is very memory hungary.
#######################################################################
##       This function is an R port of the Splus function aft.fun( ) 
##       from Jin ZheZhen. The use of l1fit( ) is replaced by rq( ) from 
##       R quantreg library.   And the resampling part is cut.
##       In the computations, I always get same result as Splus function
##       aft.fun( ) for the Gehan estimator but not
##       always for the log-rank estimator.
##########################################################################
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
  
  residuals<-matrix(0, nrow=n1,ncol=3+1)
  
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

  fit<- rq(yynew ~ xnew -1, method="fn")
  betag<-fit$coef
  residn<-fit$resid
  residn<-(!(residn>0))
  residn<-residn[-(length(residn))]
  betal<-betag
     for (i in 1:3){ 
      fitted<-x%*%betal
      eb<-y-fitted
      ss0b<-(n1+1-rank(eb))/n1
      ss0b1<-rep(ss0b,rep(n1,n1))
      xxdifl<-xxdif/ss0b1
      xnewl<-apply(xxdifl,2,sum)
      xnewl<-rbind(xxdifl,-xnewl)
      yyl<-c(yy/ss0b1,ynew)
      fitl<- rq(yyl ~ xnewl - 1, method="fn")
      betal<-fitl$coef
    }
  
   predmatrix<-x-t(matrix(rep(apply(x,2,mean),n1),ncol=n1))
   residuals[,1]<-y-predmatrix%*%as.matrix(betag)
   residuals[,2]<-y-predmatrix%*%as.matrix(betal)
         
  object<-list(beta=rbind(betag,betal),residuals=residuals, 
               message=fit$message)
  class(object)<-"AFT"
  object
}
