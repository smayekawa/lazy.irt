#' Numerical Derivatives of ICRF w.r.t. theta
#'
#' @param param Item Parameter Data Frame
#' @param theta Discrete theta values
#' @param log = 1 to calculate log derivatives.
#' @param second = 1 to calculate second derivatives
#' @param eps eps for JacobianMat
#'
#'
#' @details
#' This function uses lazy.mat::JacobianMat to calculate the first derivative.
#' Except for the second derivatives, the analytic functions are available
#' in \code{dirf}.
#'
#' @return
#' dICRF  length(theta) x sum(ncat) matrix
#'
#'
#' @examples
#' res1 <- dicrf_num( paramS1, -3:3 )
#' res2 <- dirf( paramS1, -3:3, print=0, plot=0 )$dICRF
#' Print(apply(abs(res1-res2),2,max))
#'
#' # second derivative of P
#' res3 <- dicrf_num( paramS1, -3:3, second=1 )
#'
#'
#' \dontrun{
#'
#' thmin <- -4
#' thmax <- 4
#' npoints <- 121
#' theta <- seq(thmin,thmax,length=npoints)
#'
#' param <- paramA1[1:7,]
#'
#'
#' log <- 1
#' res1 <- dicrf_num( param, theta, log=log, second=1 )
#'
#' # all
#' matplot(theta,res1, type="l", main="2nd deriv. of logP")
#'
#' # 0-th category: all negative
#' fromP <- get_range(param$ncat )[,1]
#' matplot(theta,res1[,fromP], type="l", main="2nd deriv. of logP: 0-th cat")
#'
#' # all negative except for the 1st category of 3PLM
#' matplot(theta,res1[,-c(4,8)], type="l"
#'         , main="2nd deriv. of logP: all but 1st of 3PLM")
#'
#' # 1st category of 3PLM
#' matplot(theta,res1[,c(4,8)], type="l", main="2nd deriv. of logP: 1st of 3PLM")
#'
#' }
#'
#' @export
#'
#'

dicrf_num <- function( param, theta, log=0, second=0, eps=1e-6 ){
 # numerical version of dirf
 # Shin-ichi Mayekawa
 # 20180111
 # log: 20180128
 # second: 20180129
 # roxgen2 example: 20180203, 20221206,08
 #

 icrf <- function( x, param, log ){
  res=irf( param, x, plot=0, print=0 )$ICRF
  if( log ) res=log(res)
  return( res )
 } # end of icrf


 dicrf <- function( x, param, log ){
  res=dirf( param, x, log=log, plot=0, print=0 )$dICRF
  return( res )
 } # end of dicrf


 nitem=nrow(param)
 ncolP=sum(param$ncat)
 npoint=length(theta)
 res=matrix(0,npoint,ncolP)
 rownames(res)=format(theta,digits=3)

 if( second ){
  for( i in 1:npoint ){
   # res[i,]=jacobian( icrf, theta[i], param=param )
   res[i,]=JacobianMat( theta[i], dicrf, ..eps..=eps, param=param, log=log )
  }
 }
 else{
  for( i in 1:npoint ){
   # res[i,]=jacobian( icrf, theta[i], param=param )
   res[i,]=JacobianMat( theta[i], icrf, ..eps..=eps, param=param, log=log )
  }
 }
 colnames(res)=colnames(icrf(theta[1],param,log))
 return( res )

} # end of dicrf_num







# second derivative of log P
param=paramS1
pres=fitP2G_ls( param[3,], print=0,plot=0 )
param[4,]=pres$paramP
param[4,"name"]="Q3p"
theta=seq(-4,4,length=31)

log=1
res1=dicrf_num( param, theta, log=log, second=1 )

# all
matplot(theta,res1, type="l", main="2nd deriv. of logP")

# 0-th category
fromP=get_range(param$ncat )[,1]
matplot(theta,res1[,fromP], type="l", main="2nd deriv. of logP: 0-th cat")

# all negagive exept for the 1st category of 3PLM
matplot(theta,res1[,-4], type="l"
        , main="2nd deriv. of logP: all but cat 1 of 3PLM")

# 1st category of 3PLM
matplot(theta,res1[,4], type="l", main="2nd deriv. of logP: cat 1 of 3PLM")




comments(
 '
# library(numDeriv)



param=paramS1
pres=fitP2G_ls( param[3,], print=0,plot=0 )
param[4,]=pres$paramP
param[4,"name"]="Q3p"
theta=seq(-4,4,length=31)

log=1
res1=dicrf_num( param, theta, log=log )

res2=dirf( param, theta, plot=2, print=0, log=log )

dd=res2$dICRF
Print(max(abs(res1-dd)))








param=paramS2
theta=seq(-4,4,length=31)

res1=dicrf_num( param, theta, log=1 )

res2=dirf( param, theta, plot=2, print=0, log=1 )

dd=res2$dICRF
Print(max(abs(res1-dd)))





icrf( 0, param )

icrf( -1:1, param )

dirf( param, -1:1, print=0, plot=0 )$dICRF

JacobianMat( -1:1, icrf, param=param )

jacobian( icrf, 0, param=param )

jacobian( icrf, -1:1, param=param )







 param=paramS1
 pres=fitP2G_ls( param[3,], plot=1 )
 param[4,]=pres$paramP
 param[4,"name"]="Q31p"

 theta=seq(-4,4,length=121)
 res0=irf( paramS1, theta, plot=1, smallP=0 )
 lP=log(res0$ICRF)
 matplot(theta,lP[,res0$fromP[1]:res0$toP[1]],type="l")
 matplot(theta,lP[,res0$fromP[2]:res0$toP[2]],type="l")
 matplot(theta,lP[,res0$fromP[3]:res0$toP[3]],type="l")
 matplot(theta,lP[,res0$fromP[4]:res0$toP[4]],type="l")


 res2=dirf( param, theta, plot=1, print=0, log=1 )





')
