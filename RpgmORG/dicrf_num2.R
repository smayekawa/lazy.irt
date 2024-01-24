#' Numerical Derivatives of ICRF w.r.t. theta
#'
#' @param param Item Parameter Data Frame
#' @param theta Discrete theta values
#' @param log = 1 to calculate log derivatives.
#' @param eps eps for JacobianMat
#'
#'
#' @details
#' This function uses JacobianMat to calculate the first derivative.
#'
#' @return
#' dICRF  length(theta) x sum(ncat) matrix
#'
#' @examples
#' res1=dicrf_num( paramS1, -1:1 )
#' res2=dirf( paramS1, -1:1, print=0, plot=0 )$dICRF
#' Print(max(abs(res1-res2)))
#'
#' @export
#'
#'

dicrf_num <- function( param, theta, log=0, eps=1e-6 ){
 # numerical version of dirf
 # Shin-ichi Mayekawa
 # 20180111
 # log: 20180128
 #

 icrf <- function( x, param, log ){
  res=irf( param, x, plot=0, print=0 )$ICRF
  if( log ) res=log(res)
  return( res )
 } # end of icrf

 nitem=nrow(param)
 ncolP=sum(param$ncat)
 npoint=length(theta)
 res=matrix(0,npoint,ncolP)
 rownames(res)=format(theta,digits=3)
 for( i in 1:npoint ){
  # res[i,]=jacobian( icrf, theta[i], param=param )
  res[i,]=JacobianMat( theta[i], icrf, ..eps..=eps, param=param, log=log )
 }
 colnames(res)=colnames(icrf(theta[1],param,log))
 return( res )

} # end of dicrf_num








param=paramS1
pres=fitP2G_ls( param[3,], print=0,plot=0 )
param[4,]=pres$paramP
param[4,"name"]="Q3p"
theta=seq(-4,4,length=31)

res1=dicrf_num( param, theta, log=1 )

res2=dirf( param, theta, plot=2, print=0, log=1 )

dd=res2$dICRF
Print(max(abs(res1-dd)))




comments(
 '
# library(numDeriv)






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
