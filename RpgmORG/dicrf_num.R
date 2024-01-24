#' Numerical Derivatives of ICRF w.r.t. theta
#'
#' @param param Item Parameter Data Frame
#' @param theta Discrete theta values
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

dicrf_num <- function( param, theta, eps=1e-6 ){
 # numerical version of dirf
 # Shin-ichi Mayekawa
 # 20180111
 #

 icrf <- function( x, param ){
  res=( irf( param, x, plot=0, print=0 )$ICRF )
  return( res )
 } # end of icrf

 nitem=nrow(param)
 ncolP=sum(param$ncat)
 npoint=length(theta)
 res=matrix(0,npoint,ncolP)
 rownames(res)=format(theta,digits=3)
 for( i in 1:npoint ){
  # res[i,]=jacobian( icrf, theta[i], param=param )
  res[i,]=JacobianMat( theta[i], icrf, ..eps..=eps, param=param )
 }
 colnames(res)=colnames(icrf(theta[1],param))
 return( res )

} # end of dicrf_num




library(numDeriv)

param=paramS2
theta=seq(-4,4,length=31)

res1=dicrf_num( param, theta )

res2=dirf( param, theta, plot=0, print=0 )$dICRF

Print(max(abs(res1-res2)))

comments(
 '
icrf( 0, param )

icrf( -1:1, param )

dirf( param, -1:1, print=0, plot=0 )$dICRF

JacobianMat( -1:1, icrf, param=param )

jacobian( icrf, 0, param=param )

jacobian( icrf, -1:1, param=param )

')
