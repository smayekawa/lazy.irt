#' Estimation of Theta
#'
#' @param Uc The compressed item response data
#' @param param The parameter data frame
#' @param method = "EAP" or "ML"   (ML is not yet availble.)
#' @param theta Discrete theta values
#' @param thd pdf of theta normalized to sum to unity
#' @param thmin Minimum value of discrete thata value
#' @param thmax Maximum value of discrete thata value
#' @param npoints # of discrete points for theta
#' @param thmean The prior mean of normal theta distribution
#' @param thsd   The prior standard deviation of normal theta distribution
#' @param print = 1 to print the result
#'
#' @details
#' The core part of EAP is as follows:
#' \preformatted{
#' H=exp( logP\%*\%t(U)+log(thd) )
#' H=H\%*\%diag(1/colSums(H))
#' thetahat=t( t(theta)\%*\%H )
#' }
#' where \code{logP} is the log of irf, \code{U} is the un-compressed
#' item response, and \code{thd} is the normalized prior theta pdf. \cr
#' The second line normalizes the posterior distribution of theta
#' given \code{U} stored in \code{H} matrix. \cr
#' The third line calculates the posterior mean.
#'
#' @examples
#'
#' # item response at equally spaced 21 theta points in [-3,3]
#' set.seed(1701)
#' Uc <- gendataIRT( 1, paramS1, npoints=21, thmin=-3, thmax=3, compress=1 )$U
#' res=esttheta( Uc, paramS1, print=1 )
#' plot(rownames(Uc),res$thetahat)
#'
#'
#' # 100 normally distriuted theta
#' set.seed(1701)
#' Uc <- gendataIRT( 1, paramS1, npoints=100, compress=1, thdist="rnorm" )$U
#' res=esttheta( Uc, paramS1, print=1 )
#' plot(rownames(Uc),res$thetahat)
#'
#' @export
#'

esttheta <- function( Uc=NULL, param=NULL, method="EAP"
                      , theta = NULL, thd=NULL
                      , thmin = -4, thmax = 4, npoints = 21
                      , thmean=0, thsd=1
                      , print=0 ){
 # estimation of theta
 # Shin-ichi Mayekawa
 # 20151113
 #

 if( is.null(Uc) ){
  stop("Uc must be given.")
 }
 if( is.null(param) ){
  stop("param must be given.")
 }

 # from parameter data set
 nitem=nrow(param)
 type=param$type
 ncat=param$ncat

 if( nitem != ncol(Uc) ){
  cat("\nerror1:(esttheta) # of items in ", param, " = ", nitems
      , " but # of columns in ", Uc, " = ", ncol(Uc), "\n" )
  return()
 }

 # uncompress data
 N=nrow(Uc)
 U=dummy_expand( Uc )$U

 if( sum(ncat) != ncol(U) ){
  cat("\nerror1:(esttheta) Sum of ncat in ", param, " = ", nitems
      , " but # of columns in U = ", ncol(U), "\n" )
  return()
 }


 if( print ){
  cat("\nEstimation of Theta\n")
  cat(" # of items =" , nitem, "\n")
  cat(" # of subjects = ", N, "\n")
  cat("\n")
  cat(" estimatio method =", method, "\n")
  Print( thmin, thmax, npoints, thmean, thsd )
  }

 if( method == "EAP" ){

  # generate theta
  if( is.null(theta) ){
   theta=seq(thmin,thmax,length.out=npoints)
  }
  else{
   if( is.matrix(theta) ) theta=as.vector(theta)
   npoints=length(theta)
  }

  if( is.null(thd) ){
   thd=dnorm( theta, thmean, thsd )
  }
  thd=thd/sum(thd)

  # irf
  irfres=irf( param, theta, zero=1, print=0 )
  P=irfres$ICRF
  fromP=irfres$fromP
  toP=irfres$toP

  # Print(npoints, P, fromP, toP )

  logP=log(P)

  H=exp( logP%*%t(U)+log(thd) )
  H=H%*%diag(1/colSums(H))
  thetahat=t( t(theta)%*%H )

  if( print ){
   Print( Uc, thetahat)
  }

  res=named_list( thetahat, method, thmin, thmax, npoints, thmean, thsd )

  return( res )

 }
 else{

  cat( "error1:(esttheta)"
      , "\n\n Method other than 'EAP' is Not Yet Available.\n\n")
  return()

 }


} # end of esttheta




# item response at equally spaced 21 theta points in [-3,3]
set.seed(1701)
Uc <- gendataIRT( 1, paramS1, npoints=21, thmin=-3, thmax=3, compress=1 )$U
res=esttheta( Uc, paramS1, print=1 )


# 100 normally distriuted theta
set.seed(1701)
Uc <- gendataIRT( 1, paramS1, npoints=100, compress=1, thdist="rnorm" )$U
res=esttheta( Uc, paramS1, print=1 )

plot(rownames(Uc),res$thetahat)

Print(unique(Uc))












