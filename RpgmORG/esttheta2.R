#' Estimation of Theta
#'
#' @param Uc The compressed item response data.
#' @param U  The uncompressed item response data.
#' @param param The parameter data frame.
#' @param method = "EAP" or "ML"   (ML is not yet availble.)
#' @param theta Discrete theta values.
#' @param thd pdf of theta normalized to sum to unity. \cr
#' When NULL, N(0,1)  will be used. \cr
#' When thd == 1, locally uniform prior will be used. \cr
#' @param thmin Minimum value of discrete thata value.
#' @param thmax Maximum value of discrete thata value.
#' @param npoints # of discrete points for theta.
#' @param thmean The prior mean of normal theta distribution.
#' @param thsd   The prior standard deviation of normal theta distribution.
#' @param print = 1 to print the result.
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
#' \cr\cr
#' When \code{thd} is a scalar, the locally uniform theta distribution is used.
#'
#' @return A list of \cr
#' thetahat Estimated theta valuesr \cr
#' theta Theta values for EAP \cr
#' thd Theta disrtibution for EAP\cr
#' method  \cr
#' thmin \cr
#' thmax \cr
#' npoints \cr
#' thmean \cr
#' thsd \cr
#'
#' @examples
#'
#' # compressed item response at equally spaced 101 theta points in [-3,3]
#' set.seed(1701)
#' resg <- gendataIRT( 1, paramS3, npoints=101, compress=1 )
#' Uc <- resg$U
#' theta <- resg$theta
#' res <- esttheta( Uc, param=paramS3, print=1 )
#' plot(theta,res$thetahat)
#'
#' # uncompressed item response at equally spaced 101 theta points in [-3,3]
#' set.seed(1701)
#' resg <- gendataIRT( 1, paramS3, npoints=101, compress=0 )
#' U <- resg$U
#' theta <- resg$theta
#' res <- esttheta( U=U, param=paramS3, print=1 )
#' plot(theta,res$thetahat)
#'
#' # 100 normally distriuted theta
#' set.seed(1701)
#' resg <- gendataIRT( 1, paramS3, npoints=100, compress=1, thdist="rnorm" )
#' Uc <- resg$U
#' theta <- resg$theta
#' res <- esttheta( Uc=Uc, param=paramS3, print=1 )
#' plot(theta,res$thetahat)
#'
#' @export
#'

esttheta <- function( Uc=NULL, U=NULL, param=NULL, method="EAP"
                      , theta = NULL, thd=NULL
                      , thmin = -4, thmax = 4, npoints = 21
                      , thmean=0, thsd=1
                      , print=0 ){
 # estimation of theta
 # Shin-ichi Mayekawa
 # 20151113
 # U: 20161116
 #

 # error
 if( is.null(param) ){
  stop("param must be given.")
 }

 # from parameter data set
 nitem=nrow(param)
 type=param$type
 ncat=param$ncat

 # item response data
 if( !is.null(Uc) ){
  # Uc is given.

  # error
  if( !is.null(U) ){
   cat("\nerror1(esttheta): Both U and Uc are given.\n")
   return()
  }
  if( nitem != ncol(Uc) ){
   cat("\nerror1:(esttheta) # of items in ", param, " = ", nitems
       , " but # of columns in ", Uc, " = ", ncol(Uc), "\n" )
   return()
  }
  # uncompress data
  U=dummy_expand( Uc )$U
  rm(Uc)

 }
 else if( is.null(U) ){
  cat("\nerror1(esttheta): Either U or Uc must be given.\n")
  return()
 }

 # Here, U is given.

 N=nrow(U)

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
  else if( all(thd == 1) ) thd=rep(1,npoints)
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
   Print(thetahat)
  }

  res=named_list( thetahat, theta, thd
                  , method, thmin, thmax, npoints, thmean, thsd )

  return( res )

 }
 else{

  cat( "error1:(esttheta)"
      , "\n\n Method other than 'EAP' is Not Yet Available.\n\n")
  return()

 }


} # end of esttheta




# compressed item response at equally spaced 101 theta points in [-3,3]
set.seed(1701)
resg <- gendataIRT( 1, paramS3, npoints=101, compress=1 )
Uc <- resg$U
theta <- resg$theta
res <- esttheta( Uc, param=paramS3, print=1 )
plot(theta,res$thetahat)

# uncompressed item response at equally spaced 101 theta points in [-3,3]
set.seed(1701)
resg <- gendataIRT( 1, paramS3, npoints=101, compress=0 )
U <- resg$U
theta <- resg$theta
res <- esttheta( U=U, param=paramS3, print=1 )
plot(theta,res$thetahat)



# 100 normally distriuted theta
set.seed(1701)
resg <- gendataIRT( 1, paramS3, npoints=100, compress=1, thdist="rnorm" )
Uc <- resg$U
theta <- resg$theta
res <- esttheta( Uc=Uc, param=paramS3, print=1 )
plot(theta,res$thetahat)

res <- esttheta( Uc=Uc, param=paramS3, print=1, thd=1 )
plot(theta,res$thetahat)











