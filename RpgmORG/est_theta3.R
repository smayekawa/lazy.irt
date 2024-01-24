#' Estimation of Theta
#'
#' @param Uc The compressed item response data.
#' @param U  The uncompressed item response data.
#' @param param The parameter data frame.
#' @param method = "EAP", "MAP" or "ML"
#' @param theta Discrete theta values.
#' @param thd pdf of theta normalized to sum to unity. \cr
#' When NULL, N(0,1)  will be used. \cr
#' When thd == 1, locally uniform prior will be used.
#' @param thdist = "NORMAL" or "UNIFORM"
#' @param thmin Minimum value of discrete thata value.
#' @param thmax Maximum value of discrete thata value.
#' @param npoints # of discrete points for theta.
#' @param thmean The prior mean of normal theta distribution.
#' @param thstd   The prior standard deviation of normal theta distribution.
#' @param print = 2 to print the result.
#' @param plot = 1 to show the histogram of theta.
#'
#' @details
#' The core part of EAP is as follows:
#' \preformatted{
#'  ULP=U\%*\%t(logP)
#'  eULP=exp(ULP)*matrix(1,nrow(ULP))\%*\%t(thd)
#'  rseULP=rowSums(eULP)
#'  H=eULP/rseULP
#'  thetahat=H\%*\%theta
#' }
#' where \code{logP} is the log of irf matrix (ntheta x sum(ncat)), \cr
#' \code{U} is the un-compressed item response matrix( N x sum(ncat)), \cr
#' and \code{thd} is the normalized prior theta pdf vector (ntheta). \cr
#' Therefore, \code{H} is N x ntheta matrix of normalized posterior pdf \cr
#' First two lines calculate the unnormalized joint pdf of U and theta. \cr
#' The next line line normalizes the posterior distribution of theta
#' given \code{U} stored in \code{H} matrix. \cr
#' The last line calculates the posterior mean.
#' \cr\cr
#' When \code{thd} is a scalar, the locally uniform theta distribution
#'  is used.
#' \cr\cr
#' method = "ML" or method ="MAP" uses native optimze function and
#' may be very slow.
#'
#' @return A list of \cr
#' thetahat Estimated theta valuesr \cr
#' theta Theta values for EAP \cr
#' thd Theta disrtibution for EAP\cr
#' method  \cr
#' thdist \cr
#' thmin \cr
#' thmax \cr
#' npoints \cr
#' thmean \cr
#' thstd \cr
#'
#' @examples
#'
#' # compressed item response at equally spaced 21 theta points in [-3,3]
#' set.seed(1701)
#' resg <- gendataIRT( 1, paramS3, npoints=21, compress=1 )
#' Uc <- resg$U
#' theta <- resg$theta
#' thetaEAP <- est_theta( Uc, param=paramS3, print=1 )$thetahat
#' thetaMAP <- est_theta( Uc, param=paramS3, print=1, method="MAP" )$thetahat
#' thetaML <- est_theta( Uc, param=paramS3, print=1, method="ML" )$thetahat
#' cor( cbind(theta,thetaEAP,thetaMAP,thetaML) )
#' plot(theta, thetaEAP)
#'
#' # uncompressed item response at equally spaced 101 theta points in [-3,3]
#' set.seed(1701)
#' resg <- gendataIRT( 1, paramS3, npoints=101, compress=0 )
#' U <- resg$U
#' theta <- resg$theta
#' res <- est_theta( U=U, param=paramS3, print=1 )
#' plot(theta,res$thetahat)
#'
#' # 100 normally distriuted theta
#' set.seed(1701)
#' resg <- gendataIRT( 1, paramS3, npoints=100, compress=1, thdist="rnorm" )
#' Uc <- resg$U
#' theta <- resg$theta
#' res <- est_theta( Uc=Uc, param=paramS3, print=1 )
#' plot(theta,res$thetahat)
#'
#' @export
#'

est_theta <- function( Uc=NULL, U=NULL, param=NULL, method="EAP"
                      , theta = NULL, thd=NULL, thdist="NORMAL"
                      , thmin = -4, thmax = 4, npoints = 21
                      , thmean=0, thstd=1
                      , print=0, plot=0 ){
 # estimation of theta
 # Shin-ichi Mayekawa
 # 20151113
 # U: 20161116
 # revised algorithm: 20161202
 # renamed from esttheta: 20121202
 # make U a matrix: 2
 # ML and MAP: 20161213dnc
 # thsd -> thstd: 20161213dnc
 # bugfix: 20161214
 # 720 -> 600: 20171107
 # plot: 20171107
 # other bug fix: 20171107
 # ncat in dummy_expand: 20200103cot
 #


 mllh_theta <- function( theta, u ){
  # minus log likelihood of theta
  # Shin-ichi Mayekawa
  # 20161213dnc
  #
  # theta 1 x 1
  # u        nitem x 1
  #
  # param, method, thmean, and thstd are from calling environment.
  #

  # irf
  logP=log( irf( param, theta, zero=1, print=0 )$ICRF )
  ULP=sum( u*logP )
  if( method == "MAP" ) ULP=ULP - 0.5*(theta-thmean)^2/thstd^2

  return( -ULP )

 } # mllh_theta


 method=toupper(method)
 thdist=toupper(thdist)


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
   cat("\nerror1(est_theta): Both U and Uc are given.\n")
   return()
  }
  if( nitem != ncol(Uc) ){
   cat("\nerror1:(est_theta) # of items in param = ", nitem
       , " but # of columns in ", Uc, " = ", ncol(Uc), "\n" )
   return()
  }
  # uncompress data
  U=dummy_expand( Uc, ncat=ncat )$U
  rm(Uc)

 }
 else if( is.null(U) ){
  cat("\nerror1(est_theta): Either U or Uc must be given.\n")
  return()
 }

 # Here, U is given.

 # make U a matrix
 if( is.data.frame(U) ) U=as.matrix(U)

 # replace NA by 0
 U[is.na(U)]=0


 N=nrow(U)

 if( sum(ncat) != ncol(U) ){
  Print(ncat,U,Uc)
  cat("\nerror1:(est_theta) Sum of ncat in  param = ", sum(ncat)
      , " but # of columns in U = ", ncol(U), "\n" )
  cat(" Responses to the highest category may not be observed.\n")
  return()
 }


 if( print ){
  cat("\nEstimation of Theta\n")
  cat(" # of items =" , nitem, "\n")
  cat(" # of subjects = ", N, "\n")
  cat("\n")
  cat(" estimatio method =", method, "\n")
  Print( thmin, thmax, npoints, thmean, thstd )
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
   if( thdist == "NORMAL" ){
    thd=dnorm( theta, thmean, thstd )
   }
   else thd=rep(1,npoints)
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

  # H=exp( logP%*%t(U)+log(thd) )
  # H=H%*%diag(1/colSums(H))
  # thetahat=t( t(theta)%*%H )

  ULP=U%*%t(logP)
  # For each row, make the minimum value -745 hoping that
  #  -745 < ULP < 709
  rowmin=apply(ULP,1,min)
  rowmax=apply(ULP,1,max)
  ULP=ULP-rowmin  - 600       # instead of 745 -> 720 -> 600
# rowmin1=apply(ULP,1,min)
# rowmax1=apply(ULP,1,max)
# Print(rowmin,rowmax,rowmin1,rowmax1)
  ULP[ULP > 709]=709
  eULP=exp(ULP)*matrix(1,nrow(ULP))%*%t(thd)
  rseULP=rowSums(eULP)
  rseULP[rseULP < 1e-307]=1e-307
  H=eULP/rseULP
  thetahat=H%*%theta
  thetahat=as.vector(thetahat)


 }
 else if( method == "ML" ){

  # item difficulty
  b=param[,grep("^p[[:digit:]]",colnames(param))][,-1]
  b=as.matrix(b)
  meanb=rowSums(b)/(ncat-1)
  # easiest and most difficult items
  loce=which.min(meanb)
  locd=which.max(meanb)
  locd1=ncat[locd]
  # Print(ncat,b,meanb, loce,locd)

  # 0 and ncat-1 category location
  rr=get_range(ncat)



  thetahat=numeric(N)
  for( i in 1:N ){

   if( print >= 1 ){
    if( i %% 100 == 0 ){
     cat(".")
    }
   }

   ui=U[i,]

   if( all( ui[rr[,1]] == 1 ) ) ui[rr[loce,1:2]]=0.5
   if( all( ui[rr[,2]] == 1 ) ) ui[rr[locd,(locd1-1):locd1]]=0.5

   thetahat[i]=optimize( mllh_theta, c(-5,5), u=ui )$minimum
  }

 }
 else if( method == "MAP" ){

  thetahat=numeric(N)
  for( i in 1:N ){

   if( print >= 1 ){
    if( i %% 100 == 0 ){
     cat(".")
    }
   }

   ui=U[i,]
   thetahat[i]=optimize( mllh_theta, c(-5,5), u=ui )$minimum
  }

 }
 else{

  cat( "error1:(est_theta)"
       , "\n\n Method other than 'EAP' is Not Yet Available.\n\n")
  return()

 }


 if( print >= 2 ){
  cat("\n")
  Print(thetahat)
 }

 if( plot > 0 ){
  hist( thetahat )
 }


 res=named_list( thetahat, theta, thd
                 , method, thmin, thmax, npoints, thdist, thmean, thstd )

 return( res )


} # end of est_theta



param=paramA1[5,]
Uc=matrix(0:2,3,1)
est_theta( Uc=Uc, param=param )

param=paramA1[5,]
Uc=matrix(0:2,3,1)
est_theta( Uc, param=param )

U=matrix(c(0,1,0),1,3)
est_theta( U=U, param=param )


comments(
 '


# uncompressed item response at equally spaced 101 theta points in [-3,3]
set.seed(1701)
resg <- gendataIRT( 1, paramS3, npoints=101, compress=0 )
U <- resg$U
theta <- resg$theta
res <- est_theta( U=U, param=paramS3, print=1, method="EAP" )
plot(theta,res$thetahat)







# compressed item response at equally spaced 101 theta points in [-3,3]
set.seed(1701)
resg <- gendataIRT( 1, paramS3, npoints=101, compress=1 )
Uc <- resg$U
theta <- resg$theta
res <- est_theta( Uc, param=paramS3, print=1 )
plot(theta,res$thetahat)

# uncompressed item response at equally spaced 101 theta points in [-3,3]
set.seed(1701)
resg <- gendataIRT( 1, paramS3, npoints=101, compress=0 )
U <- resg$U
theta <- resg$theta
res <- est_theta( U=200*U, param=paramS3, print=1 )
plot(theta,res$thetahat)



# 100 normally distriuted theta
set.seed(1701)
resg <- gendataIRT( 1, paramS3, npoints=100, compress=1, thdist="rnorm" )
Uc <- resg$U
theta <- resg$theta
res <- est_theta( Uc=Uc, param=paramS3, print=1 )
plot(theta,res$thetahat)

res <- est_theta( Uc=Uc, param=paramS3, print=1, thd=1 )
plot(theta,res$thetahat)




')






