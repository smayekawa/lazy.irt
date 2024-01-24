#' Calculate Observed Score Distribution (simple version of \code{obscore})
#'
#' @param param Parameter data frame
#' @param weight Weight data frame
#' @param thmin Minimum value of discrete thata value
#' @param thmax Maximum value of discrete thata value
#' @param npoints # of discrete points for theta
#' @param thdist Type of theta distribution \cr
#' = 0 to use uniform,  \cr = 1 to use N(0,1)
#' @param print = 1 to print the result
#'
#' @details
#' This function calculates the marginal distribution of the test score. \cr
#' \preformatted{
#'  1. Calculate the conditional distribution of X given theta
#'  2. Calculate the joint distribution of X and theta.
#'  3. Calcukate the marginal distribution of X.
#' }
#' In step 1 above, \code{icrf} is calculated by \code{lazy.irt::irf}
#' and the distribution of the weighted sum of scored multinomial distributions
#' will be calculated by \code{lazy.irt::sumsmnw}.
#'
#' @examples
#' res <- obscore_s( paramS2, weightS21, print=1 )
#'
#' @export
#'

obscore_s <- function( param, weight
                       , thmin=-4, thmax=4, npoints=31, thdist=1, print=0 ){
  # Simplest verssion of lazy.irt::obscore.
  # Shin-ichi Mayekawa
  # 20230717
  # roxygen2: 20231031dnc
  #

  # constants
  ncat=param$ncat
  nitems=nrow(param)
  w=weight$w
  v=weight[,regexpr("^v", colnames(weight)) > 0, drop=0]

  # generate thata and prior theta dist
  theta=seq(thmin,thmax,length.out=npoints)
  thname=format(theta,digits=2)
  if( thdist == 0 ) Pt=matrix(1/npoints,npoints)
  else if( thdist == 1 ){
    Pt=exp(-0.5*theta^2);
    Pt=Pt/sum(Pt)
  }

  # calculate icrf
  temp=irf( param, theta, weight, print=0, debug=0, plot=0 )
  icrf=temp$ICRF
  fromP=temp$fromP
  toP=temp$toP
  maxscore_t=temp$maxscore_t
  score=0:maxscore_t
  rm(temp)

  # conditional dist of x given theta
  # Px_t=matrix(0,maxscore_t+1,npoints)
  Px_t=NULL
  for( k in 1:npoints ){
    Pk=matrix(NA,max(ncat),nitems)
    for( j in 1:nitems ){
      Pk[1:ncat[j],j]=t( icrf[k,fromP[j]:toP[j]] )
    }
    Px_t=cbind( Px_t
                , sumsmnw( Pk, t(v), w, compress=0
                           , print=0, plot=0, debug=0 )[,2] )
  }
  rownames(Px_t)=score; colnames(Px_t)=thname

  # joint
  Pxt=Px_t%*%Diag(Pt)
  Pxt[is.nan(Pxt)]=NA

  # marginal x
  Px=matrix(rowSums(Pxt),,1)

  if( print ){
    printm( param, weight )
    Print( score, Px )
  }

  return( named_list(score,Px) )

} # end of obscore_s


