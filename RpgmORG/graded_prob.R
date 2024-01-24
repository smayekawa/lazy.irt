#' Calculate the probability contents associated with the integer score
#'
#' @param x The integer score vector
#' @param mean The mena of the underlying continuous variable
#' @param std The standard deviationi of the underlying continuous variable
#' @param min The minimum value of x
#' @param max The maximum value of x
#' @param truncate = 1 to truncate the range of X between [min,max]
#' @param method = 1 to use pnorm(x[i+1])-pnorm(x[i]) \cr
#' = 0 to set the probabilities proportional to dnorm(x)
#' @param print = 1 to print the result
#' @param plot = 1 to plot the result
#'
#'
#' @details
#' Let Y be distributed as N( mean, std ).
#' This function evaluates the probabilities:
#'   p[1] = Pr( min < X <= x[1] ), p[2] =  Pr( x[1] < X <= x[2] ),  ...
#'   p[length(x)] =  Pr( x[n] < X <= max ) \cr\cr
#'
#' If method = 1, \cr
#'  p[i] = pnorm( b[i+1], mean, std ) - pnorm( b[i], mean, std )  \cr
#' else \cr
#'  p[i] = dnorm( x[i], mean, std ) \cr
#'
#' If method = 2 or method = 1 and truncate = 1, \cr
#' the probabilites will be normalized so that  sum(p) = 1.
#' \cr
#' When truncate = 1, if the interval (min, max) is not wide enough to
#' cover the whole range of X, the resulting probabilities may not
#' sum to unity.
#'
#' @examples
#' # The rightmost/leftmost categoris are wider than the rest.
#' # truncation has no effect
#' graded_prob( 0:4,  2, 2,  truncate=0, print=1  )
#' graded_prob( 0:4,  2, 2,  truncate=1, print=1  )
#' graded_prob( 0:4,  2, 2,  method=2, print=1  )
#'
#' # The rightmost/leftmost categoris have the same length as the rest.
#' graded_prob( 0:4,  2, 2,  min=-0.5, max=4.5, truncate=0  )
#' graded_prob( 0:4,  2, 2,  min=-0.5, max=4.5, truncate=1  )
#'
#' @return
#'  p A vector (same size as x)
#'
#' @export
#'

graded_prob <- function( x, mean, std, min=NULL, max=NULL, truncate=0
                         , method=1, print=0, plot=0 ){
 # calculate prob for discrete x
 # Shin-ichi Mayekawa
 # 20180629dnc,30
 # method: 20180703dnc,04dnc
 #

 # range of x
 if( is.null(min) ) min=mean-50*std
 if( is.null(max) ) max=mean+50*std


 if( method == 2 ){
  p=dnorm( x, mean, std )
  p=p/sum(p)
  dd=1; brk=NA; cp=NA; pp=NA
 }
 else{
  if( truncate ) dd=pnorm( max, mean, std )-pnorm( min, mean, std )
  else dd=1

  brk=midp2brk( x, min=min, max=max )

  cp=pnorm( brk, mean, std )
  pp=cumsuminv(cp)
  p=pp[-1]
  p=p/dd
 }

 if( print >= 1 ){
  Print(mean, std, truncate, method)
  Print( x, brk, cp, pp, p, fmt="8.3" )
  Print(sum(p), dd)
 }

 if( plot == 1 ){
  title=paste("mean =", mean,", std =", std,",  method =", method
              , ", truncate =", truncate )
  bp = barplot( height = p, names.arg = x
                , main=title, space = 0, las = 2 )
 }

 return( p )

} # end of graded_prob


