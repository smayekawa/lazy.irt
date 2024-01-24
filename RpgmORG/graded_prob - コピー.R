#' Calculate the probability contents associated with the integer score
#'
#' @param x The integer score vector
#' @param mean The mena of the underlying contimuous variable
#' @aliases std The standard deviationi of the underlying contimuous variable
#' @param min The minimum value of x
#' @param max The maximum value of x
#' @param truncate = 1 to truncate the range of X between [min,max]
#'
#' @details
#' Let Y be distributed as N( mean, std ).
#' This function evaluates the probabilities:
#'   p[1] = Pr( min < X <= x[1] ), p[2] =  Pr( x[1] < X <= x[2] ),  ...
#'   p[length(x)] =  Pr( x[n] < X <= max ) \cr\cr
#'
#' If truncate = 1, the probabilites will be normalized so that
#'  sum(p) = 1.
#' \cr
#' When truncate = 1, if the interval (min, max) is not wide enough to
#' cover the whole range of X, the resulting probabilities may not
#' sum to unity.
#'
#' @examples
#' # The rightmost/leftmost categoris are wider than the rest.
#' # truncation has no effect
#' graded_prob( 0:4,  2, 2,  truncate=0  )
#' graded_prob( 0:4,  2, 2,  truncate=1  )
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
                         , print=0 ){
 # calculate prob for discrete x
 # Shin-ichi Mayekawa
 # 20180629dnc,30
 #

 # range of x
 if( is.null(min) ) min=mean-50*std
 if( is.null(max) ) max=mean+50*std

 if( truncate ) dd=pnorm( max, mean, std )-pnorm( min, mean, std )
 else dd=1

 brk=midp2brk( x, min=min, max=max )

 cp=pnorm( brk, mean, std )
 pp=cumsuminv(cp)
 p=pp[-1]
 p=p/dd

 if( print >= 1 ){
  Print(mean, std, truncate)
  Print( x, brk, cp, pp, p, fmt="8.3" )
  Print(sum(p), dd)
 }

 return( p )

} # end of graded_prob


%//

# The rightmost/leftmost categoris are wider than the rest.
# truncation has no effect
graded_prob( 0:4,  2, 2,  truncate=0  )
graded_prob( 0:4,  2, 2,  truncate=1  )

# The rightmost/leftmost categoris have the same length as the rest.
graded_prob( 0:4,  2, 2,  min=-0.5, max=4.5, truncate=0  )
graded_prob( 0:4,  2, 2,  min=-0.5, max=4.5, truncate=1  )

