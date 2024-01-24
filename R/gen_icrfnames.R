#' Generate the Row Names of vec(icrf) as the Combination of Theta Values
#' and Category Names.
#'
#' @param theta Numeric vector of theta values.
#' @param ncat  # of categories
#' @param zero = 1 to include the 0-th category
#' @param cat.first  = 1 to chage the category fist. \cr
#' = 2 to change the category first and place category name before theta. \cr
#' This is for vec(icrf, byrow=1 ).
#' @param digits # of digits after the decimal point of theta.
#' @param catroot The root character for category.
#' @param sep The caracter separating theta and category.
#' @param compress = 1 to remove all the blanks in the result
#'
#' @examples
#' gen_icrfnames(seq(-4,4,length=7), 3 )
#' gen_icrfnames(seq(-4,4,length=7), 3, cat.first=1 )
#'
#' @details
#' This function creates a character vector consisting of \cr
#'   \code{ c( thetavalue1:catvalue1,  thetavalue2:catvalue1, ... ) }. \cr
#' If cat.first=1, the result is \cr
#'  \code{ c( thetavalue1:catvalue1,  thetavalue1:catvalue2, ... ) }, \cr
#' and if cat.first=2, \cr
#'  \code{ c( catvalue1:thetavalue1,  catvalue2:thetavalue1, ... ) }, \cr
#'
#' @return
#' A vector consisting of rownames of vec(icrf) or vec(icrf, byrow=1).
#'
#' @export
#'

gen_icrfnames <- function( theta, ncat, zero=0, cat.first=0, digits=2
                           , catroot="c", sep=":", compress=0 ){
 # generate the rownames of vec(icrf)
 # Shin-ichi Mayekawa
 # Original version in dicrf_p: 20150615
 # 20180212
 #

 thname=formatC(theta,digits=digits, format="f")
 cc=(1-zero):(ncat-1)
 catname=paste( catroot, cc, sep="" )
 if( cat.first == 1 )
  rnJ=apply( expand.grid1( thname, catname, rev=1 ), 1, paste, collapse=sep )
 else if( cat.first == 2 )
  rnJ=apply( expand.grid( catname, thname ), 1, paste, collapse=sep )
 else
  rnJ=apply( expand.grid( thname, catname ), 1, paste, collapse=sep )

 if( compress ) rnJ=trim( rnJ, all=2 )
 return( rnJ )

} # end of gen_icrfnames


