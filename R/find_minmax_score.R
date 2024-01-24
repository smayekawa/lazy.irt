#' Find the possible minimum and maximum scores of the test
#' from the wei ght data frame
#'
#' @param weight The weight data frame
#' @param print = 1 to print the result
#'
#' @return
#' minmax A vector containing the minimum and the maximum score of the test.\cr
#' minmax_i A matrix containing the minimum and the maximum score of each item
#' before applying the item weight \code{w}.
#'
#' @examples
#' find_minmax_score( weightS21, 1 )
#'
#' @export
#'

find_minmax_score <- function( weight, print=0 ){
 # find the min and max score from weight data frame
 # Shin-ichi Mayekawa
 # 20230716cot
 # use ncat: 20230718
 # rogygen2 fix: 20231030dnc
 #

 # constants
 nitem=nrow(weight)
 ncat=weight$ncat
 w=weight$w
 v=as.matrix( weight[,regexpr("^v", colnames(weight)) > 0, drop=0] )
 # remove excess v
 for( j in 1:nitem ){
  if( ncat[j] < ncol(v)) v[j,(ncat[j]+1):ncol(v)]=NA
 }
 iname=weight$name

 # max score
 maxscore_i=as.matrix(apply(v,1,max,na.rm=1),nitems)
 rownames(maxscore_i)=iname
 maxscore=sum(w*maxscore_i)
 minscore_i=as.matrix(apply(v,1,min,na.rm=1),nitems)
 rownames(minscore_i)=iname
 minscore=min(w*minscore_i)
 minmax=c(minscore,maxscore)
 names(minmax)=c("minscore","maxscore")
 minmax_i=cbind(minscore_i,maxscore_i)
 colnames(minmax_i)=c("min","max")

 if( print ){
  Print( minmax )
  Print( minmax_i, w, v, ncat)
 }

 return( named_list( minmax, minmax_i) )

} # end of find_minmax_score


