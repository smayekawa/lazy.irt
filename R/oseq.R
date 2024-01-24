#' IRT Observed Score Equating
#'
#' This function performs the IRT Observed Score Equating of
#' two test scores. \cr
#' Japanese help file: (\link[lazy.irt]{oseq_JPH})
#'
#' @param param1 The item parameter data frame for test1
#' @param param2 The item parameter data frame for test2
#' @param weight1 The weight data frame for test1, or NULL
#' @param weight2 The weigth data frame for test2, or NULL
#' @param thmin The minimum value of theta to calculate
#' the marginal score distribution.
#' @param thmax  The maximum value of theta to calculate
#' the marginal score distribution.
#' @param npoints # of theta points in \code{[thmin,thmax]}
#' @param thdist = 1 to use normal theta, = 0 to use uniform theta.
#' @param smooth1	# of times to smooth cdf1.
#' @param bandwid1 bandwidth for running average smooth of cdf1.
#' @param smooth2	# of times to smooth cdf2.
#' @param bandwid2 bandwidth for running average smooth of cdf2.
#' @param round = 1 to round the result to integer
#' @param print = 1 or 2 to print the result
#' @param plot = 1 to plot the result
#'
#'
#' @details
#' This function equates the test score of test1, x1, to the score of test2,
#' x2. by the IRT Observed Score Equating.\cr
#' The resulting score will be named as x2_1.
#'
#' In irt observed score equating, the conditional distribution of X
#' given theta, Px_t, is calculated first.
#' Then, using the marginal distribution of theta, P_t,
#' the joint distribution of X and theta, Pxt, is calculated.
#' Finally, the marginal distribution of X, Px, is calculated from Pxt. \cr
#' Given (x1, Px1)  and (x2, Px2), equipercentile equating will be used
#' to equate x1 to x2.
#'
#'
#' Note that, when \code{round=1}, the result may NOT be symmetric.\cr
#' Run the example below with \code{by_x=1} and \code{round=1}.
#'
#' The marginal distribution of the test score will be calculated by
#' \link[lazy.irt]{obscore_s} which is a simpler form of 
#' \link[lazy.irt]{obscore}.
#'
#' The equipercentile equating will be performed by \link[lazy.irt]{coseq}.
#'
#'
#' @return A list of \cr
#' x1 The test score of test1. \cr
#' x2 The test score of test2. \cr
#' x2_1 The test2 equivalent score of x1. \cr
#' p1 The marginal distribution of test1 score. \cr
#' p2 The marginal distribution of test2 score.
#'
#'
#' @examples
#' res=oseq( paramS1, paramS2, print=2 )
#' res=oseq( paramS1, paramS2, weight1=weightS12, weight2=weightS21, print=2 )
#'
#' # The effect of item weight.
#' res=oseq( paramS1, paramS1, weight2=weightS12, print=3 )
#'
#' # checking if symmetric
#' # equate test1 to test2
#' res1 <- oseq( paramS1, paramS2 )
#' # equate test2 to test1
#' res2 <- oseq( paramS2, paramS1 )
#' # merge result
#' r1 <- data.frame(x1=res1$x1, y1=res1$x2_1)
#' r2 <- data.frame(y2=res2$x1, x2=res2$x2_1)
#' rx <- merge( r1, r2, by.x="x1", by.y="x2", all=TRUE)
#' rx$y <- ifelse( is.na(rx$y1), rx$y2, rx$y1 )
#' ry <- merge( r1, r2, by.x="y1", by.y="y2", all=TRUE)
#' ry$x <- ifelse( is.na(ry$x1), ry$x2, ry$x1 )
#' rxr <- round(rx,4); ryr=round(ry,4)
#' printm(rxr,ryr)
#'
#' @export
#'

oseq <- function( param1=NULL, param2=NULL, weight1=NULL, weight2=NULL
                  , thmin=-4, thmax=4, npoints=31, thdist=1
                  , smooth1 = 0, bandwid1 = 3, smooth2 = 0, bandwid2=3
                  , round=0, print=0, plot=0 ){
  # irt observed score equating
  # Shin-ichi Mayekawa
  # 20230714
  # create_weight_df: 20230715cot
  # weight1,2: 20230716cot
  # delete obscore_s: 20230721
  # JPH: 20231030
  # x2 added to output: 20231030
  # eq -> coseq: 20231031
  #

  # constants
  nitem1=nrow(param1); ncat1=param1$ncat; iname1=param1$name; itype1=param1$type
  nitem2=nrow(param2); ncat2=param2$ncat; iname2=param2$name; itype2=param2$type

  # create weight data frame from item parameter data frame
  if( is.null(weight1) ) w1=create_weight_df( param1 )$weight
  else w1=weight1
  if( is.null(weight2) ) w2=create_weight_df( param2 )$weight
  else w2=weight2

  minmax1=find_minmax_score( w1 )$minmax
  minmax2=find_minmax_score( w2 )$minmax

  if( print ){
    cat("\nIRT Observed Score Equating \n")
    Print( nitem1,minmax1, nitem2, minmax2)
    if( print >= 3 ){
      printm(param1,w1, minmax1)
      printm(param2,w2, minmax2)
    }
  }

  # marginal distribution of test scores
  temp=obscore_s( param1, w1
                , thmin=thmin, thmax=thmax, npoints=npoints, thdist=thdist )
  x1=temp$score
  p1=c(temp$Px)
  names(p1)=x1
  cdf1=cumsum(p1)/sum(p1)
  temp=obscore_s( param2, w2
                , thmin=thmin, thmax=thmax, npoints=npoints, thdist=thdist )
  x2=temp$score
  p2=c(temp$Px)
  names(p2)=x2
  cdf2=cumsum(p2)/sum(p2)

  if( print >= 3 ){
   cat("\n Marginal Distribution of Test Scores\n")
   cat(" theta distribution is ", thdist,"\n")
   Print(x1,p1,x2,p2)
  }

  # equipercentile equating
  res=coseq( score1=x1, cdf1=cdf1, score2=x2, cdf2=cdf2
        , smooth1=0, smooth2=0
        , print=0, plot=0 )$ctable0
  x2_1=res[,2]

  if( round ) x2_1=round(x2_1,0)

  if( print > 1){
    cat("\nConversion table from test1 to test2\n")
    Print(x1, x2_1)
  }

  if( plot ){
    plot( x1, x2, x2_1, type="b", main="IRT Observed Score Equating")
  }

  return( named_list( x1, x2_1, minmax1, minmax2, p1, p2) )

} # end of oseq


