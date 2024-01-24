#' IRT True Score Equating
#'
#' @param param1 The item parameter data frame for test1
#' @param param2 The item parameter data frame for test2
#' @param weight1 The weight data frame for test1, or NULL
#' @param weigth2 The weigth data frame for test2, or NULL
#' @param by_x The increment of x
#' @param round = 1 to round the result to integer
#' @param method = 0 to use interpolation when calculating inv trf \cr
#' = 1 to use native \code{uniroot} function (slow).
#' @param interpol Interpolation option: See \code{lazy.tools::interpol}
#' @param thmin The minimum value of theta for interpolation.
#' @param thmax  The maximum value of theta for interpolation.
#' @param npoints # of theta points in \code{[thmin,thmax]} for interpolation.
#' @param print = 1 or 2 to print the result
#' @param plot = 1 to plot the result
#'
#'
#' @details
#' @details
#' This function equates the test score of test1, x1, to the score of test2,
#' x2. \cr
#' The resulting score will be named as x2_1.
#'
#' The score of test1 at \cr
#' \code{x1 <- seq(minscore,maxscore, by=by_x)} \cr
#' will be equated to the score of test2.
#'
#' Note that, when \code{round=1}, the result may NOT be symmetric.\cr
#' Run the example below with \code{by_x=1} and \code{round=1}.
#'
#' @return A list of \cr
#' theta The value of theta corresponding to x. \cr
#' x1 The test score of test1 \cr
#' x2_1 The test2 equivalent of x
#'
#' @examples
#' res <- tseq( paramS1, paramS2, print=2 )
#' res <- tseq( paramS1, paramS2, weight1=weightS12, weight2=weightS21, print=2 )
#'
#' # The effect of item weight.
#' res <- tseq( paramB1, paramB1, weight1=weightB11, weight2=weightB12, print=3 )
#'
#' # checking if symmetric
#' # equate test1 to test2
#' res1 <- tseq( paramS1, paramS2, by_x=0.5 )
#' # equate test2 to test1
#' res2 <- tseq( paramS2, paramS1, by_x=0.5 )
#' # merge result
#' r1 <- data.frame(x1=res1$x1, y1=res1$x2_1)
#' r2 <- data.frame(y2=res2$x1, x2=res2$x2_1)
#' rx <- merge( r1, r2, by.x="x1", by.y="x2", all=TRUE)
#' rx$y <- ifelse( is.na(rx$y1), rx$y2, rx$y1 )
#' ry <- merge( r1, r2, by.x="y1", by.y="y2", all=TRUE)
#' ry$x <- ifelse( is.na(ry$x1), ry$x2, ry$x1 )
#' printm(rx,ry)
#'
#' # comparison of methods
#' res1 <- tseq( paramS1, paramS1, method=1, by_x=0.1, weight2=weightS12 )
#' res2 <- tseq( paramS1, paramS1, method=0, by_x=0.1, weight2=weightS12 )
#' # Print(res1$x1, res1$x2_1, res1$x2_1-res2$x2_1, fmt="6.1, 6.4 9.6")
#' summary( abs(res1$x2_1-res2$x2_1) )
#'
#'
#' @export
#'

tseq <- function( param1=NULL, param2=NULL, weight1=NULL, weight2=NULL
                , method=0, interpol="spline"
                , thmin=-5, thmax=5, npoints=121
                , round=0, by_x=1, print=0, plot=0 ){
 # irt true score equating
 # Shin-ichi Mayekawa
 # 20230714
 # create_weight_df: 20230715cot
 # weight1,2: 20230716cot
 # y -> x2_1: 20230717
 # print: 20230721
 #



 # constants
 nitem1=nrow(param1)
 ncat1=param1$ncat
 iname1=param1$name
 itype1=param1$type
 nitem2=nrow(param2)
 ncat2=param2$ncat
 iname2=param2$name
 itype2=param2$type

 # create weight data frame from item parameter data frame
 if( is.null(weight1) ) w1=create_weight_df( param1 )$weight
 else w1=weight1
 if( is.null(weight2) ) w2=create_weight_df( param2 )$weight
 else w2=weight2

 minmax1=find_minmax_score( w1 )$minmax
 minmax2=find_minmax_score( w2 )$minmax

 if( print ){
  cat("\nIRT True Score Equating \n")
  Print( nitem1,minmax1, nitem2, minmax2)
  if( print >= 3 ){
   printm(param1,w1, minmax1)
   printm(param2,w2, minmax2)
  }
 }


 # test1 score
 x1=seq(minmax1[1],minmax1[2],by_x)
 # theta at x by invtrf
 th=invtrf( param1, x1, weight=w1, method=method, interpol=interpol
          , thmin=thmin, thmax=thmax, npoints=npoints )$theta

 # exclude the scores at min and max
 locmin=which(th == -9999)
 locmax=which(th == 9999)
 th1=th[-c(locmin,locmax)]
 theta=x1
 theta[locmin]=NA
 theta[locmax]=NA
 theta[-c(locmin,locmax)]=th1

 # trf of test 2 at th1
 y1=irf( param2, th1, weight=w2, print=0 )$TRF
 x2_1=x1
 x2_1[locmin]=minmax2[1]
 x2_1[locmax]=minmax2[2]
 x2_1[-c(locmin,locmax)]=y1

 if( round ) x2_1=round(x2_1,0)

 if( print > 1){
   cat("\nConversion table from test1 to test2\n")
   Print(theta, x1, x2_1)
 }

 if( plot ){
  plot( x1, x2_1, type="b", main="IRT True Score Equating" )
 }

 return( named_list( x1, x2_1, theta, locmin, locmax, minmax1, minmax2) )

} # end of tseq





# comparison of methods
res1=tseq( paramS1, paramS1, method=1, by_x=0.1, weight2=weightS12 )
res2=tseq( paramS1, paramS1, method=0, by_x=0.1, weight2=weightS12 )
# Print(res1$x1, res1$x2_1, res1$x2_1-res2$x2_1, fmt="6.1, 6.4 9.6")
summary( abs(res1$x2_1-res2$x2_1) )
res3=tseq( paramS1, paramS1, method=0, by_x=0.1
            , thmin=-4, thmax=4, npoints=201, weight2=weightS12 )
# Print(res1$x1, res1$x2_1, res1$x2_1-res3$x2_1, fmt="6.1, 6.4 9.6")
summary( abs(res1$x2_1-res3$x2_1) )






comments(
'



# checking if symmetric
# equate test1 to test2
res1=tseq( paramS1, paramS2, by_x=1, round=1 )
# equate test2 to test1
res2=tseq( paramS2, paramS1, by_x=1, round=1 )
# merge result
r1=data.frame(x1=res1$x, y1=res1$y)
r2=data.frame(y2=res2$x, x2=res2$y)
rx=merge( r1, r2, by.x="x1", by.y="x2", all=TRUE)
rx$y=ifelse( is.na(rx$y1), rx$y2, rx$y1 )
ry=merge( r1, r2, by.x="y1", by.y="y2", all=TRUE)
ry$x=ifelse( is.na(ry$x1), ry$x2, ry$x1 )
printm(rx,ry)




 res1=tseq( paramS1, paramS2, round=1
            , by_x=0.5, weight1=weightS12, weight2=weightS21 )
 res2=tseq( paramS2, paramS1, round=1
            , by_x=0.5, weight2=weightS12, weight1=weightS21 )
 r1=data.frame(x1=res1$x, y1=res1$y)
 r2=data.frame(y2=res2$x, x2=res2$y)
 rx=merge2( r1, r2, by.x="x1", by.y="x2", all=TRUE)
 ry=merge2( r1, r2, by.x="y1", by.y="y2", all=TRUE)





 res1=tseq( paramS1, paramS2, round=1
            , by_x=0.5, weight1=weightS12, weight2=weightS21 )
 res2=tseq( paramS2, paramS1, round=1
            , by_x=0.5, weight2=weightS12, weight1=weightS21 )
 r1=data.frame(x1=res1$x, y1=res1$y)
 r2=data.frame(y2=res2$x, x2=res2$y)
 rx=merge2( r1, r2, by.x="x1", by.y="x2", all=TRUE)
 ry=merge2( r1, r2, by.x="y1", by.y="y2", all=TRUE)



'
)
