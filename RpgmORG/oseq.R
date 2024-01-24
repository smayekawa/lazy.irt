#' IRT Observed Score Equating
#'
#' @param param1 The item parameter data frame for test1
#' @param param2 The item parameter data frame for test2
#' @param weight1 The weight data frame for test1, or NULL
#' @param weigth2 The weigth data frame for test2, or NULL
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
#' \code{lazy.irt::sumsmnw}.
#'
#' The equipercentile equating will be performed by \code{lazy.tools::eq}.
#'
#'
#' @return A list of \cr
#' theta The value of theta corresponding to x. \cr
#' x1 The test score of test1 \cr
#' x2_1 The test2 equivalent score of x1
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
  #


  obscore_s <- function( param, weight
                         , thmin=-4, thmax=4, npoints=31, thdist=1 ){
    # Simplest verssion of lazy.irt::obscore.
    # Shin-ichi Mayekawa
    # 20230717
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

    return( named_list(score,Px) )

  } # end of obscore_s


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
  p1=temp$Px
  cdf1=cumsum(p1)/sum(p1)
  temp=obscore_s( param2, w2
                , thmin=thmin, thmax=thmax, npoints=npoints, thdist=thdist )
  x2=temp$score
  p2=temp$Px
  cdf2=cumsum(p2)/sum(p2)

  if( print >= 3 ){
   cat("\n Marginal Distribution of Test Scores\n")
   cat(" theta distribution is ", thdist,"\n")
   Print(x1,p1,x2,p2)
  }

  # equipercentile equating
  res=eq( score1=x1, cdf1=cdf1, score2=x2, cdf2=cdf2
        , smooth1=0, smooth2=0
        , print=0, plot=0 )$ctable0
  x2_1=res[,2]

  if( round ) x2_1=round(x2_1,0)

  if( print > 1){
    Print(x1, x2_1)
  }

  if( plot ){
    plot( x1, x2_1, type="b", main="IRT Observed Score Equating")
  }

  return( named_list( x1, x2_1, minmax1, minmax2) )

} # end of oseq



res=oseq( paramS1, paramS2, weight1=weightS12, weight2=weightS21, print=3 )





comments(
'



res1=oseq( paramB1, paramB1, weight1=weightB11, weight2=weightB12
           , print=2, plot=1 )

res2=tseq( paramB1, paramB1, weight1=weightB11, weight2=weightB12
           , print=2, plot=1 )







# checking if symmetric
# equate test1 to test2
res1=oseq( paramS1, paramS2, by_x=1, round=1 )
# equate test2 to test1
res2=oseq( paramS2, paramS1, by_x=1, round=1 )
# merge result
r1=data.frame(x1=res1$x, y1=res1$y)
r2=data.frame(y2=res2$x, x2=res2$y)
rx=merge( r1, r2, by.x="x1", by.y="x2", all=TRUE)
rx$y=ifelse( is.na(rx$y1), rx$y2, rx$y1 )
ry=merge( r1, r2, by.x="y1", by.y="y2", all=TRUE)
ry$x=ifelse( is.na(ry$x1), ry$x2, ry$x1 )
printm(rx,ry)




 res1=oseq( paramS1, paramS2, round=1
            , by_x=0.5, weight1=weightS12, weight2=weightS21 )
 res2=oseq( paramS2, paramS1, round=1
            , by_x=0.5, weight2=weightS12, weight1=weightS21 )
 r1=data.frame(x1=res1$x, y1=res1$y)
 r2=data.frame(y2=res2$x, x2=res2$y)
 rx=merge2( r1, r2, by.x="x1", by.y="x2", all=TRUE)
 ry=merge2( r1, r2, by.x="y1", by.y="y2", all=TRUE)





 res1=oseq( paramS1, paramS2, round=1
            , by_x=0.5, weight1=weightS12, weight2=weightS21 )
 res2=oseq( paramS2, paramS1, round=1
            , by_x=0.5, weight2=weightS12, weight1=weightS21 )
 r1=data.frame(x1=res1$x, y1=res1$y)
 r2=data.frame(y2=res2$x, x2=res2$y)
 rx=merge2( r1, r2, by.x="x1", by.y="x2", all=TRUE)
 ry=merge2( r1, r2, by.x="y1", by.y="y2", all=TRUE)



'
)
