#' Observed Score Equating
#'
#' This performs the observed score equating of
#' the test score 1 to the test score 2.
#'
#' @param score1 a vector consisting of the score of test1
#' @param freq1  a vector consisting of the frequency counts at score1
#' @param cdf1   a vector consisting of the cumulative frequencies at score1
#' @param score2 a vector consisting of the score of test2
#' @param freq2  a vector consisting of the frequency counts at score2
#' @param cdf2   a vector consisting of the cumulative frequencies at score2
#'
#' @param lim1 min and max score of test1
#' @param lim2 min and max score of test2
#'
#' @param smooth1 # of times to smooth cdf1
#' @param bandwid1 bandwidth for running average smooth of cdf1
#' @param smooth2 # of times to smooth cdf2
#' @param bandwid2 bandwidth for running average smooth of cdf2
#'
#'
#' @param method = 1 to use linear equating \cr
#'               = 2 to use equi-percentile equating
#'
#' @param interpol_method = "constant", "linear" or "spline"
#'
#' @param title  title string
#' @param nolow  = 1 to avoid lowering scores
#' @param round  = 0 not to round the result to integer. (not yet available)
#'
#' @param print = 1 to print the result
#' @param plot = 1 to plot the conversion table \cr
#'             = 2 to plot the cdf \cr
#'             = 3 to plot the smoothed cdf
#'
#' @details
#' Equipercentile equating of test1 score \code{x1} to test2 score \code{x2}
#'  is defined as\cr
#' \preformatted{
#'  x21 = invF2( F1(x1) )
#' }
#' where \code{F1} is the distribution function of \code{x1} and
#' \code{invF2} is the inverse of the distribution function of \code{x2}.
#' In this function, \code{invF2} is calculated by interpolating
#' \code{( F2(x2), x2 )} at \code{F1(x1)} with or without smoothing.
#'
#' Equipercentile equating is essentially the same as native \code{qqplot}. \cr
#' Therefore, the following two codes produce similar results: \cr
#' \preformatted{
#'  eq( score1=score1, freq1=freq1, score2=score2, freq2=freq2 )
#'  qqplot( expand_freqdist( score1, freq1 )
#'              , expand_freqdist( score2, freq2 ), type="l" )
#' }
#' Note that, since native \code{qqplot} cannot handle case weight,
#' \code{lazy.toold::expand_freqdist} is used to recover the raw data
#' from frequency table.
#'
#'
#' \code{cdf} has priority over \code{freq}.
#'
#'
#'
#' @return a list of the following:\cr
#' ctable: the conversion table consisting of (score, score21, freq1) \cr
#' mands:  the summary stat of the converted score dist: (score21,freq1) \cr
#' newfreq: the frequency distribution of the converted score (score21, freq21)
#' \cr
#' sdist1 and sdist2: input and smoothed score distributions \cr
#' cntr: a list of control parameters used.
#'
#' \code{ctable} shows that
#' test score score1[i] of test1 is equivalent to
#' test score score21[i] of test2.
#'
#' @examples
#'
#' seed <- 1701
#' set.seed(seed)
#'
#' scoredist1 <- gen_test_score( 500, 0,10, beta=c(2,4), plot=1 )
#' scoredist2 <- gen_test_score( 1000, 0,15, beta=c(4,2), plot=1 )
#' reseq <- eq( score1=scoredist1[,1], freq1=scoredist1[,2]
#'            , score2=scoredist2[,1], freq2=scoredist2[,2]
#'            , smooth1=3, smooth2=3, method=2, plot=3 )
#'
#' @export
#'

eq <- function( score1, freq1, cdf1=NULL, score2, freq2, cdf2=NULL
              , lim1=NULL, lim2=NULL
              , smooth1=0, bandwid1=3, smooth2=0, bandwid2=3
              , method=2, interpol_method="linear"
              , title="", nolow=0, round=1, print=1, plot=1 ){
 # Observed Score Equating of Test1 to Test2
 # Shin-ichi Mayekawa
 # 20160722,23,27
 # 20180810dnc,12,14dnc
 # cdf priority: 20230717
 # mandd -> mands: 20230717
 #

 # limits
 if( is.null(lim1) ){
  lim1=c(min(score1),max(score1))
 }
 if( is.null(lim2) ){
  lim2=c(min(score2),max(score2))
 }
 minscore=min(lim1,lim2)
 maxscore=max(lim1,lim2)

 # calculate cdf if not given: when method==1 they are not used.
 if( is.null(cdf1) ){
   cdf1=cumsum(freq1)/sum(freq1)
 }
 else{
   freq1=cumsuminv(cdf1)
 }
 if( is.null(cdf2) ){
   cdf2=cumsum(freq2)/sum(freq2)
 }
 else{
   freq2=cumsuminv(cdf2)
 }
 cdf01=cdf1; cdf02=cdf2

 #
 len1=length(score1)
 n1=sum(freq1)
 len2=length(score2)
 n2=sum(freq2)

 # summary stat
 mands=rbind(t(mands(score1,freq1)[c(1,2,4)]),t(mands(score2,freq2)[c(1,2,4)]))
 rownames(mands)=c("test1","test2")
 colnames(mands)=c("mean","std","n")
 Print(mands)

 if( print >= 1 ){
  cat("\nTest1 Score Distribution\n")
  prmatrix( cbind(score1, freq1, cdf1), rowlab=rep("",len1) )
  cat("\nTest2 Score Distribution\n")
  prmatrix( cbind(score2, freq2, cdf2), rowlab=rep("",len2) )
 }


 if( method == 1 ){

  # linear
  m1=mands[1,1]; s1=mands[1,2]
  m2=mands[2,1]; s2=mands[2,2]

  score210=s2/s1*(score1-m1)+m2
  score21=round( score210 )

  # Print(score1,score21,freq1)

 }
 else if( method == 2  ){

  # equi-percentile

  if( plot >= 3 ){
   title1=paste(title, "\n", "CDF of two tests", sep="")
   plot( score1, cdf1, type="l", ylab=""
         , xlab="", xlim=c(minscore,maxscore), ylim=c(0,1), lty=1 )
   par( new=TRUE )
   plot( score2, cdf2, type="l", xlim=c(minscore,maxscore), ylim=c(0,1)
         , ylab="CDF", xlab="score", main=title1, lty=2 )
  }

  # smooth cdf if required
  if( smooth1 > 0 ){
   cdf1=smoothra( cdf1, bandwid=bandwid1, maxiter=smooth1 )
  }
  if( smooth2 > 0 ){
   cdf2=smoothra( cdf2, bandwid=bandwid2, maxiter=smooth2 )
  }

  if( smooth1 > 0 | smooth2 > 0 ){

   if( print >= 2 ){
    cat("\nTest1 Score Distribution after Smooth\n")
    prmatrix( cbind(score1, freq1, cdf01, cdf1), rowlab=rep("",n1) )
    cat("\nTest2 Score Distribution after Smooth\n")
    prmatrix( cbind(score2, freq2, cdf02, cdf2), rowlab=rep("",n2) )
   }

   if( plot >= 2 ){
    title1=paste(title, "\n", "Smoothed CDF of two tests", sep="")
    plot( score1, cdf1, type="l", ylab=""
          , xlab="", xlim=c(minscore,maxscore), ylim=c(0,1), lty=1 )
    par( new=TRUE )
    plot( score2, cdf2, type="l", xlim=c(minscore,maxscore), ylim=c(0,1)
          , ylab="CDF", xlab="score", main=title1, lty=2 )
   }
  }

  # Print(score1, freq1, cdf1, score2, freq2, cdf2 )

  # equipercentile equating
  score210=interpol( cdf2, score2, cdf1,  method=interpol_method )[,2]
  score21=round( score210 )

 } # end of method == 2


 # nolow
 if( nolow ){
  loc=which(score21 < score1)
  score21[loc]=score1[loc]
 }


 # conversion table
 min21=min(score2,score21); max21=max(score2,score21)
 s21=min21:max21
 f21=numeric(length(s21))
 ctable=cbind(score1,score21,freq1)
 ctable0=cbind(score1,score210,freq1)

 # pdf from cdf
 pdf1=cumsuminv(cdf1)
 pdf2=cumsuminv(cdf2)

 # freq dist of score21
 xt=xtabs( freq1~score21 )
 loc=loceq(min21:max21,score21)
 f21[loc]=xt
 newfreq=cbind(s21,f21)
 colnames(newfreq)=c("score21","freq21")

 # new summary stat
 ms21=unlist(mandd(s21,f21)[c(1,3,4)])
 mands=rbind(mands,ms21)
 rownames(mands)[3]="2from1"

 # Print(ms21)
 # Print(unlist(mandd(score21,pdf1)[c(1,3,4)]))



 if( print ){
  cat("\nSummary Statistics\n")
  print( mands)
  if( print >= 1 ){
   title1=paste(title, "\n", "The Conversion Table (rounded)", sep="")
   cat(title1,"\n")
   prmatrix( ctable, rowlab=rep("",nrow(ctable)) )
   title1=paste(title, "\n", "The Conversion Table", sep="")
   cat(title1,"\n")
   prmatrix( ctable0, rowlab=rep("",nrow(ctable)) )
  }
 }

 if( plot ){
  title1=paste(title, "\n", "Frequency Distribution of the Converted Score"
               , sep="")
  barplot( height=f21, names.arg=s21, space=0, main=title1 )
  title1=paste(title, "\nThe Conversion Table: from Test1 to Test2", sep="")
  plot( score1, score21, type="b", xlim=lim1, ylim=lim2
        , main=title1 )
 }

 cntr=named_list( method, smooth1, smooth2, bandwid1, bandwid2=3
                 , nolow, round )
 sdist1=cbind(score1,freq1, cdf01, cdf1, pdf1 )
 sdist2=cbind(score2,freq2, cdf02, cdf2, pdf2 )
 res=named_list( ctable, ctable0, mands, newfreq, sdist1, sdist2, method, cntr )
 return( res )

} # end of eq








scoredist1 <- gen_test_score( 1500, 0,10, beta=c(3,4), plot=0 )
scoredist2 <- gen_test_score( 1000, 0,15, beta=c(4,2), plot=0 )
reseq <- eq( score1=scoredist1[,1], freq1=scoredist1[,2]
             , score2=scoredist2[,1], freq2=scoredist2[,2]
             , smooth1=3, smooth2=3, method=2, plot=0 )



comments(
'


seed=1701
set.seed(seed)

minscore=0; maxscore=50
n=500

scoredist1=gen_test_score( n, minscore, maxscore, beta=c(9,7), plot=0 )
scoredist2=gen_test_score( n, minscore, maxscore, normal=c(25,8), plot=0 )
reseq=eq( score1=scoredist1[,1], freq1=scoredist1[,2]
          , score2=scoredist2[,1], freq2=scoredist2[,2]
          , smooth1=3, smooth2=3, method=2, plot=3 )








seed=1701
set.seed(seed)


minscore=0
maxscore=100

n=500


scoredist1=gen_test_score( n, minscore, maxscore, beta=c(9,8), plot=0 )
resd1=density(scoredist1[,1], weights=scoredist1[,2]/sum(scoredist1[,2])
              , from=minscore, to=maxscore, n=101)
barplot(height=scoredist1[,2]/sum(scoredist1[,2]), names.arg=scoredist1[,1]
        , space=0
        , main="test1")
lines(resd1, col="blue")
scoredist2=gen_test_score( n, minscore, maxscore, normal=c(50,12), plot=1 )
resd2=density(scoredist2[,1], weights=scoredist2[,2]/sum(scoredist2[,2])
              , from=minscore, to=maxscore)
barplot(height=scoredist2[,2]/sum(scoredist2[,2]), names.arg=scoredist2[,1]
        , space=0
        , main="test2")
lines(resd2, col="blue")

# Print(scoredist1,scoredist2)


reseq=eq( score1=scoredist1[,1], freq1=scoredist1[,2]
    , score2=scoredist2[,1], freq2=scoredist2[,2]
    , smooth1=0, smooth2=3
    , method=2, plot=3 )


ct=reseq$ctable




reseq=eq( score1=scoredist1[,1], freq1=resd1$y
          , score2=scoredist2[,1], freq2=resd2$y
          , smooth1=0, smooth2=0
          , method=2, plot=2, print=9 )
')



