#' Find a Transformation g of the Observed Score X such that
#' Y=g(X) has a Flat Standard Error of Measurement.
#'
#' @param out_obscore The result of obscore function \cr
#' This has the priority over the set of the arguments of obscore function.
#' @param sigma The standard error of the transformed score
#' @param t_ty The interval of the true score to be used in the calculation
#'
#' @param param Item Parameter Data Frame for obscore
#' @param weight Weight data frame for obscore
#' @param npoints # of discrete points for theta for obscore
#' @param thmin Minimum value of discrete thata value for obscore
#' @param thmax Maximum value of discrete thata value for obscore
#' @param thdist Type of theta distribution for obscore \cr
#' = 0 to use uniform,  = 1 to use N(0,1)
#' @param alpha small prob for quantile and confidence interval for obscore
#' @param compress = 1 to remove zero-probability weighted total observed
#'      scores for obscore
#'
#' @param print > 1 to print result
#' @param plot > 1 to plot result
#' @param debug = 1 to print intemediate result
#'
#' @return A list of the following: \cr
#' t: The value of the true score \cr
#' stdx_t: SEM of X at T \cr
#' u: The transformed true score: Y=g(X) and u=g(t) \cr
#' gdash: The derivative of g \cr
#' stdy_u: SEM of Y at u \cr
#' lengtht: length of t \cr
#' sigma: New SEMvalue specified \cr\cr
#' out_obscore: The output from the obscore function.
#'
#' @details
#'
#' Let \code{stdx(t)} be the standard error of measurement of X at t. \cr
#' This can be calculated as \code{stdx}_t by the obscore function. \cr
#' The standard deviation of \code{Y=g(X)} at t can be approximated by \cr
#'  \code{ g-dash(t)*stdx(t) }\cr
#' and we want it to be a constant (sigma). \cr
#' Therefore, \cr
#'  \code{ g-dash(t) = sigma / stdx(t)) } \cr
#' and the \code{g} function can be recovered by integrating the above
#' \code{g-dash}.
#'
#' @examples
#'
#' # tiny set of binary items
#' param=paramB1
#' maxscore=sum(param$ncat-1)
#' param$p1=1
#' out_obscore <- obscore( param )
#' res=flatten_SEM( out_obscore, sigma=1, t_by=0.5, plot=1, print=1 )
#'
#' # binary and polytomous items
#' res2=flatten_SEM( param=paramS1, sigma=1, t_by=0.1, plot=1, print=1 )
#'
#'
#' @export
#'

flatten_SEM <- function( out_obscore=NULL, sigma=1, t_by=0.5
                         , param=NULL, weight=NULL
                         , npoints=31, thmin=-4, thmax=4, thdist=1, alpha=0.1
                         , compress=0, print=1, plot=0, debug=0){
 # Find a transformation of the observed score to the one with flat SEM.
 # Shin-ichi Mayekawa
 # 20170115,20dnc
 # break point of x to u: 20170222
 # comments: 20170307dnc
 #

 # Y = g(X)

 get_gdash <- function( t ){
  # derivative of u=g(t): g'(t)=gd0=sigma/stdx_TRF
  res=interpol( TRF, gd0, t )[,2]
  return( res )
 } # gdash


 if( is.null(out_obscore) & is.null(param) ){
  cat("\n\nerror1:(flatten_SEM) out_obscore or param must be given.\n\n")
  return( NULL )
 }

 if( is.null(out_obscore) ){
  out_obscore=obscore( param=param, weight=weight
                       , npoints=npoints, thmin=thmin, thmax=thmax
                       , thdist=thdist, alpha=alpha, compress=compress
                       , print=max(0,print-2), plot=max(0,plot-2)
                       , debug=debug )
 }

 # Here, obscore is already run
 theta_stat=out_obscore$theta_stat
 theta=theta_stat$theta
 TRF=theta_stat$TRF
 stdx_TRF=theta_stat$stdx_t
 minscore=out_obscore$minscore_t
 maxscore=out_obscore$maxscore_t
 # Print(theta,TRF,stdx_TRF)
 #
 pdfname=out_obscore$pdfname
 wdfname=out_obscore$wdfname
 npoints=out_obscore$npoints
 thmin=out_obscore$thmin
 thmax=out_obscore$thmax
 thdist=out_obscore$thdist
 nitems=out_obscore$nitems

 # from_obscore=named_list(pdfname,wdfname,npoints,thmin,thmax,thdist,nitems)

 # Equally spaced true score and SEM values
 t=seq( minscore, maxscore, t_by )
 stdx_t=interpol( TRF, stdx_TRF, t )[,2]
 lengtht=length(t)

 # reciprocal of stdx_TRF to be integrated (assuming len(TRF) > lent(t))
 gd0=sigma/stdx_TRF

 # find the integral of gdash: (t,u) gives the transformation u=g(t).
 u=Integrate( get_gdash, t, 0 )

 # sem of Y=g(X) at t
 gdash=get_gdash( t )
 stdy_u=gdash*stdx_t

 # discretize t to get u
 u05=seq(min(floor(u))+.5,max(floor(u))+.5)
 brk=interpol( u,t, u05 )[,2]
 brk_x2u=cbind(u05,brk)
 colnames(brk_x2u)=c("u","brk")

# plot( t, u, type="l" )
# points( brk, u05 )

 if( plot ){
  plot( t, u, type="l"
        , main=paste("Transformation g of the Observed Score:  sigma ="
             , sigma, sep=""), sub="Y = g(X)  and  U = g(T)" )
  points( brk, u05 )
  plot( u, stdy_u, type="l"
        , main=paste("SEM of Transformed Score:  sigma ="
                , sigma,sep=""), sub="Y = g(X)  and  U = g(T)")
  plot( t, stdx_t, type="l"
        , main="SEM of the Observed Score" )
  plot( TRF, stdx_TRF, type="l"
        , main="SEM of the Observed Score from obscore" )
 }

 if( print ){
  cat("\nTransformation of the Observed Score\n\n")
  Print( sigma, t_by, lengtht )
  Print( nitems, minscore, maxscore)
  Print( thdist, thmin, thmax, npoints )
  Print( t, stdx_t, u, gdash, stdy_u )
  Print( brk_x2u )
 }

 res=named_list( t, stdx_t, u, gdash, stdy_u, lengtht, sigma, brk_x2u
                , out_obscore )
 return( res )

} # flatten_SEM







comments('


param=paramS1
maxscore=sum(param$ncat-1)
param$p1=1
out_obscore <- obscore( param )
res=flatten_SEM( out_obscore, sigma=1, t_by=0.1, plot=1, print=1 )




res2=flatten_SEM( param=pp, sigma=1, t_by=0.5, plot=1, print=1, npoints=151 )


res3=graded_info( res2$out_obscore, brk=res2$brk_x2u[,2], plot=1 )


')









