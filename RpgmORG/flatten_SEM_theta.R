#' Find a Transformation g of the Thetahat Based Observed Score X such that
#' Y=g(X) has a Flat Standard Error of Measurement.
#'
#' @param out_obscore The result of obscore function \cr
#' This has the priority over the set of the arguments of obscore function.
#' @param sigma The standard error of the transformed score
#' @param by_s The interval for continuous S.
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
#' theta: The value of theta \cr
#' stdx_theta: SEM of X at theta \cr
#' s: The transformed true score: Y=g(X) and s=g(t) \cr
#' gdash: The derivative of g \cr
#' stdy_s: SEM of Y at s \cr
#' lengtht: length of t \cr
#' sigma: New SEM value specified \cr
#' brk_x2u: Break points of X to create S.\cr
#' brk_x2uc: Break points of X to create almost condinuous S.\cr\cr
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
#' \code{g-dash}. \cr
#' This \code{g} is the vaiance-stabilizing transformation.
#'
#' Notes: \cr
#'  Recommended to use \code{npoints=151, thmin=-4, thmax=4} or larger
#'  for obscore function.
#'
#' @examples
#'
#' # tiny set of binary items
#' param=paramB1
#' maxscore=sum(param$ncat-1)
#' param$p1=1
#' out_obscore <- obscore( param )
#' res=flatten_SEM_theta( out_obscore, sigma=1, plot=1, print=1 )
#'
#' # binary and polytomous items
#' res2=flatten_SEM_theta( param=paramS1, sigma=1, plot=1, print=1 )
#'
#'
#' @export
#'

flatten_SEM_theta <- function( out_obscore=NULL, sigma=1, by_s=0.1
                         , param=NULL, weight=NULL
                         , npoints=131, thmin=-4, thmax=4, thdist=1, alpha=0.1
                         , compress=0, print=1, plot=0, debug=0){
 # Find a transformation of the theta hat to the one with flat SEM.
 # Shin-ichi Mayekawa
 # 20170115,20dnc
 # break point of x to s: 20170222
 # comments: 20170307dnc
 # get rid of equally spaced t: 20170308
 # use of midp2brk: 20170308
 # use min and max of TCC, output ncat: 20170310
 # pint ncat: 20170314
 # u -> s: 20170417
 # modified version of flatten_SEM
 # thetahat: 20171214
 #

 # Y = g(X)

 get_gdash <- function( t ){
  # derivative of s=g(t): g'(t)=gd0=sigma/stdthetahat(t)
  res=interpol( theta, gd0, t )[,2]
  return( res )
 } # gdash


 if( is.null(out_obscore) & is.null(param) ){
  cat(
   "\n\nerror1:(flatten_SEM_theta) out_obscore or param must be given.\n\n")
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
 stdthetahat=1/sqrt( theta_stat$info_LO ) # This is SEM of thetahat
 minscore=out_obscore$minscore_t
 maxscore=out_obscore$maxscore_t
 # Print(theta,TRF,stdthetahat)
 #
 pdfname=out_obscore$pdfname
 wdfname=out_obscore$wdfname
 npoints=out_obscore$npoints
 thmin=out_obscore$thmin
 thmax=out_obscore$thmax
 thdist=out_obscore$thdist
 nitems=out_obscore$nitems

 # from_obscore=named_list(pdfname,wdfname,npoints,thmin,thmax,thdist,nitems)

 # # Equally spaced true score and SEM values
 # t=seq( minscore, maxscore, t_by )
 # stdx_t=interpol( TRF, stdthetahat, t )[,2]
 th=theta
 stdx_theta=stdthetahat
 lengtht=length(th)

 # reciprocal of stdthetahat to be integrated (assuming len(TRF) > lent(t))
 gd0=sigma/stdthetahat

 # find the integral of gdash: (t,s) gives the transformation s=g(t).
 s=Integrate( get_gdash, th, thmin )

 # sem of Y=g(X) at t
 gdash=get_gdash( th )
 stdy_s=gdash*stdx_theta

 # discretize t to get s
 umidp=seq(0,ceiling(max(s)), 1)
 ubrk=midp2brk( umidp )
 brk=interpol( s,th, ubrk )[,2]
 brk[length(brk)]=max( brk[length(brk)], thmax)
 # Print(umidp,ubrk,brk)
 brk_x2u=cbind(c(NA,umidp),brk)
 colnames(brk_x2u)=c("s","brk")
 ncat=length(umidp)

 # continuous
 umidpc=seq(0,ceiling(max(s)), by_s)
 ubrkc=midp2brk( umidpc )
 brkc=interpol( s,th, ubrkc )[,2]
 brkc[length(brkc)]=max( brkc[length(brkc)], thmax)
 # Print(umidpc,ubrkc,brkc)
 brk_x2uc=cbind(c(NA,umidpc),brkc)
 colnames(brk_x2uc)=c("uc","brkc")
 ncatc=length(umidpc)

 # plot( t, s, type="l" )
 # points( brk, u05 )

 if( plot ){
  plot( th, s, type="l"
        , main=paste("Transformation g of the Observed Score:  sigma ="
             , sigma, ",  ncat =", ncat,sep="")
        , sub="Y = g(X)  and  S = g(T)" )
  points( brk, ubrk )
  # plot( t, stdx_theta, type="l"
  #      , main="SEM of the Observed Score" )
 }

 if( print ){
  cat("\nTransformation of the Observed Score\n\n")
  Print( sigma, lengtht, by_s, ncat, ncatc )
  Print( nitems, minscore, maxscore)
  Print( thdist, thmin, thmax, npoints )
  Print( th, stdx_theta, s, gdash, stdy_s )
  Print( brk_x2u )
 }

 res=named_list( th, stdx_theta, s, gdash, stdy_s, lengtht, sigma, brk_x2u
               , brk_x2uc, ncat, ncatc, out_obscore )
 return( res )

} # flatten_SEM_theta







comments('


param=pp
maxscore=sum(param$ncat-1)
param$p1=1
out_obscore <- obscore( param, thmin=-3, thmax=3, npoints=31 )
res=flatten_SEM_theta( out_obscore, sigma=1, plot=1, print=1 )




res2=flatten_SEM_theta( param=pp, sigma=1, plot=1, print=1, npoints=151 )


res3=graded_info( res2$out_obscore, brk=res2$brk_x2u[,2], plot=1 )


')









