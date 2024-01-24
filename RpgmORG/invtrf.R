#' Inverse Function of trf
#'
#' @param param Item parameter data frame
#' @param x The value of x, the observed test score, or NULL
#' @param xmin The minimum value of x, or NULL
#' @param xmax The maximum value of x, or NULL
#' @param by_x  The increment of x values
#' @param weight The weight data frame or the vector of item weight, or NULL
#' @param theta The value of theta for method=0, or NULL
#' @param thmin The minimum value of theta for method=0, or NULL
#' @param thmax The maximum value of theta for method=0, or NULL
#' @param npoints The # of theta points for method=0
#' @param method = 1 to use native \code{uniroot} function (slow) \cr
#' = 0 to use interpolation with option \code{lazy.tools::interpol} (fast).
#' @param interpol The interpolation method. See \code{lazy.tools::interpol}.
#' @param print = 2 to print the result
#' @param plot = 1 to plot the result
#'
#' @details
#' The test response function, trf, is calculated as:
#' \preformatted{
#'  trf(theta) = sum_{j=1}^n w_j x sum_{k=0}^{ncat[i]-1} v_{kj} P_{kj}(theta)
#'  where
#'  w_j is the item weight for item j, (default=1)
#'  v_{kj} is the item category weight for the k-th category of item j.
#'  The default value of v_{kj} = k, k=0, 1, 2, ... ncat[i]-1
#' }
#' When values \code{x=NULL} the value of \code{x} will be calculated
#' using the weight information.
#'
#' If x is too small or too large, -9999 or 9999 will be returned.
#'
#' @return
#' A list of \cr
#' \code{x} The value of x.\cr
#' \code{theta} The value of theta. \cr
#' \code{locminmax} A list containing the locations of -9999 or 9999.
#'
#' @examples
#' res <- invtrf( paramS1, method=1, print=2 )
#' res <- invtrf( paramS1, method=0, print=2 )
#' res <- invtrf( paramS1, weight=weightS12, method=1, print=2 )
#' res <- invtrf( paramS1, weight=weightS12, method=0, print=2 )
#'
#' res <- invtrf( paramS1, by_x=0.5, print=2 )
#' res <- invtrf( paramS1, weight=c(1,1,2,2), print=2 )
#'
#' @export
#'

invtrf <- function( param, x=NULL, xmin=NULL, xmax=NULL, by_x=1
                  , weight=NULL, theta=NULL, thmin=-5, thmax=5, npoints=121
                  , method=0, interpol="spline"
                  , print=0, plot=0 ){
 # calculation of theta from x
 # Shin-ichi Mayekawa
 # 20230714
 # comment: 20230716cot
 #

 # constants
 nitem=nrow(param)
 ncat=param$ncat
 iname=param$name
 itype=param$type


 # weight to be used
 # If weight is a valid weight data frame, pass it to irf as it is.
 # else if weight is a scalar or vector, convert it to a valid weight df,
 # else if weight is NULL, pass it to irf as it is.
 if( is.data.frame(weight) ){
  if( !("w" %in% colnames(weight)) ){
   cat("\nerror1:(invtrf) weight data frame does not have 'w' variable.\n")
   return( NULL )
  }
 }
 else{
  # weight is a vector, scalar or NULL
  if( is.null(weight) ) weight=1
   if( is.vector(weight) ) weight=matrix(weight,nitem)
   w=weight; rownames(w)=iname; colnames(w)="w"
   v=matrix(0:(max(ncat)-1),nitem,max(ncat),byrow=1)
   rownames(v)=iname; colnames(v)=paste("v",0:(max(ncat)-1),sep="")
   for( i in 1:nitem ){
    if( ncat[i]+1 <= ncol(v) ) v[i,(ncat[i]+1):ncol(v)]=NA
   }
   weight=cbind(ncat,w,v); rownames(weight)=iname
   colnames(weight)=c("ncat","w",colnames(v))
   weight=data.frame(name=iname,type=itype, weight)
 }

 # max score
 w=weight$w
 v=weight[,regexpr("^v", colnames(weight)) > 0, drop=0]
 maxscore_i=as.matrix(apply(v,1,max,na.rm=1),nitems)
 rownames(maxscore_i)=iname
 maxscore=sum(w*maxscore_i)
 minscore_i=as.matrix(apply(v,1,min,na.rm=1),nitems)
 rownames(minscore_i)=iname
 minscore=min(w*minscore_i)

 # x range
 if( is.null(x) ){
  if( is.null(xmin) ) xmin=minscore
  if( is.null(xmax) ) xmax=maxscore
  x=seq(xmin,xmax,by_x)
  if( x[length(x)] != xmax ) x=c(x,xmax)
 }
 x=c(x)
 nx=length(x)


 if( print ){
    cat("\nInverse Function of trf at X\n\n")
    cat("# of x points =", nx, "\n")
    Print(nitem,minscore,maxscore)
    Print(method, interpol)
 }


 if( method == 1 ){

  # use uniroot
  # function whose root is sought.
  invtrf <- function( theta ){
     xx - irf( param, theta=theta, weight=weight, print=0 )$TRF
  } # end of invtrf
  theta=numeric(nx)

  xxmin=irf( param, theta=thmin, weight=weight, print=0 )$TRF
  xxmax=irf( param, theta=thmax, weight=weight, print=0 )$TRF

  for( i in 1:nx ){
   xx=x[i]
   if( xx <= xxmin ) th=-9999
   else if( xx >= xxmax ) th=9999
   else{
    th=uniroot(invtrf, c(thmin, thmax))$root
   }
   theta[i]=th
  }

 } # end of method=1
 else{

  # use interpol
  # theta to be used
  if( is.null(theta) ) th=seq(thmin,thmax,length.out=npoints)
  else th=theta
  trf=c( irf( param, theta=th, weight=weight, print=0 )$TRF )

  mintrf=min(trf); maxtrf=max(trf)
  locmin=which(x < mintrf); locmax=which( x > maxtrf)

  theta=numeric(nx)
  theta[locmin]=-9999; theta[locmax]=9999
  xx=x[-c(locmin,locmax)]
 #Print(xx,x,theta,th,trf)
  th0=interpol( trf, th, xx )[,2]
  theta[-c(locmin,locmax)]=th0

 } # end of method=2


 locminmax=list(locmin=which(theta==-9999), locmax=which(theta==9999))



 if( print >= 2 ){
  Print(x, theta, fmt="5.1 12.4")
 }

 if( plot ){
  plot( x[-unlist(locminmax)], theta[-unlist(locminmax)]
        , main=paste("method =", method))
 }
 return( named_list( x, theta, locminmax ) )

} # end of invtrf



res=invtrf( paramS1, npoints=21, method=1, by_x=.1 )


comments(
'


res=invtrf( paramS1, weight=weightS1, npoints=21 )
res=invtrf( paramS1, npoints=21 )
res=invtrf( paramS1, weight=2, npoints=21 )
res=invtrf( paramS1, weight=c(1,1,2,2), npoints=21 )



'
)
