#' Calculation of the Information Function associated with the
#' Graded Observed Score.
#'
#' @param out_obscore  Output from obscore function
#' @param ncat # of categories to be used
#' @param method = 1 to equal interval on X \cr
#' = 2 to equal interval on percent
#' @param brk Break points to be used or NULL.
#' This has priority over ncat.
#' @param scorey  The value of Y   or  0 to length(brk)-1.
#' @param print = 1 to print result
#' @param plot = 1 to plot information functions\cr
#' = 2 to plot posterior distribution of theta given Y.
#'
#'
#' @return  A list of: \cr
#'  theta Discrete theta points defined in obscore function. \cr
#'
#'  info Information function (LO) defined in obscore function. \cr
#'
#'  infoX Information function associated with X defined in obscore function.
#' \cr
#'  infoY Information function associated with Y \cr
#'
#'  TRFy_t_t
#'  stdy_t
#'
#'  Py_t Distribution of Y given theta \cr
#'
#'  Pt_y Distribution of theta given Y \cr
#'
#'  meant_y Posterior mean of theta given Y \cr
#'
#'  stdt_y Posterior std of theta given Y \cr
#'
#' @details
#' The graded score, Y, 0 <= Y <= ncat-1,  will be calculated on the basis of
#' the (weighted) observed score X. \cr
#' Then, the probability distribution of Y given theta will be calculated
#' by summing the probability distribution of X given theta. \cr
#' Finally, the information function associated with the graded score Y will be
#' calculated as the ratio of the slope of TRF of Y squared to the
#' conditional variance of Y given theta. \cr
#' The slope of TRF will calculated numerically.
#'
#'
#' @examples
#' # Define the observed score X using the category and item weights given
#' # in weightsS21, and calculate the score distribution etc.
#' out_obscore <- obscore( paramS2, weight=weightS21, npoints=21, print=0 )
#' # On the basis of the observed score X calculated above using weightS21,
#' # categorize X into ncat categories to create new score Y.
#' res <- graded_info( out_obscore, ncat=5, method=1, plot=1 )
#' res <- graded_info( out_obscore, ncat=9, method=1, plot=1 )
#'
#' @export
#'

graded_info <- function( out_obscore, ncat=5, method=1
                         , brk=NULL, scorey=NULL, print=1, plot=0 ){
 # calculate the information of the graded score
 # Shin-ichi Mayekawa
 # 20150313,14
 # 20150617
 # plot std vs theta: 20170223
 # 1:ncat -> 0:(ncat-1)
 # bugfix: 20170227
 # comments: 20170307
 # bugfix when brk is given: 20170308
 # 100*SEM, new numerical deriv: 20170309
 #
 #

 # How cut and breaks work: example of ncat=3
 #
 #    y=0     brk[1]  <=  x  <= brk[2]   1st level of cut
 #    y=1     brk[2]  <   x  <= brk[3]   2nd level of cut
 #    y=2     brk[3]  <   x  <= brk[4]   3rd level of cut
 #

 inv_tcc <- function( x, theta, tcc ){
  # inverse of tcc at x
  res=interpol( tcc, theta, x )[,2]
  return( res )
 } # inv_tcc


 if( !is.null( out_obscore ) ){
  # from obscore output
  # obs_stat:
  #  score  Px  meant_x stdt_x  qth_L  qth_U  qth_wid2 ci_L  ci_U   ci_wid2
  obs_stat=out_obscore$obs_stat
  # theta_stat:
  # theta Pt TRF slope_TRF stdx_t info r_sqr_info info_LO r_sqr_info_LO
  # qt_L qt_U   poststd
  theta_stat=out_obscore$theta_stat
  theta=theta_stat$theta
  Pt=theta_stat$Pt
  infoX=theta_stat$info
  info=theta_stat$info_LO
  stdx_t=theta_stat$stdx_t
  tcc=theta_stat$TRF
  scorex=obs_stat$score
  Px=obs_stat$Px
  Px_t=out_obscore$Px_t
  thname=formatC(theta,digits=2, format="f")
 }
 else{
  cat("\n\n error1:(graded) out_obscore must be specified. \n\n\n")
  return()
 }


 npoint=ncol(Px_t)
 ncat0=length(scorex)
 minscorex=scorex[1]; maxscorex=scorex[ncat0]

 # break points
 if( !is.null(brk) ){
  ncat=length(brk)-1
  method=0
  xmidp=brk2midp( brk )
 }
 else{
  if( method == 1 ){
   # equal interval
   xmidp=seq( minscorex,maxscorex,len=ncat )
   # brk[1]=brk[1]-(scorex[2]-scorex[1])/2
   # brk[ncat+1]=brk[ncat+1]+(scorex[ncat0]-scorex[ncat0-1])/2
   brk=midp2brk( xmidp )
  }
  else{
   # brk is a set of quantiles at equally spaced points in (0,1)
   brk=wquantile( scorex, Px, seq(0,1,len=ncat+1) )
   brk[1]=minscorex-0.5; brk[ncat+1]=maxscorex+0.5
   xmidp=brk2midp( brk )
  }
 }

 if( is.null(scorey) )  scorey=0:(ncat-1)
 minscorey=scorey[1]; maxscorey=scorey[ncat]


 brk_theta=inv_tcc( brk[-c(1,ncat+1)], theta, tcc )
 brk_theta=c(NA,brk_theta,NA)

 if( print ){
  cat( "\n\n Calculation of the Information Function\n"
       , " associated with the Graded Observed Score\n")
  Print( minscorex, maxscorex, ncat0 )
  Print( minscorey, maxscorey, ncat, method )
  cat( "\n Break Points for conversion from X to Y\n" )
  cat("\n ", brk,"\n")
  Print(scorex,xmidp,brk,brk_theta, scorey)
 }

 # categorized score
 cc=unclass( cut( scorex, brk, include.lowest=TRUE ) )
 Print( scorex,cc,table(cc) )

 # score distribution of y, the graded score from x
 Py_t=matrix(0,ncat,npoint)
 rownames(Py_t)=scorey; colnames(Py_t)=thname
 for( k in 1:ncat ){
  loc=which( cc == k )
  Py_t[k,]=colSums(Px_t[loc,,drop=0])
 }
 sumP=colSums(Py_t)

 if( print ){
  cat( "\n Conditional Distribution of Y given theta\n" )
  Print(Py_t, fmt=".2")
  Print("  Checking...", min(sumP),max(sumP))
 }


 # TRF test response function etc
 TRFy_t=c( scorey%*%Py_t )
 stdy_t=c( sqrt( (scorey^2)%*%Py_t - TRFy_t^2 ) )

 comments('
 # slope of TRFy_t by difference: numerical differentiation
 dTRFy_t=TRFy_t
 for( k in 2:(npoint-1) ){
  dTRFy_t[k]=(TRFy_t[k+1]-TRFy_t[k-1])/(theta[k+1]-theta[k-1])
 }
 dTRFy_t[1]=(TRFy_t[2]-TRFy_t[1])/(theta[2]-theta[1])
 dTRFy_t[npoint]=
  (TRFy_t[npoint]-TRFy_t[npoint-1])/(theta[npoint]-theta[npoint-1])
 ')

 dTRFy_t=TRFy_t
 for( k in 2:(npoint-1) ){
  dTRFy_t[k]=( (TRFy_t[k]-TRFy_t[k-1])/(theta[k]-theta[k-1]) +
              (TRFy_t[k+1]-TRFy_t[k])/(theta[k+1]-theta[k]) )/2
 }
 dTRFy_t[1]=(TRFy_t[2]-TRFy_t[1])/(theta[2]-theta[1])
 dTRFy_t[npoint]=
  (TRFy_t[npoint]-TRFy_t[npoint-1])/(theta[npoint]-theta[npoint-1])




 # information function associated with y
 infoY=(dTRFy_t^2)/(stdy_t^2)
 rsqrinfoY=sqrt(1/infoY)

 if( print ){
  cat( "\n Information Functions\n" )
  Print( theta, info, infoX, infoY, fmt="6.3" )
 }


 if( plot ){
  title=paste("Plot of E(Y|theta) and SEM of Y:  ncat = ",ncat
  , " (method = ", method, ")", sep="")
  title2=paste("score range:  X(",minscorex,",",maxscorex,"),  Y("
               ,minscorey,",",maxscorey,")", sep="")
  temp=cbind(TRFy_t,stdy_t)
  matplot( theta, temp, type="l", main=title, sub=title2)
  legend( 2,max(temp)-0.2, legend=c("E(Y|theta)","SE(theta)")
          , cex=1, lwd=1, col=1:2, merge=TRUE )

  title=paste("Plot of SEM of Y and E(Y|theta):  ncat = ",ncat
              , " (method = ", method, ")", sep="")
  title2=paste("score range:  X(",minscorex,",",maxscorex,"),  Y("
               ,minscorey,",",maxscorey,")", sep="")
  plot( TRFy_t, stdy_t, type="l", main=title, sub=title2)
  legend( 2,max(temp)-0.2, legend=c("E(Y|theta)","SE(theta)")
          , cex=1, lwd=1, col=1:2, merge=TRUE )

  title=paste("Plot of SEM of X and E(X|theta)"
              , " (method = ", method, ")", sep="")
  title2=paste("score range:  X(",minscorex,",",maxscorex,"),  Y("
               ,minscorey,",",maxscorey,")", sep="")
  plot( tcc, stdx_t, type="l", main=title, sub=title2)
  legend( 2,max(temp)-0.2, legend=c("E(Y|theta)","SE(theta)")
          , cex=1, lwd=1, col=1:2, merge=TRUE )

  # plot info
  temp=cbind(info,infoX,infoY)
  title=paste("Comparison of Information Functions:  ncat = ", ncat
              , " (method = ", method, ")", sep="")
  matplot( theta,temp, type = "l", ylab="information"
           , lty=1, col=1:5, cex=1, lwd=1
           , main=title, sub=title2 )
  legend( 2,max(info), legend=c("info","infoX","infoY")
          , cex=1, lwd=1, col=1:5, merge=TRUE )

  temp=cbind(stdx_t,stdy_t)
  title=paste("Comparison of SEMs:  ncat = ",ncat
              , " (method = ", method, ")", sep="")
  matplot( theta,temp, type = "l", ylab="std of score"
           , lty=1, col=2:5, cex=1, lwd=1
           , main=title, sub=title2 )
  legend( 2,max(temp), legend=c("X","Y")
          , cex=1, lwd=1, col=2:5, merge=TRUE )

  temp=cbind(100*stdx_t/(maxscorex-minscorex),100*stdy_t/maxscorey)
  title=paste("Comparison of relative SEMs: 100*SEM/maxscore:  ncat = ", ncat
              , " (method = ", method, ")", sep="")
  matplot( theta,temp, type = "l", ylab="std of score"
           , lty=1, col=2:5, cex=1, lwd=1
           , main=title
           , sub=title2 )
  legend( 2,max(temp), legend=c("X","Y")
          , cex=1, lwd=1, col=2:5, merge=TRUE )

 }

 # posterior distribution of theta given y
 Pt_y=Py_t*matrix(Pt,ncat,npoint,byrow=1)
 Pt_y=(1/rowSums(Pt_y))*Pt_y
 meant_y=Pt_y%*%theta
 stdt_y=sqrt( Pt_y%*%(theta^2) - meant_y^2 )

 if( print ){
  cat( "\n Conditional Distribution of theta given Y\n" )
  Print(Pt_y,fmt=".2")
  Print(scorey, meant_y, stdt_y)
 }

 if( plot >= 2 ){
  matplot( theta, t(Pt_y), type = "l", lty=1, col=1:5, cex=1, lwd=1
           , main="Posterior Distribution of Theta given Y", ylab="prob")
  legend( -3,max(Pt_y), legend=as.character(scorey)
          , cex=1, lwd=1, col=1:5, merge=TRUE )
 }


 res=named_list( theta, info, infoX, infoY, TRFy_t, stdy_t
           , dTRFy_t, Py_t, Pt_y, meant_y, stdt_y )

 return( res )

} # end of graded_info




#res3=graded_info( res2$out_obscore, brk=res2$brk_x2u[,2], method=2, plot=1 )



comments(
'


pp=rbind(paramB1,paramB1,paramB1)
pp$name=1:nrow(pp)

out_obscore=obscore( pp, weight=NULL, npoints=121, print=0, plot=3 )

res=graded_info( out_obscore, ncat=3, method=1, plot=1 )
res=graded_info( out_obscore, ncat=5, method=1, plot=1 )
res=graded_info( out_obscore, ncat=9, method=1, plot=1 )





# param data frame EQ5:
np=100
p1=1
p2=0
p3=0
name=paste("Q",1:np,sep="")
type="B"
ncat=2
param_EQ=data.frame( name,type,ncat,p1,p2,p3, stringsAsFactors = 0 )

out_obscore=obscore( param_EQ, weight=NULL, npoints=121, print=0, plot=3
, thmin=-2, thmax=2 )

res=graded_info( out_obscore, ncat=3, method=1, plot=1 )
res=graded_info( out_obscore, ncat=5, method=1, plot=1 )
res=graded_info( out_obscore, ncat=9, method=1, plot=1 )













')




