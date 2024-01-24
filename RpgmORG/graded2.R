graded_info <- function( out_obscore, ncat=5, method=1, brk=NULL, print=1 ){
 # calculate the information of the graded score
 # Shin-ichi Mayekawa
 # 20150313,14
 # 20150617
 #
 #

 # How cut and breaks work: example of ncat=3
 #
 #    y=1     brk[1]  <  x  <= brk[2]
 #    y=2     brk[2]  <  x  <= brk[3]
 #    y=3     brk[3]  <  x  <= brk[4]
 #

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
  infox=theta_stat$info
  info=theta_stat$info_LO
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
 if( !is.null(brk) ) ncat=length(brk)+1
 else{
  if( method == 1 ){
   # equal interval
   brk=seq( minscorex,maxscorex,len=ncat+1 )
   brk[1]=brk[1]-(scorex[2]-scorex[1])/2
   brk[ncat+1]=brk[ncat+1]+(scorex[ncat0]-scorex[ncat0-1])/2
  }
  else{
   # brk is a set of quantiles at equally spaced points in (0,1)
   brk=wquantil( scorex, Px, seq(0,1,len=ncat+1) )
   brk[1]=minscorex-0.5; brk[ncat+1]=maxscorex+0.5
  }
 }
 #ncat=length(brk)-1

 if( print ){
  cat( "\n\n Calculation of the Information Function\n"
       , " associated with the Graded Observed Score\n")
  Print( minscorex, maxscorex, ncat0, ncat, method )
  cat( " Break Points for conversion from X to Y\n" )
  cat("\n ", brk,"\n")
 }

 # categorized score
 cc=unclass( cut( scorex, brk, include.lowest=TRUE ) )
 #Print( scorex,cc,table(cc) )

 # score distribution of y, the graded score from x
 scorey=1:ncat
 Py_t=matrix(0,ncat,npoint)
 rownames(Py_t)=scorey; colnames(Py_t)=thname
 for( k in 1:ncat ){
  loc=which( cc == k )
  Py_t[k,]=colSums(Px_t[loc,,drop=0])
 }

 if( print ){
  cat( "\n Conditional Distribution of Y given theta\n" )
  Print(Py_t, fmt=".2")
 }


 # TRF test response function etc
 TRFy=c( scorey%*%Py_t )
 stdy_t=c( sqrt( (scorey^2)%*%Py_t - TRFy^2 ) )

 # slope of TRFy by difference: numerical differentiation
 dTRFy=TRFy
 for( k in 2:(npoint-1) ){
  dTRFy[k]=(TRFy[k+1]-TRFy[k-1])/(theta[k+1]-theta[k-1])
 }
 dTRFy[1]=(TRFy[2]-TRFy[1])/(theta[2]-theta[1])
 dTRFy[npoint]=
  (TRFy[npoint]-TRFy[npoint-1])/(theta[npoint]-theta[npoint-1])

 # information function associated with y
 infoy=(dTRFy^2)/(stdy_t^2)
 rsqrinfoy=sqrt(1/infoy)

 if( print ){
  cat( "\n Information Functions\n" )
  Print( theta, info, infox, infoy, fmt="6.3" )
 }


 # plot info
 temp=cbind(info,infox,infoy)
 matplot( theta,temp, type = "l", ylab="information"
        , lty=1, col=1:5, cex=1, lwd=1
        , main="Comparison of Information Functions" )
 legend( 3,max(info), legend=c("info","infoX","infoY")
       , cex=1, lwd=1, col=1:5, merge=TRUE )

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


 matplot( theta, t(Pt_y), type = "l", lty=1, col=1:5, cex=1, lwd=1
        , main="Posterior Distribution of Theta given Y", ylab="prob")
 legend( -3,max(Pt_y), legend=as.character(scorey)
         , cex=1, lwd=1, col=1:5, merge=TRUE )

} # end of graded_info



out_obscore=obscore( paramS2, weight=weightS21, npoints=21, print=0 )
graded_info( out_obscore, ncat=5, method=1 )
graded_info( out_obscore, ncat=9, method=1 )






