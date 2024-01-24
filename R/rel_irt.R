#' Calculation of Test Reliability and SEM under IRT
#'
#' @param param Item Parameter Data Frame
#' @param weight Weight data frame
#' @param maxscore The maximum score to be used to calculate var(X) etc.
#' @param npoints # of discrete points for theta
#' @param thmin Minimum value of discrete thata value
#' @param thmax Maximum value of discrete thata value
#' @param thdist Type of theta distribution \cr
#' = 0 to use uniform,  \cr = 1 to use N(0,1)
#' @param print = 1 to print result
#' @param plot = 1 or 2 to plot result
#'
#' @details
#' Test reliability of a text Z is defined as:
#' \preformatted{
#'  rel(Z) = var(T) / var(Z) = 1 - var(E)/var(Z)
#' }
#' where T and E are, resp, the true score and the error score of Z. \cr
#' For observed score: X = trf(theta) + E   and  E(E|theta) = 0,
#' \preformatted{
#' var(E) = E(var(E|theta)) = E(var(X|theta))
#' var(X) = E(var(X|theta)) + var(E(X|theta)) = var(E) + var(trf)
#' var(T) = var(trf)
#' }
#'
#' For thetahat based score: Y = thetahat = theta + error and E(error|theta)=0,
#' \preformatted{
#' var(error) = E(var(error|theta)) = var(Y|theta) = 1/info(theta)
#' var(Y) = E(var(Y|theta)) + var(E(Y|theta)) = E(1/info(theta)) + var(theta)
#' var(theta) = 1 or var(uniform)
#' }
#' \cr
#' The standard error of measurement, SEM, is defined as, \cr
#' For X
#' \preformatted{
#'  SEMx(X) = sqrt( var(E|X) )  or  SEMx(theta) = sqrt( var(E|trf(theta)) )
#' }
#' For Y,
#' \preformatted{
#'  SEMt(Y) = sqrt( var(error|theta) ) or
#'  SEMt(X) = sqrt( var(error|trf(theta)) )
#' }
#' Their averages are denoted, resp, as \cr
#' aSEM or aSEM_theta.
#' \cr
#' Note that, SEMx(X) and SEMx(theta) are convex and
#' SEMt(Y) and SEMt(theta) are concave. \cr
#'
#' Note that the value depends on the range and the shape
#' of the theta distribution.
#' \cr
#' Try changing npoints, thmin and thmax.
#'
#' @examples
#' res <- rel_irt( paramA1, thmin=-4, thmax=4, plot=2 )
#' res <- rel_irt( paramA1, thmin=-4, thmax=4, thdist=0, plot=2 )
#'
#' @return
#' A list of: \cr\cr
#' thmin, thmax, npoints, thdist \cr
#' SigmaX2 Observed Score Variance \cr
#' SigmaT2 True Score Variance \cr
#' SigmaT22 Shoulld be the same as above \cr
#' SigmaE2 Error Score Variance \cr
#' aSEM Average Standard Error of Measurement \cr
#' rel The reliability coefficient \cr
#' \cr
#' SigmaX_theta2 Variance of Thetahat \cr
#' SigmaT_theta2 Variance of Theta\cr
#' SigmaE_theta2 Average Variance of Thetahat \cr
#' aSEM_theta Average Standard Error of Measurement of Thetahat \cr
#' rel_theta The reliability coefficient of Thetahat based Score \cr
#' \cr
#' SigmaX_xs2 Variance of X in [0,maxscore] \cr
#' SigmaT_xs2 Variance of T in [0,maxscore] \cr
#' SigmaE_xs2 Average Variance of E in [0,maxscore]  \cr
#' aSEM_xs Average Standard Error of Measurement of X in [0,maxscore] \cr
#' rel_xs The reliability coefficient \cr
#' \cr
#' SigmaX_ys2 Variance of Thetahat in [0,maxscore] \cr
#' SigmaT_ys2 Variance of Theta in [0,maxscore] \cr
#' SigmaE_ys2 Average Variance of Thetahat given theta in [0,maxscore]  \cr
#' aSEM_ys Average Standard Error of Measurement of Thetahat in [0,maxscore]
#'  \cr
#' rel_ys The reliability coefficient \cr
#' \cr
#'  SEM_x      SEM of X at True Score \cr
#'  SEM_t      SEM of X and Y(thetahat) at Theta \cr
#'  SEM_x100   Rescaled SEM of X and Y(thetahat) at Rescaled True Score \cr
#'  SEM_t100   Rescaled SEM of X and Y(thetahat) at theta \cr
#'
#'
#' @references
#' Lord, F. M. (1980).
#' Applications of Item Response Theory To Practical Testing Problems.
#' Lawrence Erlbaum Associates. Mahwah, New Jergey
#'
#'
#' Samejima, F. (1994) Estimation of Reliability Coefficients Using
#' the Test Information Function and Its Modifications.
#' Applied Psychological Measurement. 18, 229-244.
#'
#' Toyoda. H. (1989) The Methods for Estimating the Reliability Coefficient
#' under Item Response Model. (in Japanese)
#' Japanese Journal of Educational Psychology. vol. 37, 283-285.
#'
#' @export
#'

rel_irt <- function( param, weight=NULL, maxscore=100, npoints=51
                     , thmin=-4, thmax=4, thdist=1, print=1, plot=1 ){
 # calculation of test reliability and average SEM under irt
 # Shin-ichi Mayekawa
 # mostly, except from obscore
 # 20171205
 # locally optimal weight: 20171208,09
 # above moved to info_func: 20171209
 # 0-100 score: 20171212dnc, 1214
 # roxgen comment: 20180112,0430cot
 # Y -> theta and SE -> SEM in plot: 20180430cot
 # plot irf: 20230511
 # plot SEM vs Original True Score: 20230511
 # roxygen2 comments and TCC -> trf: 20230714
 #

 # names
 pdfname=as.character(substitute(param))
 wdfname=as.character(substitute(weight))

 param=checkparam( param, "ALL", "obscore" )
 if( is.null(param) ){
  cat("\n\nInput Item Parameter ERROR.\n\n")
  return()
 }

 # constants
 nitems=nrow(param)
 iname=param$name
 nc=ncol(param)
 ncat=param$ncat
 type=param$type

 # item and category weight
 if( is.null(weight) ){
  # natural weight
  w=matrix(1,nitems); rownames(w)=iname; colnames(w)="w"
  v=matrix(0:(max(ncat)-1),nitems,max(ncat),byrow=1)
  rownames(v)=iname; colnames(v)=paste("v",0:(max(ncat)-1),sep="")
  for( i in 1:nitems ){
   if( ncat[i]+1 <= ncol(v) ) v[i,(ncat[i]+1):ncol(v)]=NA
  }
  weight=data.frame(iname, type, ncat,w,v); rownames(weight)=iname
  colnames(weight)=c("name","type","ncat","w",colnames(v))
 }
 w=weight$w
 v=weight[,5:ncol(weight), drop=0]
 if( length(w) != nitems | nrow(v) != nitems |  ncol(v) != max(ncat)  ){
  cat("\nerror1(obscore): Weight data frame does not conform. \n")
  return()
 }

 maxscore_i=apply(weight[,5:ncol(weight)],1,max,na.rm=1)
 maxscore_t=sum(w*maxscore_i)
 minscore_i=apply(weight[,5:ncol(weight)],1,min,na.rm=1)
 minscore_t=sum(w*minscore_i)

 # generate thata and prior theta dist
 theta=seq(thmin,thmax,length.out=npoints)
 thname=format(theta,digits=2)
 if( thdist == 0 ) Pt=matrix(1/npoints,npoints)
 else if( thdist == 1 ){
  Pt=exp(-0.5*theta^2);
  Pt=Pt/sum(Pt)
 }
 thmean=sum(theta*Pt)
 thvar=sum(theta^2*Pt)-thmean^2
 thdisttype=c("uniform","normal(0,1)")
 Pt_fr=formatC( Pt, format="f", digits=2, width=4 )
 names(Pt_fr)=thname

 if( print ){
  cat("\n\nTest Reliability and Average SEM\n")
  cat(" parameter data frame name =", pdfname
      , ",  item category weight data frame name =", wdfname,"\n\n")
  cat(" # of item parameters =",nitems,"\n")
  cat(" range of theta = [",thmin,",",thmax,"] with", npoints
      ,"discrete points\n")
  cat( " theta distribution type = ", thdisttype[thdist+1],"\n")
  cat( "   mean theta =", thmean, ",   var theta =", thvar, "\n")
  cat("\nMaximum Observed Score =", maxscore_t, "\n")

 }



 # calculate icrf and trf=cond. mean of X given theta
 temp=irf( param, theta, weight, print=0, debug=0, plot=0 )
 icrf=temp$ICRF
 irf=temp$IRF
 trf=temp$TRF
 fromP=temp$fromP
 toP=temp$toP
 vecv=temp$vecv
 rm(temp)


 # first derivative of icrf
 temp=dirf( param, theta, weight=weight, print=0 )
 dicrf=temp$dICRF
 slope_trf=temp$dTRF
 dirf=temp$dIRF
 rm(temp)


 # conditional variance of the weighted total score X;
 # This should be equal to the one calculated using Px_t.
 varx_t=matrix(0,npoints,1)
 varuvw_t=matrix(0,npoints,nitems)
 varuv_t=matrix(0,npoints,nitems)
 for( k in 1:npoints ){
  Pk=icrf[k,,drop=0]
  dum=0
  for( j in 1:nitems ){
   Pjk=Pk[fromP[j]:toP[j]]
   DD=diag(Pjk)-outer(Pjk,Pjk,"*")
   vj=matrix(vecv[fromP[j]:toP[j]],ncat[j])
   vDDv=t(vj)%*%DD%*%vj
   varuv_t[k,j]=vDDv
   vjwj=vj*w[j]
   vDDv=t(vjwj)%*%DD%*%vjwj
   dum=dum+vDDv
   varuvw_t[k,j]=vDDv
  }
  varx_t[k]=dum
 }
 stdx_t=sqrt(varx_t)


 # variance of X: var(X) = E( var(X|theta) ) + var( E(X|theta) )
 #                         Average Error Var + True Score Var
 SigmaX2=sum(Pt*varx_t) + ( sum(Pt*trf^2)-sum(Pt*trf)^2 )
 # reliability etc
 SigmaE2=sum(Pt*varx_t)
 SigmaT2=SigmaX2-SigmaE2
 SigmaT22=sum(trf^2*Pt)-sum(trf*Pt)^2
 # average standard error of measurement
 aSEM=sqrt(SigmaE2)
 rel=SigmaT2/SigmaX2

 if( print ){
  cat("\nVariances and the Test Reliability of the Observed Score\n")
  Print( SigmaX2, SigmaT2, SigmaT22, SigmaE2, fmt="8.3")
  Print( aSEM, rel, fmt="8.3")
 }



 # stdx_t as the function of True Score
 # interpol (trf,stdx_t) at trf=score
 score=0:maxscore_t
 stdx_x=interpol( trf, stdx_t, score )[,2]






 # information function with the locally optimal category weights
 info_LO=rowSums( (dicrf^2)/icrf )
 var_thetahat=1/info_LO
 SE_thetahat=sqrt(var_thetahat)


 # variance of theta
 SigmaT_theta2=thvar
 # expectation of conditional variance of thetahat given theta
 SigmaE_theta2=sum(Pt*var_thetahat)
 # var(thetahat) = E( var(thetahat|theta) ) + var( E(thetahat|theta) )
 SigmaX_theta2=SigmaE_theta2 + SigmaT_theta2

 # average standard error of estimation of theta
 aSEM_theta=sqrt(SigmaE_theta2)
 rel_theta=SigmaT_theta2/SigmaX_theta2

 if( print ){
  cat("\nVariances and the Test Reliability of the Theta-hat based Score\n")
  Print( SigmaX_theta2, SigmaT_theta2, SigmaE_theta2, fmt="8.3")
  Print( aSEM_theta, rel_theta, fmt="8.3")
 }


 # Rescaling: range conversion
 # x: (trf(thmin), trf(thmax)) -> (0, maxscore)
 # y: (thmin, thmax) -> (0, maxscore)

 # for X
 score100=0:maxscore
 ddx=maxscore/(trf[npoints]-trf[1])
 trf100=(trf-trf[1])*ddx
 SE_x=ddx*stdx_t
 SE_xs=interpol( trf100, SE_x, score100 )[,2]
 aSEM_xs=aSEM*ddx

 # for Y (thetahat)
 # Note that, since  this is a linear transformation of theta,
 # the shape of se does not change.
 ddy=maxscore/(thmax-thmin)
 y=(theta-thmin)*ddy
 SE_y=ddy*SE_thetahat
 SE_ys=interpol( y, SE_y, score100 )[,2]
 aSEM_ys=aSEM_theta*ddy

 # variance etc
 SigmaX_xs2=SigmaX2*ddx^2
 SigmaT_xs2=SigmaT2*ddx^2
 SigmaE_xs2=SigmaE2*ddx^2
 rel_xs=SigmaT_xs2/SigmaX_xs2

 SigmaX_ys2=SigmaX_theta2*ddy^2
 SigmaT_ys2=SigmaT_theta2*ddy^2
 SigmaE_ys2=SigmaE_theta2*ddy^2
 rel_ys=SigmaT_ys2/SigmaX_ys2

 if( print ){
  cat("\nVariances and the Test Reliability of X and Y scaled to [0,"
      , maxscore,"]\n", sep="")
  Print( SigmaX_xs2, SigmaT_xs2, SigmaE_xs2, aSEM_xs, rel_xs, fmt="6.2")
  Print( SigmaX_ys2, SigmaT_ys2, SigmaE_ys2, aSEM_ys, rel_ys, fmt="6.2")
 }



 # Print(score,SE_xs,SE_ys)

 if( plot ){


  # plot irf
  matplot( theta, irf, type="l"
           , ylim=c(0,max(irf)), ylab="irf", main="Item Response Functions")
  leg=paste("average SEM ="
            , formatC(aSEM, digits=2, width=6, format="f"), sep="")
  leg=paste( leg, ",  rel = "
             , formatC(rel, digits=2, width=4, format="f" ), sep="" )
  legend("topleft", legend=leg, col = 1)

  # plot trf +- SEM vs theta

  upper=trf+stdx_t
  lower=trf-stdx_t
  title="Test Response Function with SEM"
  matplot( theta, cbind(trf,upper,lower), type="l", lty=c(1,2,2)
     , ylim=c(0,maxscore_t), col=1, main=title
     , xlab="theta", ylab="X", lw=c(3,1,1), yaxp=c(0,maxscore_t,maxscore_t) )
  leg=paste("trf:  average SEM ="
            , formatC(aSEM, digits=2, width=6, format="f"), sep="")
  leg=paste( leg, ",  rel = "
             , formatC(rel, digits=2, width=4, format="f" ), sep="" )
  legend("topleft", legend=leg, col = 1, lty = 1:2)

  # plot trf +- SEM vs trf
  upper=trf+stdx_t
  lower=trf-stdx_t
  upper=interpol( trf, upper, score )[,2]
  lower=interpol( trf, lower, score )[,2]
  tt=interpol( trf, trf, score )[,2]
  title="True Score with SEM"
  matplot( tt, cbind(tt,upper,lower), type="l", lty=c(1,2,2)
       , ylim=c(0,maxscore_t), col=1, main=title
       , xlab="T", ylab="X", lw=c(3,1,1)
       , yaxp=c(0,maxscore_t,maxscore_t), xaxp=c(0,maxscore_t,maxscore_t) )
  leg=paste("T:  average SEM ="
            , formatC(aSEM, digits=2, width=6, format="f"), sep="")
  leg=paste( leg, ",  rel = "
             , formatC(rel, digits=2, width=4, format="f" ), sep="" )
  legend("topleft", legend=leg, col = 1, lty = 1:2)


  # plot of SEM vs theta
  plot( theta, stdx_t, type="l", ylim=c(0,max(stdx_t)*1.1)
        , main="Plot of SEM as the function of Theta")
  leg=paste("SEM:  average SEM ="
            , formatC(aSEM, digits=2, width=6, format="f"), sep="")
  leg=paste( leg, ",  rel = "
             , formatC(rel, digits=2, width=4, format="f" ), sep="" )
  legend("bottom", legend=leg, col = 1, lty = 1:2)

  # plot of SEM vs true score
  plot( score, stdx_x, type="l", ylim=c(0,max(stdx_x)*1.1)
        , xlim=c(0,maxscore_t), xaxp=c(0,maxscore_t,maxscore_t)
        , main="Plot of SEM as the function of True Score")
  leg=paste("SEM:  average SEM ="
            , formatC(aSEM, digits=2, width=6, format="f"), sep="")
  leg=paste( leg, ",  rel = "
             , formatC(rel, digits=2, width=4, format="f" ), sep="" )
  legend("bottom", legend=leg, col = 1, lty = 1:2)


  if( plot >= 2 ){

    # plot of RESCALED  SEM vs theta and score
    leg=paste("X:  average SEM ="
              , formatC(aSEM_xs, digits=2, width=6, format="f"), sep="")
    leg=paste( leg, ",  rel = "
               , formatC(rel_xs, digits=2, width=4, format="f" ), sep="" )
    leg=c(leg, paste("theta:  average SEM ="
                  , formatC(aSEM_ys,digits=2, width=6, format="f")
                  , ",  rel = ", formatC(rel_ys,digits=2, width=4, format="f")
                  , sep="") )

    # plot SEM vs theta
    title1=paste("Plot of Rescaled SEMs as the function of Theta:"
                 , " theta distribution =", thdisttype[thdist+1], sep="" )
    matplot( theta, cbind(SE_x,SE_y), type="l", main=title1, ylab="SEM", col=1
             , ylim=c(0,max(c(SE_x,SE_y))) )
    legend("top", legend=leg, col = 1, lty = 1:2)

    # plot of SEM vs score
    title2=paste( "Plot of Rescaled SEMs as the function of Score: "
                  , " theta distribution =", thdisttype[thdist+1], sep="" )
    matplot( score100, cbind(SE_xs,SE_ys), type="l"
             , main=title2, ylab="SEM", col=1
             , ylim=c(0,max(c(SE_x,SE_y))) , xaxp=c(0,maxscore,maxscore) )
    legend("top", legend=leg, col = 1, lty = 1:2)

  }

 }


  SEM_x=cbind(score,stdx_x)
  SEM_t=cbind(theta, stdx_t, SE_thetahat)
  SEM_x100=cbind(score100,SE_xs,SE_ys)
  SEM_t100=cbind(theta ,SE_x,SE_y)

  res=named_list( thmin, thmax, npoints, thdist, maxscore
                  , SigmaX2, SigmaT2, SigmaT22, SigmaE2, aSEM, rel
                  , SigmaX_theta2, SigmaT_theta2, SigmaE_theta2
                  , aSEM_theta, rel_theta
                  , SigmaX_xs2, SigmaT_xs2, SigmaE_xs2, aSEM_xs, rel_xs
                  , SigmaX_ys2, SigmaT_ys2, SigmaE_ys2, aSEM_ys, rel_ys
                  , SEM_x, SEM_t, SEM_x100, SEM_t100
  )

 return( res )

} # end of rel_irt

