#' Calculation of Test Reliability and Average SEM under IRT
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
#' @param plot = 1 to plot result
#'
#' @details
#' For observed score: X = TCC(theta) + E: \cr
#' var(E) = E(var(E|theta)) = E(var(X|theta)) \cr
#' var(X) = E(var(X|theta)) + var(E(X|theta)) \cr
#'        = var(E) + var(TCC) \cr
#' var(T) = var(TCC) \cr
#' \cr
#' For thetahat based score: Y = theta + error \cr
#' vae(error) = E(var(error|theta)) = var(thetahat|theta) = 1/info \cr
#' var(thetahat) = E(var(thetahat|theta)) + var(E(thetahat|theta)) \cr
#'               = E(1/info) + var(theta)  \cr
#' var(theta) = 1 \cr
#'
#' \cr
#' Note that the value depends on the range and the shape
#' of the theta distribution.
#' \cr
#' Try changing npoints, thmin and thmax.
#'
#' @examples
#' res <- rel_irt( paramB1, thmin=-3, thmax=3 )
#'
#' res <- rel_irt( paramB1, thmin=-3, thmax=3, thdist=0 )
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
#' aSEM_xs Average Standard Error of Measurement ofX in [0,maxscore] \cr
#' rel_xs The reliability coefficient \cr
#' \cr
#' SigmaX_ys2 Variance of Thetahat in [0,maxscore] \cr
#' SigmaT_ys2 Variance of Theta in [0,maxscore] \cr
#' SigmaE_ys2 Average Variance of Thetahat given theta in [0,maxscore]  \cr
#' aSEM_ys Average Standard Error of Measurement of Thetahat in [0,maxscore]
#'  \cr
#' rel_ys The reliability coefficient \cr
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

 }



 # calculate icrf and trf=cond. mean of X given theta
 temp=irf( param, theta, weight, print=0, debug=0, plot=0 )
 icrf=temp$ICRF
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


 # range conversion
 # x: (tcc(thmin), tcc(thmax)) -> (0, maxscore)
 # y: (thmin, thmax) -> (0, maxscore)

 score=0:maxscore
 ddx=maxscore/(trf[npoints]-trf[1])
 x=(trf-trf[1])*ddx
 SE_x=ddx*stdx_t
 SE_xs=interpol( x, SE_x, score )[,2]
 aSEM_xs=aSEM*ddx

 ddy=maxscore/(thmax-thmin)
 y=(theta-thmin)*ddy
 SE_y=ddy*SE_thetahat
 SE_ys=interpol( y, SE_y, score )[,2]
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
  cat("\nVariances and the Test Reliability of the Two Scores in [0,"
      , maxscore,"]\n", sep="")
  Print( SigmaX_xs2, SigmaT_xs2, SigmaE_xs2, aSEM_xs, rel_xs, fmt="6.2")
  Print( SigmaX_ys2, SigmaT_ys2, SigmaE_ys2, aSEM_ys, rel_ys, fmt="6.2")
 }



 # Print(score,SE_xs,SE_ys)

 if( plot ){
  # plot
  matplot( theta, icrf[,seq(2,2*nitems,2)], type="l"
           , ylim=c(0,1), ylab="icrf", main="Item Response Functions")

  upper=trf+stdx_t
  lower=trf-stdx_t
  title="Test Response Function with SEM"
  matplot( theta, cbind(trf,upper,lower), type="l", lty=c(1,2,2)
           , ylim=c(0,maxscore_t), col=1, main=title
           , xlab="", ylab="", lw=c(3,1,1), yaxp=c(0,maxscore_t,maxscore_t) )



  # plot
  leg=paste("X:  average SEM ="
            , formatC(aSEM_xs, digits=2, width=6, format="f"), sep="")
  leg=paste( leg, ",  rel = "
             , formatC(rel_xs, digits=2, width=4, format="f" ), sep="" )
  leg=c(leg, paste("theta:  average SEM ="
                   , formatC(aSEM_ys,digits=2, width=6, format="f")
        , ",  rel = ", formatC(rel_ys,digits=2, width=4, format="f")
        , sep="") )
  title1="Plot of SEMs as the function of Theta"
  title2=paste( "Plot of SEMs as the function of Score: "
                , " theta distribution =", thdisttype[thdist+1], sep="" )
  matplot( theta, cbind(SE_x,SE_y), type="l", main=title2, ylab="SEM", col=1
           , ylim=c(0,max(c(SE_x,SE_y))) )
  legend("top", legend=leg, col = 1, lty = 1:2)
  matplot( score, cbind(SE_xs,SE_ys), type="l", main=title2, ylab="SEM", col=1
           , ylim=c(0,max(c(SE_x,SE_y))) )
  legend("top", legend=leg, col = 1, lty = 1:2)


  res=named_list( thmin, thmax, npoints, thdist, maxscore
                  , SigmaX2, SigmaT2, SigmaT22, SigmaE2, aSEM, rel
                  , SigmaX_theta2, SigmaT_theta2, SigmaE_theta2
                  , aSEM_theta, rel_theta
                  , SigmaX_xs2, SigmaT_xs2, SigmaE_xs2, aSEM_xs, rel_xs
                  , SigmaX_ys2, SigmaT_ys2, SigmaE_ys2, aSEM_ys, rel_ys
  )

 }

 return( res )

} # end of rel_irt





comments(
 '
 res <- rel_irt( paramB1, thdist=1 )
 res <- rel_irt( paramB1, thdist=0 )


 res <- rel_irt( paramB1, thmin=-2, thmax=2 )
 res <- rel_irt( paramB1, thmin=-3, thmax=3 )
 res <- rel_irt( paramB1, thmin=-4, thmax=4 )
 res <- rel_irt( paramB1, thmin=-5, thmax=5 )
 res <- rel_irt( paramB1, thmin=-6, thmax=6 )





 param=paramB1[1:3,]
 thmin=-3; thmax=3; npoint=51
 theta=seq(thmin,thmax,length=npoint)
 res <- rel_irt( param, npoints=npoint, thmin=thmin, thmax=thmax )
 v_LO=res$v_LO
 v_LO1=v_LO[,seq(2,2*nrow(param),2)]
 matplot(theta,v_LO1,type="l")
 matplot(theta,res$w_LO,type="l")
 Print(res$w_LO,v_LO1)

 param=paramB1[1:3,]
 param$p1=c(0.5, 1.0, 1.5)
 # param$p2=1
 thmin=-3; thmax=3; npoint=11
 theta=seq(thmin,thmax,length=npoint)
 res <- rel_irt( param, npoints=npoint, thmin=thmin, thmax=thmax )
 v_LO=res$v_LO
 v_LO1=v_LO[,seq(2,2*nrow(param),2)]
 matplot(theta,v_LO1,type="l")
 matplot(theta,res$w_LO,type="l")
 Print(res$w_LO,v_LO1)







 param=paramB1[1:3,]
 param$type="B3"
 #param$p1=c(0.5, 1.0, 1.5)
 param$p2=1
 param$p3=c(0.0, 0.001, 0.8)
 thmin=-3; thmax=3; npoint=21
 theta=seq(thmin,thmax,length=npoint)
 res <- rel_irt( param, npoints=npoint, thmin=thmin, thmax=thmax )



 '
)


res <- rel_irt( param, thmin=-3, thmax=3 )
# res <- rel_irt( paramB1, thmin=-4, thmax=4 )
