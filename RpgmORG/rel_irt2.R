#' Calculation of Test Reliability and Average SEM under IRT
#'
#' @param param Item Parameter Data Frame
#' @param weight Weight data frame
#' @param npoints # of discrete points for theta
#' @param thmin Minimum value of discrete thata value
#' @param thmax Maximum value of discrete thata value
#' @param thdist Type of theta distribution \cr
#' = 0 to use uniform,  \cr = 1 to use N(0,1)
#' @param print = 1 to print result
#'
#' @details
#' For observed score X = TCC(theta) + E: \cr
#' var(E) = E(var(E|theta)) = E(var(X|theta)) \cr
#' var(X) = E(var(X|theta)) + var(E(X|theta)) \cr
#'        = var(E) + var(TCC) \cr
#' var(T) = var(TCC) \cr
#' \cr
#' For thetahat based score \cr
#' var(thetahat|theta) = 1/info \cr
#' var(thetahat) = E(var(thetahat|theta)) + var(E(thetahat|theta)) \cr
#'               = E(var(thetahat|theta)) + \cr
#' var(theta) = known \cr
#'
#' \cr
#' Note that the value depends on the range and the shape
#' of the theta distribution.
#' \cr
#' See the following examples.
#'
#' @examples
#' res <- rel_irt( paramB1, thmin=-2, thmax=2 )
#' res <- rel_irt( paramB1, thmin=-3, thmax=3 )
#' res <- rel_irt( paramB1, thmin=-4, thmax=4 )
#' res <- rel_irt( paramB1, thmin=-5, thmax=5 )
#'
#' res <- rel_irt( paramB1, thmin=-5, thmax=5, thdist=0 )
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
#' v_LO
#' Locally Optimal Weight scaled so that v[0,j]=0 and v[1,1]=1.
#'
#' @export
#'

rel_irt <- function( param, weight=NULL
                     , npoints=51, thmin=-4, thmax=4, thdist=1, print=1 ){
 # calculation of test reliability and average SEM under irt
 # Shin-ichi Mayekawa
 # mostly, except from obscore
 # 20171205
 # locally optimal weight: 20171208,09
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

 # locally optimal item weight
 w_LO=dirf/varuv_t

 # locally optimal categorie weight: scaled so that v[0,j]=0 and v[1,1]=1.
 # v_LO=dicrf/icrf
 # for( j in 1:nitems ){
 # v_LO[,fromP[j]:toP[j]]=v_LO[,fromP[j]:toP[j]] - v_LO[,fromP[j]]
 # }
 # v1=v_LO[,2]
 # for( j in 1:nitems ){
 #  v_LO[,fromP[j]:toP[j]]=v_LO[,fromP[j]:toP[j]] / v1
 # }

 # locally optimal categorie weight: scaled so that v[0,j]=0 and v[1,1]=1.
 dd=dicrf/icrf
 v_LO=mapply( function(x,y){dd[,x:y]-dd[,x]}, fromP, toP, SIMPLIFY=FALSE )
 v_LO=matrix(unlist(v_LO),npoints)
 v1=v_LO[,2]
 v_LO=mapply( function(x,y){v_LO[,x:y]/v1}, fromP, toP, SIMPLIFY=FALSE  )
 v_LO=matrix(unlist(v_LO),npoints)

 # information function with the locally optimal category weights
 info_LO=rowSums( (dicrf^2)/icrf )
 rsqrinfo_LO=sqrt(1/info_LO)
 rinfo_LO=1/info_LO

 # variance of theta
 SigmaT_theta2=thvar
 # expectation of conditional variance of thetahat given theta
 SigmaE_theta2=sum(Pt*rinfo_LO)
 # var(thetahat) = E( var(thetahat|theta) ) + var( E(thetahat|theta) )
 SigmaX_theta2=SigmaE_theta2 + SigmaT_theta2

 # average standard error of estimation of theta
 aSEM_theta=sqrt(SigmaE_theta2)
 rel_theta=1-SigmaE_theta2/SigmaT_theta2

 if( print ){
  cat("\nVariances and the Test Reliability of the Theta-hat based Score\n")
  Print( SigmaX_theta2, SigmaT_theta2, SigmaE_theta2, fmt="8.3")
  Print( aSEM_theta, rel_theta, fmt="8.3")
 }




 res=named_list( thmin, thmax, npoints, thdist
                 , SigmaX2, SigmaT2, SigmaT22, SigmaE2, aSEM, rel
                 , SigmaX_theta2, SigmaT_theta2, SigmaE_theta2
                 , aSEM_theta, rel_theta, w_LO, v_LO)

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

'
)


param=paramB1[1:3,]
param$type="B3"
#param$p1=c(0.5, 1.0, 1.5)
param$p2=1
param$p3=c(0.0, 0.001, 0.8)
thmin=-3; thmax=3; npoint=21
theta=seq(thmin,thmax,length=npoint)
res <- rel_irt( param, npoints=npoint, thmin=thmin, thmax=thmax )
v_LO=res$v_LO
v_LO1=v_LO[,seq(2,2*nrow(param),2)]
matplot(theta,v_LO1,type="l")
Print(res$w_LO,v_LO1)

matplot(theta,res$w_LO,type="l")


