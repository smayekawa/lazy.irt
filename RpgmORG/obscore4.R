obscore <- function( param, weight=NULL
                   , npoints=31, thmin=-4, thmax=4, thdist=1, alpha=0.1
                   , compress=0, print=1, plot=0, debug=0 ){
 # calculation of observed score distribution
 # Shin-ichi Mayekawa
 # 120210,11,12
 # 120213
 # negetive and non integer scores: 120216,20,21
 # checkparam and natural weight: 120302
 # bugfix: 20150313
 # info_LOW: 20150322 IC, 20150330
 # matplot bugfix: 20150613
 # legend: 20150616
 # bugfix for LOW: 20150617
 #
 #
 # Args:
 #    param    parameter data frame name
 #    weight   item and category data frame name
 #             see the descriptions in read.param or read.weight.
 #             If weight==NULL a set of narural weights will be used.
 #
 #    npoints  # of theta points
 #    thmin    min value of theta
 #    thmax    max value of theta
 #    thdist   = 0 to use uniform,  = 1 to use N(0,1)
 #    alpha    small prob for quantile and confidence interval
 #    compress = 1 to remove zero-probability weighted total observed scores
 #    print    = 1 to verbose
 #    plot     = 1 to produce several plot
 #
 #
 #
 # Values:
 #
 #  list( theta_stat, obs_stat, Px_t, Pt_x,  etcetc )
 #
 #   where
 #
 #  theta_stat as data frame
 #    theta       theta points
 #    Pt          prior probability distribution of theta
 #    TRF         test response function
 #    slope_TRF   slope of TRF
 #    stdx_t      standard deviation of X (observed score) given theta
 #    info        information function defind as (slope_TRF)^2 / (stdx_t)^2
 #    info_LOW    information function with locally optimal item weight
 #    info_LO     information function with locally optimal category weight
 #    qt_L        upper quantile of X given theta
 #    qt_U        lower quantile of X given theta
 #    poststd     posterior std of theta given X
 #                   as a function of posterior mean
 #
 #  obs_stat as data frame
 #    score      domain of X (observed score)
 #    Px         marginal probability  of X
 #    post_mean  posterior mean of theta given X
 #    post_std   posterior standard deviation of theta given X
 #    ci_L       upper limit of conficence interval for theta given X
 #    ci_U       lower limit of conficence interval for theta given X
 #    ci_hwid    half the width of CI
 #
 #
 #  Px_t           score x theta  conditional prob of X given theta
 #  Pt_x           score x theta  conditional prob of theta given X
 #  Px             score x 1      marginal dist of X
 #
 #
 # Needs:
 #  irf,  icrfB, icrfG,  icrfPN
 #  dirf,  dicrfB, dicrfG,  dicrfPN
 #  sumsmnw,  sumsmnw12
 #  checkparam
 #


 fuzzy.fmt <- function( A, fmt, fuzz=-1, fzchar=".", revr=0, revc=0, print=0 ){
  # format a matrix according to fmt and fuzz
  # Shin-ichi Mayekawa
  # 120211
  # avoid list: 120307
  #
  #  Args:
  #   A         a matrix to be printed
  #   fmt       fmt to be used in sprintf
  #   fuzz      Those elements of A blow fuzz will be printed using fzchar
  #   revr = 1  to reverse the rows
  #   revc = 1  to reverse the columns
  #   print = 1 to print the result
  #

  if( is.list(A) ){
   cat("\n\nerror:(fuzzy.fmt) Input cannot be a list.\n\n")
   return()
  }

  A[which(abs(A)<fuzz)]=NA
  Af=sprintf(fmt,A)
  Af[grep("NA *",Af)]=fzchar
  if( is.matrix(A) ){
   Af=matrix( Af,nrow(A),ncol(A) )
   rownames(Af)=rownames(A); colnames(Af)=colnames(A)
  }
  if( revr ){
   Afr=Af; Afr=Af[nrow(A):1,,drop=0]
   Af=Afr
  }
  if( revc ){
   Afr=Af; Afr=Af[,ncol(A):1,drop=0]
   Af=Afr
  }

  if( print > 0 ) {
   print( Af, quote=0 )
  }
  return( Af )

 } # end of fuzzy.fmt





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
 v=weight[,5:ncol(weight)]

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


 # calculate score vector using real cat. weight and dummy probs
 Pk=matrix(NA,max(ncat),nitems)
 for( j in 1:nitems ){ Pk[1:ncat[j],j]=1/ncat[j] }
 temp=sumsmnw( Pk, t(v), w, compress=compress
               , print=0, plot=0, debug=0 )
 score=temp[,1,drop=0];  nscore_t=nrow(score);
 minscore_t=score[1]; maxscore_t=score[nscore_t]
 temp=temp[,2]
 locmiss=which(temp == 0); nmiss=length(locmiss)
 rm(temp)


 if( debug > 0 ) Print(iname,nitems,nc,ncat,w, "/",maxscore_i, maxscore_t)

 if( print > 0 ){
  cat("\n\nCalculation of")
  cat(" the Marginal Distribution of the Weighted Observed Score \n")
  cat(" and the Posterior Mean/Std of Theta\n\n")
  param=cbind(param,maxscore_i)
  weight=cbind(weight,maxscore_i)
  cat(" parameter data frame name =", pdfname
      , ",  item category weight data frame name =", wdfname,"\n\n")
  cat(" # of item parameters =",nitems,"\n")
  cat(" range of theta = [",thmin,",",thmax,"] with", npoints
      ,"discrete points\n")
  cat(" alpha for quantile and CI =", alpha,"\n")
  cat(" range of observed score = [",minscore_t,",", maxscore_t,"]\n")
  cat("  (total of", nscore_t," discrete score points with ",nmiss
      ," zero probability points)\n")
  cat("\n item parameters\n")
  print( param )
  cat("\n item and item category weight \n")
  print( weight )
  cat("\n")
 }


 # calculate icrf and trf=cond. mean of X given theta
 temp=irf( param, theta, weight, print=0, debug=0, plot=0 )
 icrf=temp$ICRF
 trf=temp$TRF
 fromP=temp$fromP
 toP=temp$toP
 vecv=temp$vecv
 rm(temp)

 if( debug > 0 ) Print(icrf, fromP,toP)



 # slope of trf by difference: numerical differentiation
 slopetrf=trf;
 for( k in 2:(npoints-1) ){
  slopetrf[k]=(trf[k+1]-trf[k-1])/(theta[k+1]-theta[k-1])
 }
 slopetrf[1]=(trf[2]-trf[1])/(theta[2]-theta[1])
 slopetrf[npoints]=
  (trf[npoints]-trf[npoints-1])/(theta[npoints]-theta[npoints-1])
 if( debug > 0 ) Print(theta,slopetrf)

 # first derivative of icrf
 temp=dirf( param, theta, weight=weight, print=0 )
 dicrf=temp$dICRF
 slope_trf=temp$dTRF
 dirf=temp$dIRF
 rm(temp)

 if( print >= 1 ){
  c=cor(slopetrf,slope_trf, use="pairwise.complete.obs")
  cat("\n Correlation coefficient between two slopes ="
      , format(c,digits=6),"\n")
 }
 if( debug ) Print(slopetrf,slope_trf, digits=2)

 rm(slopetrf)

 # conditional variance of the weighted total score X;
 # This should be equal to the one calculated using Px_t.
 stdx_t=matrix(0,npoints,1)
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
  stdx_t[k]=dum
 }
 stdx_t=sqrt(stdx_t);

 # information function associated with the obaserved weighted score
 info=(slope_trf^2)/(stdx_t^2)
 rsqrinfo=sqrt(1/info)
 if( debug > 0 ) Print(theta, stdx_t, info, rsqrinfo)

 # information function with the locally optimal category weights
 info_LO=rowSums( (dicrf^2)/icrf )
 rsqrinfo_LO=sqrt(1/info_LO)

 # information function with the locally optimal item weights
 info_LOW=rowSums( (dirf^2)/varuv_t )
 rsqrinfo_LOW=sqrt(info_LOW)


 ################# This is the core of this program #######################
 # calculate cond dist of x given theta as the sum of scored multinomials #
 ##########################################################################
 Px_t=matrix(0,maxscore_t+1,npoints)
 Px_t=NULL
 for( k in 1:npoints ){
  Pk=matrix(NA,max(ncat),nitems)
  for( j in 1:nitems ){
   Pk[1:ncat[j],j]=t( icrf[k,fromP[j]:toP[j]] )
  }
  if( debug > 0 ) Print(k,Pk,t(v),digits=2)
  Px_t=cbind( Px_t
              , sumsmnw( Pk, t(v), w, compress=compress
                         , print=0, plot=0, debug=0 )[,2] )
 }
 rownames(Px_t)=score; colnames(Px_t)=thname
 if( debug ) Print(minscore_t, maxscore_t, nscore_t, score)

 if( any( abs(colSums(Px_t)-1) > 1e-8 ) ){
  cat("\n\nwarning(obscore) Px_t does not sum to 1.",colSums(Px_t),"\n" )
 }

 if( print > 0 ){
  cat("\n Conditional Distribution of Observed Score given Theta\n")
  Px_t_fr=fuzzy.fmt( Px_t, "%4.2f", 1e-3, revr=1)
  print(Px_t_fr,quote=0)
 }
 ##########################################################################


 # conditional means and variances of the weighted total score X given theta:
 # must be equal to trf and stdx_t^2
 meanx_t=t( t(score)%*%Px_t )
 varx_t=t( t(score^2)%*%Px_t ) - trf^2
 if( debug > 0 ) Print(meanx_t, trf, sqrt(varx_t),stdx_t, digits=2)

 # quantile of Px_t
 # adjust for continuity: add 0.5 to score and use (0,0) and (maxscore_t,1)
 cPx_t=apply( Px_t,2,cumsum )
 dd=(score[2]-score[1])/2
 score1=c(minscore_t,score+dd); score1[length(score1)]=maxscore_t
 cPx_t=rbind(rep(0,npoints),cPx_t);
 rownames(cPx_t)=score1; colnames(cPx_t)=thname
 # Print(score1, cPx_t, digits=2)
 qt=matrix(0,npoints,3)
 colnames(qt)=c(format(alpha/2,digits=2),"mean",format(1-alpha/2,digits=2))
 rownames(qt)=thname
 qt[,2]=trf
 if( debug ) Print(cPx_t, score1)
 # find the quantiles by interpolation
 for( k in 1:npoints ){
  qt[k,c(1,3)]=approx(cPx_t[,k],score1, xout=c(alpha/2, 1-alpha/2)
                      , yleft=minscore_t,yright=maxscore_t)$y
 }
 qtd=qt[,3]-qt[,1]
 if( debug > 0 ) Print(theta,qt,qtd)

 rm(score1, cPx_t)

 if( plot > 0 ){
  title=paste("TRF and ", 100*alpha/2, "% Upper and Lower Quantile")
  matplot( theta, qt, type = "l", ylab="score", col=c(2,1,3), lty=c(2,1,3)
           , xlim=c(min(theta),max(theta)), ylim=c(minscore_t,maxscore_t)
           , main=title  )

  title=paste("Information Function defined as (slope of trf)^2 / var(x|theta)"
              ,", LO and LOW")
  minP=min(cbind(info,info_LO)); maxP=max(cbind(info,info_LO))
  matplot( theta, cbind(info,info_LO,info_LOW)
           , type = "l", ylab="information"
           , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP)
           , main=title  )
  legend( range(theta)[1], max(info_LO)-.1, c("info","info_LO","info_LOW")
          , cex=0.8, col=c("black", "red", "green")
          , pch=NA, lty=c(1,2,3), title="legend" )
  title=paste("S.E. of Theta: inverse of info, info_LO, and info_LOW")
  matplot( theta, cbind(rsqrinfo,rsqrinfo_LO,rsqrinfo_LO)
           , type = "l", ylab="S.E."
           , xlim=c(min(theta),max(theta)), ylim=
           , main=title  )
  legend( theta[3], max(rsqrinfo)-1, c("info","info_LO","info_LOW")
          , cex=0.8, col=c("black", "red", "green")
          , pch=NA, lty=c(1,2,3), title="legend" )

 }

 # confidence interval
 # for each score, solve score = xU(theta) and score = xL(theta) for theta
 ci=matrix(NA,nscore_t,2)
 ci[,1]=approx(qt[,3],theta, xout=score, yleft=NA, yright=NA)$y
 ci[,2]=approx(qt[,1],theta, xout=score, yleft=NA, yright=NA)$y
 cid2=(ci[,2]-ci[,1])/2
 if( debug > 0 ) Print(score,ci,cid2)
 if( debug > 0 ) Print(qtd,slope_trf, digits=2)



 # cond dist of theta given x and mariginal x
 Pt_x=Px_t%*%Diag(Pt)
 Pt_x[is.nan(Pt_x)]=NA
 Px=matrix(rowSums(Pt_x),,1)  # marginal dist of X
 Px[is.nan(Px)]=NA
 rownames(Px)=score
 Pt_x=Diag(1/Px)%*%Pt_x
 Pt_x[is.nan(Pt_x)]=NA
 rownames(Pt_x)=score; colnames(Pt_x)=thname
 if( print > 0 ){
  cat("\n Conditional Distribution of Theta given Observed Score\n")
  Pt_x_fr=fuzzy.fmt( Pt_x, "%4.2f", 1e-3, revr=1, print=1)
 }

 #cat("\n Marginal Distribution of Observed Score\n")
 #Px_fr=fuzzy.fmt( Px, "%6.3f", 1e-3, revr=1, print=1)

 # posterior mean/std of theta given X
 meant_x=Pt_x%*%theta
 stdt_x=sqrt(Pt_x%*%(theta^2)-meant_x^2)
 locnmiss=!is.na(meant_x)

 # stdt_x as the function of theta
 poststd=approx(meant_x,stdt_x, xout=theta, yleft=NA, yright=NA)$y


 # posterior quantile of Pt_x
 # adjust for continuity:
 cPt_x=t( apply( Pt_x,1,cumsum ) )   # after t, cPt_x is score x theta
 dd=(theta[2]-theta[1])/2
 theta1=c(thmin,theta+dd); theta1[length(theta1)]=thmax
 cPt_x=cbind(rep(0,nscore_t),cPt_x);
 colnames(cPt_x)=theta1; rownames(cPt_x)=score
 if( debug ) Print(cPt_x,"/", theta1,score, digits=2)

 qth=matrix(0,nscore_t,2)
 colnames(qth)=c( format(alpha/2,digits=2), format(1-alpha/2,digits=2) )
 rownames(qth)=score

 # find the quantiles by interpolation
 for( i in 1:nscore_t ){
  if( !any( is.na(cPt_x[i,]) ) ){
   qth[i,]=approx(cPt_x[i,],theta1, xout=c(alpha/2, 1-alpha/2)
                  , yleft=thmin,yright=thmax)$y
  }
  else qth[i,]=NA
 }
 qthd2=(qth[,2]-qth[,1])/2
 if( debug > 0 ) Print(score,qth,qthd2)



 if( debug > 0 ) Print(meant_x, stdt_x, ci, cid2, fmt=".3")
 if( print > 9 ){
  cat("\nPosterior Mean etc of Theta given Observed Score\n")
  tab=cbind(Px,meant_x, stdt_x, qth, qthd2, ci, cid2)
  colnames(tab)=c("prob","post.mean","post.std","qth_L","qth_U","qth_wid/2"
                  ,"ci_L","ci_U","ci_wid/2" )
  tabf=matrix(sprintf("%6.3f",tab),nscore_t)
  colnames(tabf)=colnames(tab)
  rownames(tabf)=score
  print(tabf,quote=0)

  c=cor(stdt_x,cid2,use="pairwise.complete.obs")
  cat("\n Correlation coefficient between post.std and the width of ci ="
      , format(c,digits=3),"\n")
 }


 # mean/std of X from marginal dist of X
 xmean=sum(t(score)%*%Px)
 xstd=sqrt(sum(t(score^2)%*%Px)-xmean^2)
 if( debug > 0 ) Print(xmean,xstd)



 # use a subset defined by locnmiss when plotting with type="l"

 if( plot ){
  #  the 1st argument (height) must be a vector, or we need beside=1.
  maint=paste("Plot of the Score Distribution: # of items =", nitems)
  subt=paste(" mean =",format(xmean,digits=3)
             ,",   std =",format(xstd,digits=3))
  barplot(Px, names.arg=score, beside=1, xlab="score"
          , ylab="probability", xlim=
           , main=maint, sub=subt)
 }


 if( plot > 0 ){
  title=paste("Posterior Mean and ", 100*alpha/2, "% Upper and Lower Quantile")
  matplot( score[locnmiss], cbind(meant_x,qth)[locnmiss,]
           , type = "l", ylab="theta"
           , ylim=c(thmin,thmax), xlim=c(minscore_t,maxscore_t)
           , main=title  )
 }


 if( plot > 0 ){
  maint=paste("Posterior Mean and ", 100*(1-alpha)
              , "% CI of Theta given Observed Score")
  subt=paste("alpha =",alpha)
  matplot( score[locnmiss], cbind(meant_x,ci)[locnmiss,]
           , type = "l", ylab="theta"
           , xlim=c(minscore_t,maxscore_t), ylim=c(thmin,thmax)
           , main=maint, sub=subt  )
  maint=paste(
   "Posterior Std and the Half Width of Post Quantile and CI")
  matplot( score[locnmiss], cbind(stdt_x, qthd2, cid2)[locnmiss,]
           , type = "l", ylab="theta"
           , xlim=c(minscore_t,maxscore_t), ylim=c(0,3)
           , main=maint, sub=subt  )
 }


 if( plot > 0 ){
  maint=paste("Posterior Std and SE_LO of Theta")
  subt=paste("r =",format(c,digits=3))
  matplot( theta, cbind(poststd,rsqrinfo_LO), type = "l", ylab="se"
           , xlim=c(thmin,thmax), ylim=
            , main=maint, sub=subt  )
  matplot( theta[!is.na(poststd)]
           , cbind(poststd[!is.na(poststd)],rsqrinfo_LO[!is.na(poststd)])
           , type = "l", ylab="se"
           , xlim=, ylim=
            , main=maint, sub=subt  )
 }

 if( print > 0 ){
  c=cor(poststd,rsqrinfo_LO,use="pairwise.complete.obs")
  cat("\n Correlation coefficient between post.std and the 1/sqrt(info_LO) ="
      , format(c,digits=3),"\n")

  # checking the two ways of calculating var(X|theta)
  adif=abs(sqrt(varx_t)-stdx_t)
  cat("\n Max and mean abs. diff. between two stdx_t:  max =", max(adif)
      , ",   mean =", mean(adif),"\n\n")
 }


 # values to return

 # for each theta point
 theta_stat=as.data.frame( cbind(theta,Pt,trf,slope_trf,stdx_t,info
                                 , info_LOW, info_LO, qt[,1],qt[,3],poststd) )
 colnames(theta_stat)=c("theta","Pt","TRF","slope_TRF","stdx_t","info"
                        , "info_LOW", "info_LO","qt_L","qt_U","poststd")

 # for each score
 obs_stat=as.data.frame( cbind(score,Px,meant_x,stdt_x, qth,qthd2
                               ,ci,cid2) )
 colnames(obs_stat)=
  c("score","Px","meant_x","stdt_x","qth_L","qth_U","qth_wid2"
    ,"ci_L","ci_U","ci_wid2")


 if( print > 0 ){
  Print(as.matrix(theta_stat), obs_stat, fmt=".5")
  c=cor(theta_stat,use="pairwise.complete.obs")
  cat("\n\nCorrelation among Theta_Stat\n")
  fuzzy.fmt( c, "%5.2f", 1e-99, print=1)
  c=cor(obs_stat,use="pairwise.complete.obs")
  cat("\n\nCorrelation among Obs_Stat\n")
  fuzzy.fmt( c, "%5.2f", 1e-99, print=1)
 }

 return( list(theta_stat=theta_stat, obs_stat=obs_stat
              , Px_t=Px_t, Pt_x=Pt_x, Px=Px
              , pdfname=pdfname, wdfname=wdfname
              , npoints=npoints, thmin=thmin, thmax=thmax, thdist=thdist
              , nitems=nitems, maxscore_t=maxscore_t
              , alpha=alpha ) )


} # end of obscore



comments('

paramtest1=data.frame(
 name="Q1", type="B3", ncat=2,
 p1=1, p2=0, p3=0,p4=NA,p5=NA,p6=NA
)
paramtest11=data.frame(
 name="Q11", type="B3", ncat=2,
 p1=1, p2=0, p3=0.3,p4=NA,p5=NA,p6=NA
)

paramtest2=data.frame(
 name="Q2", type="P", ncat=4,
 p1=.7, p2=-1, p3=0, p4=1,p5=NA,p6=NA
)
paramtest21=data.frame(
 name="Q21", type="G", ncat=4,
 p1=.7, p2=-1, p3=0, p4=1,p5=NA,p6=NA
)

paramtest3=data.frame(
 name="Q3", type="P", ncat=5,
 p1=.7, p2=-1.5, p3=-1, p4=1, p5=1.5,p6=NA
)
paramtest31=data.frame(
 name="Q31", type="G", ncat=5,
 p1=.7, p2=-1.5, p3=-1, p4=1, p5=1.5,p6=NA
)

paramtest4=data.frame(
 name="Q4", type="P", ncat=6,
 p1=.7, p2=-1.5, p3=-1, p4=0, p5=1, p6=1.5
)
paramtest41=data.frame(
 name="Q41", type="G", ncat=6,
 p1=.7, p2=-1.5, p3=-1, p4=0, p5=1, p6=1.5
)

param=rbind(paramtest1,paramtest2,paramtest3,paramtest4)
#param=rbind(param, paramtest11,paramtest21,paramtest31,paramtest41 )

param[-1,5:9]=param[-1,5:9]-2


res=obscore( param, npoints=25, thmin=-4, thmax=4, thdist=1
           , compress=0, plot=1, print=2, debug=0 )
#(res)
#Print(res$theta_stat[,grep("info",colnames(res$theta_stat))])
#matplot(res$theta_stat[,"theta"],res$theta_stat[,grep("info",colnames(res$theta_stat))], type="l")

')

res=obscore( paramS2, weight=weightS22, npoints=131, plot=1 )

