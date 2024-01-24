#' Calculation of Observed Score Distribution and Posterior Distribution
#' of Theta with Various Information Functions
#'
#' This function calculates the distributioin of the weighted observed score,
#' the information function associated with it, and the posterior
#' distribution of theta given the observed score. \cr
#' Also, the information functions associated with two types of
#' locally optimam weights will be calculated.
#'
#' @param param Item Parameter Data Frame
#' @param weight Weight data frame
#' @param npoints # of discrete points for theta
#' @param thmin Minimum value of discrete thata value
#' @param thmax Maximum value of discrete thata value
#' @param thdist Type of theta distribution \cr
#' = 0 to use uniform,  = 1 to use N(0,1)
#' @param alpha small prob for quantile and confidence interval
#' @param compress = 1 to remove zero-probability weighted total observed
#'      scores
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#' @param debug = 1 to print intemediate result
#'
#'
#' @return
#' \cr
#'   list( theta_stat, obs_stat, Px_t, Pt_x,  etcetc )\cr
#' \cr
#'    where\cr
#' \cr
#'   theta_stat as data frame containing the following:\cr
#'     theta:       theta points\cr
#'     Pt:          prior probability distribution of theta\cr
#'     TRF:         test response function\cr
#'     slope_TRF:   slope of TRF\cr
#'     stdx_t:      standard deviation of X (observed score) given theta\cr
#'     info:        information function defind as (slope_TRF)^2 / (stdx_t)^2\cr
#'     info_LOW:    information function with the locally optimal item weight
#'     given categoriy weights\cr
#'     info_LO:     information function with the locally optimal category
#'      weight\cr
#'     qt_L:        upper quantile of X given theta\cr
#'     qt_U:        lower quantile of X given theta\cr
#'     poststd:     posterior std of theta given X\cr
#'                    as a function of posterior mean\cr
#' \cr
#'   obs_stat as data frame containing the following:\cr
#'     score:      domain of X (observed score)\cr
#'     Px:         marginal probability  of X\cr
#'     post_mean:  posterior mean of theta given X\cr
#'     post_std:   posterior standard deviation of theta given X\cr
#'     ci_L:       upper limit of conficence interval for theta given X\cr
#'     ci_U:       lower limit of conficence interval for theta given X\cr
#'     ci_hwid:    half the width of CI\cr
#' \cr
#' \cr
#'   Px_t           score x theta  conditional prob of X given theta\cr
#'   Pt_x           score x theta  conditional prob of theta given X\cr
#'   Px             score x 1      marginal dist of X\cr
#' \cr
#' pdfname: item parameter data frame name \cr
#' wdfname: item weight data frame name \cr
#' npoints: # of theta points \cr
#' thmin, thmax: the range of theta \cr
#' thdist: Type of theta distribution \cr
#' nitems: # of items \cr
#' minscore_t: minimum score \cr
#' maxscore_t: maximum score \cr
#' alpha: small probability value \cr
#'
#' @details
#' Note that, given category and item weights, information function is defined
#' as \cr
#' \code{
#'  ( slope of TRF at theta )^2 / (variance of X at theta)
#'  }\cr
#' where TRF and X are calculated with the given set of weights. \cr
#' This is the information function associated with the given item and
#' category weights. \cr
#' \cr
#' The information can be increased by adjusting the weights: \cr
#' However, in general, the optimal weights which maximize the information
#' function depend on the value of theta. \cr
#' Therefore, the name  "locally" optimal weight. \cr\cr
#' The optimal item weight given category weights is called as
#' the locally optimal item weight, or LOW and defined as:. \cr
#' \eqn{ w_j(\theta)= P^{*'}_j(\theta) / var( U_j^{*} | \theta) } \cr
#' where
#' \eqn{ U_j^{*} = \sum_k v_{kj} U_{kj}} and
#' \eqn{ P^{*'}_j(\theta)=\sum_k v_{kj} P_{kj}(\theta) }.
#' \cr
#' When the category weights themselves are optimized it is called as
#' the locally optimal weights, or, LO, \cr
#' which are equivalent to the basic function of Samejima(1969) defined as:.\cr
#' \eqn{ v_{kj}(\theta)=P'_{kj}(\theta) / P_{kj}(\theta) }.
#' \cr\cr
#' The information function with LO is defined as \cr
#' \eqn{ \sum_j \sum_k (P'_{kj}(\theta))^2 / P_{kj}(\theta) } \cr
#' where \eqn{P_{kj}(\theta)} is the item category response function, \cr
#' and \eqn{P'_{kj}(\theta)} is its derivative. \cr
#' This is the maximum information function for all the theta range.
#' \cr
#' The information function with LOW is defined as \cr
#' \eqn{ \sum_j  (P_j^{*'}(\theta))^2 / var(U_j^{*} | \theta) } \cr
#' where \eqn{U_j^{*} = \sum_k v_{kj} U_{kj}} is
#'  the weighted item score,\cr
#' and \eqn{P_j^{*'}(\theta)} is the derivative of the expected value of
#' \eqn{U_j^{*}} at theta. \cr
#'
#'
#'
#' @references
#' Birnbaum, A.(1968) Some Latent Traint Models.
#' In F. M. Lord and M. R. Novick, Statistical Theories of Mental Test Scores.
#'  Reading, Mass.: Addison-Wesley. \cr
#'
#' Mayekawa, S., & Arai, S. (2008).
#' Distribution of the Sum of Scored Multinomial Random Variables
#' and Its Application to the Item Response Theory.
#' In K. Shigemasu, A. Okada, T.Imaizumi, & T. Hoshino (Eds.)
#' New Trends in Psychometrics. Tokyo: University Academic Press. \cr
#'
#' Samejima,  F.  (1969). Estimation  of  a  latent  ability  using  a
#' response  pattern  of  graded  scores. Psychometrika  Monographs, 34
#' (Suppl. 4).
#'
#' @examples
#' # Define the observed raw score X
#' # and calculate the score distribution, information function, etc.
#' res <- obscore( paramS1, plot=1 )
#'
#' # Define X using the item weights w stored in weightS11
#' res <- obscore( paramS1, weight=weightS11, plot=1 )
#'
#' # Define X using the item and category weights w and v stored in weightS12.
#' res <- obscore( paramS1, weight=weightS12, plot=1 )
#'
#' @export
#'

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
 # print Pt: 20151114
 # comments etc: 20170120dnc
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
 thdisttype=c("uniform","normal(0,1)")
 Pt_fr=formatC( Pt, format="f", digits=2, width=4 )
 names(Pt_fr)=thname

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
  cat( " theta distribution type = ", thdisttype[thdist+1],"\n")
  cat(" alpha for quantile and CI =", alpha,"\n")
  cat(" range of observed score = [",minscore_t,",", maxscore_t,"]\n")
  cat("  (total of", nscore_t," discrete score points with ",nmiss
      ," zero probability points)\n")
  cat(" removal of 0-probability scores from the result = ", compress, "\n")
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
  cat("\n Marginal Theta Distribution\n")
  print(Pt_fr,quote=0)
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
              , nitems=nitems, maxscore_t=maxscore_t, minscore_t=minscore_t
              , alpha=alpha ) )


} # end of obscore









#' Distribution of the Weighted Sum of Several Independent
#' Scored Multinomial Distributions
#'
#'
#' @param P matrix of probabilities   (max # of categories x # of r.v.)
#' @param V matrix of domain values   (max # of categories x # of r.v.) or NULL
#' @param w vector of weights         (1 x # of r.v.)
#' @param ncat max # of categories       (1 x # of r.v)    or NULL
#' @param compress = 1 to remove the zero probability categories
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#' @param debug = 1 to print intemediate result
#'
#' @details
#'  ( V[,i], P[,i], w[i] ), i=1,2,...,n is the set of\cr
#'  ( domain or category weight ,probability, and weight ) for the i-th r.v.\cr
#' \cr
#' Non integer V and w will be first converted to integer by linear
#' transformation and converted back at the very end. \cr
#' \cr
#'  ncat[i] = max # of categories for the i-th r.v and\cr
#'  P[(ncat[i]+1):nrow(P),i] == NA\cr
#'\cr
#'   This program calculats the distrobution of\cr
#'     \eqn{X = sum_{i=1}^n w[i] X_i} \cr
#'    where\cr
#'     \eqn{X_i} is distributed as Scored Multinomial with (V[,i],P[,i])
#'  ,i=1,2,...,n\cr
#'     \eqn{V[1,i] <=  X_i  <= V[ncat[i],i]  or  0  <=  X_i  <=  ncat[i]}\cr
#'   That is, this program calculates the probability that\cr
#'     \eqn{Pr( X = a )} ,\cr
#'    where\cr
#'     \eqn{sum_{i=1}^n V[1,i]*w[i] <= a <=  sum_{i=1}^n V[ncat[i],i]*w[i]}\cr
#'
#' @return A matrix of (score, prob)
#'
#'
#' @examples
#' # category x variable matrix of probability:  colSums(P)=c(1,1,1...)
#' P <- matrix(c(1,2,3,2,   1,2,1,0,  1,2,0,0), 4,3)
#' P <- t(t(P)/colSums(P))
#' ncat <- c(4,3,2)
#' # category x variable matrix of category weight
#' V <- NULL
#' # variable weight
#' w <- c(.5,1,1)
#' res <- sumsmnw( P, V, w, compress=0, print=1, plot=1, ncat=ncat )
#'
#' # category x variable matrix of category weight
#' V <- matrix(c(0,1,2,3,   1,2,3,0,  1,1,0,0), 4,3)
#' # variable weight
#' w <- c(.5,1,1)
#' res <- sumsmnw( P, V, w, compress=0, print=1, plot=1, ncat=ncat )
#'
#' @references Mayekawa, S., & Arai, S. (2008).
#' Distribution of the Sum of Scored Multinomial Random Variables
#' and Its Application to the Item Response Theory.
#' In K. Shigemasu, A. Okada, T.Imaizumi, & T. Hoshino (Eds.)
#' New Trends in Psychometrics. Tokyo: University Academic Press.
#'
#' @export
#'
sumsmnw <- function( P, V=NULL, w=rep(1,ncol(P)), ncat=NULL, compress=0
                   , print=0, plot=0, debug=0 ){
 # sum of independent scored multinomial random variables
 # Shin-ichi Mayekawa
 # iml version: 000823 -- 080106
 # 120204,05
 # use V%*%diag(w) : 120205
 # bugfix: 120206
 # plot: 120207
 # print ncat: 120210
 # comments added: 130804
 # bug fix: 20150305
 # ncat: 20150305,07
 #
 #
 # Args:
 #
 #  P     matrix of probabilities   (max # of categories x # of r.v.)
 #  V     matrix of domain values   (max # of categories x # of r.v.) or NULL
 #  w     vector of weights         (1 x # of r.v.)
 #  ncat  max # of categories       (1 x # of r.v)    or NULL
 #
 #  ( V[,i], P[,i], w[i] ), i=1,2,...,n is the set of
 #  ( domain ,probability, and weight ) for the i-th r.v.
 #
 #  ncat[i] = max # of categories for the i-th r.v and
 #  P[(ncat[i]+1):nrow(P),i] == NA
 #
 #   This program calculats the distrobution of
 #     X = sum_{i=1}^n w[i] X_i
 #    where
 #     X_i is distributed as Scored Multinomial with (V[,i],P[,i]),i=1,2,...,n
 #     V[1,i] <=  X_i  <= V[ncat[i],i]   or  0  <=  X_i  <=  ncat[i]
 #   That is, this program calculates the probability that
 #     Pr( X = a ) ,
 #    where
 #     sum_{i=1}^n V[1,i]*w[i] <= a <=  sum_{i=1}^n V[ncat[i],i]*w[i]
 #
 #

 # const
 n=ncol(P)   # # of random variables to be added
 maxncat=nrow(P)
 if( is.null( ncat ) ) ncat=colSums(!is.na(P)) # # of categories for each r.v.
 else
  for( i in 1:n )
   if( ncat[i] < maxncat ) P[(ncat[i]+1):maxncat,i]=NA

 # normalize P so that colsums=1
 locna=which( is.na(P) )
 if( length(locna) > 0 ){
  P[locna]=0
 }
 P=t(t(P)/colSums(P))
 if( length(locna) > 0 ){
  P[locna]=NA
 }
 P0=P
 if( debug ) Print(ncat,P)

 # weight for each r.v.
 if( is.matrix(w) ) w=as.vector(w)

 # category score: v-score: 0,1,...,ncat
 if( is.null(V) ) V=matrix(1:nrow(P)-1,nrow(P),ncol(P))
 V0=V

 if( debug ) Print(P,V,n,ncat)

 # combine the cells with duplicate v-scores
 # There ara nuch room for improvement.
 P1=matrix(NA,nrow(P),ncol(P))
 V1=matrix(0,nrow(P),ncol(P)) # cannot have NA in V
 chksum=numeric(n)
 for( k in 1:n ){
  # sort (p[,k],v[,k]) according to v[,k] and find tie info
  pk=P[,k][1:ncat[k]]; vk=V[,k][1:ncat[k]]
  od=order(vk); pk=pk[od]; vk=vk[od]
  tie=as.matrix(table(vk)); nc=nrow(tie)
  to=cumsum(tie); from=to-tie+1
  if( debug ) Print(k,pk,vk,tie,from,to,nc)
  pk1=numeric(nc)
  for( kk in 1:nc )
   pk1[kk]=sum(pk[from[kk]:to[kk]])
  P1[1:nc,k]=pk1; V1[1:nc,k]=unique(vk)
  chksum[k]=sum(pk1)
 }
 ncat1=colSums(!is.na(P1))  # new # of categories for each r.v.
 if( debug ) Print(ncat1,P1,V1,chksum)

 P=P1; V=V1;

 # make V nonnegative integer by subtracting min and mulplying dd
 # maxv=apply(V,2,max,na.rm=1)
 # minv=apply(V,2,min,na.rm=1)
 minv=V[1,]  # min is the first value, not 0 nor NA
 maxv=V[matrix(cbind(ncat1,1:n),,2)]  # max is the last NA value

 # combine V and w into V2
 V2=V-matrix(minv,nrow(V),n,byrow=1)
 V2=V2%*%diag(w)   # if V has NA, this fails.
 if( debug ) Print(minv,maxv,V,V2)

 found=0; dd=0;  maxdec=0;
 for( d in c(1,2,3,4,5,10,15,20,25,40,50) ){
  #  if( all(d*V2-floor(d*V2) == 0) ){
  if( isTRUE(all.equal(d*V2,floor(d*V2))) ){
   V2=d*V2; found=1; dd=d
   break
  }
 }
 if( !found ){
  maxdec=-1
  for( k in 1:n ){
   if( debug ) Print(k,decp(V[,k]))
   maxdec=max(maxdec,decp(V1[,k])[,3])
  }
  if( any( maxdec > 2 ) ){
   cat("\n\nerror:(sumsmnw) category values too small.  V=", V,"\n")
   return( NULL )
  }
  dd=10^maxdec
  V2=dd*V2
 }

 V2[is.na(P)]=NA
 minv2=V2[1,]
 maxv2=V2[matrix(cbind(ncat1,1:n),,2)]
 ncat2=maxv2+1
 if( debug ) Print(ncat1,V2,maxv2)

 # Now, minv2[k]=0 <= V2[,k] <= ncat1[k]*w[k]=maxv2[k]
 # and there are ncat2[k] distinct values of V2[,k].

 # score of the sum
 score=( 0:sum(maxv2) )/dd + sum(minv*w)

 # expand P1 according to V2
 P2=matrix(0,max(maxv2)+1,n)
 rownames(P2)=0:max(maxv2)
 if( debug ) Print(P1,V2,ncat1)
 for( k in 1:n ){
  P2[V2[1:ncat1[k],k]+1,k]=P1[1:ncat1[k],k]
 }

 if( debug ){
  Print(n,ncat,ncat1, maxdec,dd,minv,maxv,minv2,maxv2,w,ncat2)
  Print(P0,P1,P2, "/", V0,V1,V2)
 }


 # main body
 p3=P2[1:ncat2[1],1]
 for( k in 2:n ){
  if( debug ) Print(k,p3,P2[1:ncat2[k],k])
  p3=sumsmnw12( p3,P2[1:ncat2[k],k],1 )[,2] # always use unweighted sum.
 }
 chksum=sum(p3)

 if( compress == 1 ){
  loc=which(p3 > 0)
  p3=p3[loc]; score=score[loc]
 }

 if( print > 0 ){
  offset=minv*w
  cat("\n\nSum of Scored Multinomial Distributions")
  cat("\n # of Scored Multinomials to be added =", n, "\n")
  cat("\nInput category weight and probability matrices and variable weight\n")
  Print(ncat, V0,P0,w)
  if( !isTRUE(all.equal( P0[!is.na(P0)],P1[!is.na(P0)] )) ){
   cat("\nAfter collecting duplicate v-values\n")
   Print(V1,P1,w)
  }
  cat("\nV <- V%*%diag(w) and After expantion of P. \n")
  Print(V2,P2,"/",offset,dd,maxdec,compress)
  cat("\nThe result\n")
  Print(score, p3, chksum, fmt="6.2, 5.3")
 }

 if( plot ){
  barplot(p3, names.arg=score, xlab="score", ylab="probability"
          , main="Plot of the Score Distribution")
 }

 res=cbind(score,p3)
 colnames(res)=c("score","p")
 rownames(res)=NULL
 return( res )


} # end of sumsmnw



#' Distribution of the Weighted Sum of Two Independent
#' Scored Multinomial Distributions
#'
#'
#' @param p1 Probability of the fist random variable
#' @param p2 Probability of the second random variable
#' @param w2 Integer weight to the second random variable
#' @param eps small number
#' @param print = 1 to print result
#'
#' @details
#'    This program calculats the distrobution of\cr
#'      X12 = X1 + w2 X2\cr
#'     where\cr
#'      Xi  is distributed as Scored Multinomial with  pi,  i=1,2.\cr
#'    That is, this program calculates the probability that\cr
#'      p12[a] = Pr( X12 = a ) , a=0,1,..., n12=(n1-1)+w2*(n2-1)\cr
#' \cr
#'  Note\cr
#'    length(pi) = ni = #'  of categories of the i-th r.v.\cr
#'               = mi + 1   where mi is the max score of the i-th r.v.\cr
#' \cr
#'
#' @return A matrix of (score, prob)
#'
#'
#' @examples
#' sumsmnw12( 1:3, 1:2, 2, print=1 )
#'
#' @references Mayekawa, S., & Arai, S. (2008).
#' Distribution of the Sum of Scored Multinomial Random Variables
#' and Its Application to the Item Response Theory.
#' In K. Shigemasu, A. Okada, T.Imaizumi, & T. Hoshino (Eds.)
#' New Trends in Psychometrics. Tokyo: University Academic Press.
#'
#' @export
#'

sumsmnw12 <- function( p1, p2, w2=1, eps=1e-8, print=0 ){
 # sum of two independent scored multinomial random variables
 # Shin-ichi Mayekawa
 # iml version: 000823 -- 080106
 # 120204
 # Z matrix: 120207
 #
 # Args:
 #   p1 and p2   scored multinomial probabilities
 #    the domains are 0:n1=length(p1) and 0:n2=length(p2), resp
 #
 #   w2          integer weight for the 2nd r.v.
 #
 #   This program calculats the distrobution of
 #     X12 = X1 + w2 X2
 #    where
 #     Xi  is distributed as Scored Multinomial with  pi,  i=1,2.
 #   That is, this program calculates the probability that
 #     p12[a] = Pr( X12 = a ) , a=0,1,..., n12=(n1-1)+w2*(n2-1)
 #
 # Note
 #   length(pi) = ni = # of categories of the i-th r.v.
 #              = mi + 1   where mi is the max score of the i-th r.v.
 #

 # const
 if( is.matrix(p1) ) p1=as.vector(p1)
 if( is.matrix(p2) ) p2=as.vector(p2)
 n1=length(p1); n2=length(p2)
 n12=(n1-1)+w2*(n2-1)+1
 score=1:n12-1

 p1=p1/sum(p1)
 p2=p2/sum(p2)

 if( w2 != floor(w2) ){
  cat("\n\nerror:(sumsmnw12) w2 must be an integer.  w2=", w2,"\n")
  return( NULL )
 }

 # main body
 Z=matrix(0,n12,n2)
 rownames(Z)=score
 for( j in 1:n2 ){
  Z[,j]=c( rep(0,w2*(j-1)), p2[j]*p1, rep(0,w2*(n2-j)) )
 }
 p12=rowSums(Z)
 chksum=sum(p12)

 if( abs(chksum-1) > eps ){
  cat("\n\nwarning(sumsmnw12): Big Trouble!! \n")
  cat(" Sum of probabilities calculated is not equa to 1.", chksum, "\n\n")
 }
 res=cbind(score,p12)
 colnames(res)=c("score","p")
 rownames(res)=rep("",n12)

 if( print > 0 ){
  cat("\n\nSum of Two Scored Multinomial Random Variables\n")
  Print(n1,p1,n2,p2,w2, digits=3)
  if( print >= 2 ){
    cat("\n Z matrix whose row sum is the desired probability.\n")
    print(Z,digits=3)
   }
  Print(score,p12,chksum)
  }

 return( res )

} # end of sumsmnw12



