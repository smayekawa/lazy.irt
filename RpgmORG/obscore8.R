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
#' @param weight Weight data frame \cr
#' or a vector containing item weights. \cr
#' or a scalar or NULL.
#' @param npoints # of discrete points for theta
#' @param thmin Minimum value of discrete thata value
#' @param thmax Maximum value of discrete thata value
#' @param thdist Type of theta distribution \cr
#' = 0 to use uniform,  \cr = 1 to use N(0,1)
#' @param method = 0 to use exact method, = 1 to use normal probability,
#' = 2 to use normal pdf. \cr
#'  See the comments on graded_prob.
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
#'     info:        information function defind as
#'     (slope_TRF)^2 / (stdx_t)^2\cr
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
#'   SigmaX2, SigmaT2, SigmaE2: (observed, true, and average error variances)
#' \cr
#'   aSEM=sqrt(SigmaE2): average Standard Error of Measurement, \cr
#'   rel: test reliability
#' \cr
#'   SigmaT_theta2, SigmaE_theta2, aSEM_theta, rel_theta
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
#' This is the information function associated with the test scores \cr
#' calculated using the given item and  category weights. \cr
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
#' which are related to the basic function of Samejima(1969)
#' defined as:.\cr
#' \eqn{ v_{kj}(\theta)=P'_{kj}(\theta) / P_{kj}(\theta)
#'  - P'_{0j}(\theta) / P_{0j}(\theta) }.
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
#' \cr
#' The test reliability coefficient is calculated as
#' \code{rel = (SigmaX^2-SigmaE^2)/SigmaX^2}, \cr
#' where SigmaE^2 is the square root of the average of stdx_t^2, \cr
#' SigmaX^2 is the variance of the observed score.
#' \cr
#'
#' @references
#' Birnbaum, A.(1968) Some Latent Traint Models.
#' In F. M. Lord and M. R. Novick, Statistical Theories of Mental Test Scores.
#'  Reading, Mass.: Addison-Wesley. \cr
#'
#' Mayekawa, S. (2018) A Method to Estimate the Ability from the Weighted
#' Total Score wich Minimizes the Standard Error of Estimation.
#' DNC Research Note. RN-18-01. \cr
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
                   , npoints=31, thmin=-4, thmax=4, thdist=1
                   , method=0, alpha=0.1
                   , compress=0, print=1, plot=0, debug=0 ){
 # calculation of observed score distribution
 # Shin-ichi Mayekawa
 # The original version was developped using sas 5.18 in 19890124-0404.
 # sas/iml version: 19910115-20080101
 # R version:
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
 # history added: 20170418dnc
 # reliability: 20170419
 # rel_theta: 20171204,05
 # simple weight input: 20180629
 # method: 20180708,13 ****************** NOT READY? 20230717
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


 fuzzy.fmt <- function( A, fmt, fuzz=-1, fzchar="."
                        , revr=0, revc=0, print=0 ){
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

 # natural weight
 w=matrix(1,nitems); rownames(w)=iname; colnames(w)="w"
 v=matrix(0:(max(ncat)-1),nitems,max(ncat),byrow=1)
 rownames(v)=iname; colnames(v)=paste("v",0:(max(ncat)-1),sep="")
 for( i in 1:nitems ){
  if( ncat[i]+1 <= ncol(v) ) v[i,(ncat[i]+1):ncol(v)]=NA
 }

 # weight data frame
 if( is.null(weight) ){
  weight=data.frame(iname, type, ncat, w, v)
  rownames(weight)=iname;
  colnames(weight)=c("name","type","ncat","w",colnames(v))
 }
 else if( is.vector(weight) ){
  if( length(weight) %in% c(1,nitems) )
   weight=data.frame(iname, type, ncat, weight, v)
  else
   weight=data.frame(iname, type, ncat, w, v)
  rownames(weight)=iname
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
  cat( "   mean theta =", thmean, ",   var theta =", thvar, "\n")
  cat( " method to calculate conditional probs = ", method, "\n")
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
 slopetrf=trf
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


 # Print( trf, stdx_t, score )


 if( method %in% 1:2 ){
  Px_t=matrix(0,maxscore_t+1,npoints)
  rownames(Px_t)=score; colnames(Px_t)=thname
  for( q in 1:npoints ){
   Px_t[,q]=graded_prob( score, trf[q], stdx_t[q], truncate=1, method=method  )
  }
 # Px_t_fr=fuzzy.fmt( Px_t, "%4.2f", 1e-3, revr=1)
 # print(Px_t_fr,quote=0)

 } # end of method 1 or 2
 else{

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

  ##########################################################################

 } # end of core part

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

 # conditional means and variances of the weighted total score X given theta:
 # must be equal to trf and stdx_t^2
 meanx_t=t( t(score)%*%Px_t )
 varx_t=t( t(score^2)%*%Px_t ) - meanx_t^2
Print(trf, meanx_t)
Print( score, Px_t, fmt="6.3")
Print(varx_t)
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

  title=
  paste("Information Function defined as (slope of trf)^2 / var(x|theta)"
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
 xmean=sum(score*Px)
 xstd=sqrt( sum(score^2*Px)-xmean^2 )
 if( debug > 0 ) Print(xmean,xstd)

 # reliability etc
 # average standard error of measurement
 SigmaE2=sum(Pt*varx_t)
 SigmaX2=xstd^2
 SigmaT2=SigmaX2-SigmaE2
 aSEM=sqrt(SigmaE2)
 rel=SigmaT2/SigmaX2
 locm=which.min(abs(theta-thmean))
 # SigmaT22=slope_trf[locm]^2*thvar
 SigmaT22=sum(trf^2*Pt)-sum(trf*Pt)^2

 # average standard error of estimation of theta
 # info_LO
 SigmaT_theta2=sum(Pt*theta^2)-sum(Pt*theta)^2
 SigmaE_theta2=sum(Pt/info_LO)
 aSEM_theta=sqrt(SigmaE_theta2)
 rel_theta=1-SigmaE_theta2/SigmaT_theta2

 if( print ){
  Print( SigmaX2, SigmaT2, SigmaT22, SigmaE2, aSEM, rel )
  Print( SigmaT_theta2, SigmaE_theta2, aSEM_theta, rel_theta  )
 }


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
  title=
   paste("Posterior Mean and ", 100*alpha/2, "% Upper and Lower Quantile")
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
  Print(theta_stat, obs_stat, fmt=".5")
  c=cor(theta_stat,use="pairwise.complete.obs")
  cat("\n\nCorrelation among Theta_Stat\n")
  fuzzy.fmt( c, "%5.2f", 1e-99, print=1)
  c=cor(obs_stat,use="pairwise.complete.obs")
  cat("\n\nCorrelation among Obs_Stat\n")
  fuzzy.fmt( c, "%5.2f", 1e-99, print=1)
 }

 return( list(theta_stat=theta_stat, obs_stat=obs_stat
              , SigmaX2=SigmaX2, SigmaT2=SigmaT2, SigmaT22=SigmaT22
              , SigmaE2=SigmaE2, aSEM=aSEM, rel=rel
              , SigmaT_theta2, SigmaE_theta2, aSEM_theta, rel_theta
              , Px_t=Px_t, Pt_x=Pt_x, Px=Px
              , pdfname=pdfname, wdfname=wdfname
              , npoints=npoints, thmin=thmin, thmax=thmax, thdist=thdist
              , nitems=nitems, maxscore_t=maxscore_t, minscore_t=minscore_t
              , method=method, alpha=alpha ) )


} # end of obscore



res <- obscore( paramB1, npoints=5, plot=1, method=0 )



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

param=paramS1

res=obscore( param, npoints=25, thmin=-4, thmax=4, thdist=1
           , compress=0, plot=1, print=2, debug=0 )
#(res)
#Print(res$theta_stat[,grep("info",colnames(res$theta_stat))])
#matplot(res$theta_stat[,"theta"],res$theta_stat[,grep("info",colnames(res$theta_stat))], type="l")


res=obscore( paramS2, weight=weightS22, npoints=131, plot=1 )


')
