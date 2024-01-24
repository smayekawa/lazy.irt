#' Calculation of Various Information Functions
#'
#' This function calculates the information functions associated with
#' the set of given weights, \cr
#' and two types of locally optimam weights.
#'
#' @param param Item Parameter Data Frame.
#' @param theta Discrete theta values
#' @param weight Weight data frame.
#' @param npoints # of discrete points for theta.
#' @param thmin Minimum value of discrete thata value.
#' @param thmax Maximum value of discrete thata value.
#' @param numderiv = 1 to use numerical derivatives.
#' @param smallP Minimum value of probability in irf and dirf functions.
#' @param print = 1 to print the summary \cr
#' = 2 to print item information functions. \cr
#' = 3 to print the locally best item and item categorie weights.
#' @param plot = 1 to plot test information functions\cr
#' = 2 to plot item information functions. \cr
#' = 3 to plot the locally best item and item categorie weights.
#' @param legend = 0 to skip printing legend.
#' @param maxinfo = The maximum value of item information function for plot.
#' @param debug = 1 to print intemediate result.
#'
#'
#' @return
#' A list of \cr
#'     theta    theta points\cr
#'     info    information function defind as (slope_TRF)^2 / (stdx_t)^2\cr
#'     fromP, toP location of each item category in item info \cr
#'     TRF test response function (tcc) \cr
#'     TRF_LO test response function (tcc) with locally best weight\cr
#'     info_LOW   information function with locally optimal item weight\cr
#'     info_LO   information function with locally optimal category weight\cr
#'     info_item_LOW  item information function
#'                 with locally optimal item weight\cr
#'     info_item_LO   item information function
#'                 with locally optimal category weight\cr
#'     w_LO the locally best item weight \cr
#'     v_LO the locally best item category weight
#'     scaled so that v[0,j]=0, j=1,2,...,nitems.
#'
#' @details
#' Note that, given category and item weights, information function is defined
#' as \cr
#'  ( slope of TRF at theta )^2 / (variance of x at theta) \cr
#' where TRF and x is calculated with the given set of weights. \cr
#' \cr
#' In general. the optimal weights depends on the value of theta.
#' Therefore, the name  locally optimal weight. \cr
#' The optimal item weight given category weights is called as
#' locally optimal item weight, or LOW. \cr
#' When the category weights themselves are optimized it is called as
#' the locally optimal weight, or, LO. \cr
#' This weight is equivalent to the basic function of Samejima(1969) .\cr
#' \cr
#' The information function with LO is defined as \cr
#' \eqn{ \sum_j \sum_k (P'_{kj}(\theta))^2 / P_{kj}(\theta) } \cr
#' where \eqn{P_{kj}(\theta)} is the item category response function, \cr
#' and \eqn{P'_{kj}(\theta)} is its derivative. \cr
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
#' Samejima,  F.  (1969). Estimation  of  a  latent  ability  using  a
#' response  pattern  of  graded  scores. Psychometrika  Monographs, 34
#' (Suppl. 4).
#'
#' @examples
#' resInfo <- info_func( paramS2, plot=2, print=1 )
#' resInfo <- info_func( paramS2, weight=weightS21, plot=1, print=0 )
#' resInfo <- info_func( paramS2, weight=weightS22, plot=1, print=0 )
#'
#' @export
#'

info_func <- function( param, theta=NULL, weight=NULL
                       , npoints=31, thmin=-4, thmax=4
                       , legend=1, maxinfo=0, numderiv=0, smallP=1e-9
                       , print=1, plot=0, debug=0 ){
 # calculation of information function
 # Shin-ichi Mayekawa
 # obscore modified: 20140329
 # legend: 20150616
 # bugfix for LOW: 20150616,17
 # plot item info: 20170627
 # v_LO and wLO: 20171209
 # trf with v_LO: 20171209
 # return trf etc: 20171219dnc
 # do not use v[1,1]=1: 20171229
 # output icrf and dicrf, etc: 20180109
 # max(trf_LO) not adjusted: 20180109
 # numderiv: 20180111
 # smallP: 20180129
 # theta as an arg: 20180201
 # basic functioin: 20180209
 #
 # Args:
 #    param    parameter data frame name
 #    weight   item and category weight data frame name
 #             see the descriptions in read.param or read.weight.
 #             If weight==NULL a set of narural weights will be used.
 #
 #    npoints  # of theta points
 #    thmin    min value of theta
 #    thmax    max value of theta
 #    print    = 1 to verbose
 #    plot     = 1 to produce several plot
 #
 #
 #
 # Values:
 #
 #  list( info, info_LO, info_LOW )
 #
 #
 #
 # Needs:
 #  irf,  icrfB, icrfG,  icrfPN
 #  dirf,  dicrfB, dicrfG,  dicrfPN
 #  sumsmnw,  sumsmnw12
 #  checkparam
 #

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
 locv=which(colnames(weight) == "v0")
 v=weight[,locv:ncol(weight)]

 maxscore_i=apply(weight[,5:ncol(weight)],1,max,na.rm=1)
 maxscore_t=sum(w*maxscore_i)
 minscore_i=apply(weight[,5:ncol(weight)],1,min,na.rm=1)
 minscore_t=sum(w*minscore_i)

 if( debug > 0 ) Print(iname,nitems,nc,ncat,w, "/",maxscore_i, maxscore_t)

 # generate thata and prior theta dist
 if( is.null(theta) ){
  theta=seq(thmin,thmax,length.out=npoints)
 }
 else{
  if( is.matrix(theta) ) theta=as.vector(theta)
  npoints=length(theta)
 }
 thname=format(theta,digits=2)


 if( print > 0 ){
  cat("\n\nCalculation of")
  cat(" the Information Function \n")
  cat(" parameter data frame name =", pdfname
      , ",  item category weight data frame name =", wdfname,"\n\n")
  cat(" # of item parameters =",nitems,"\n")
  cat(" range of theta = [",thmin,",",thmax,"] with", npoints
      ,"discrete points\n")
  cat("\n item parameters\n")
  print( param )
  cat("\n item and item category weight \n")
  print( weight )
  cat("\n")
 }


 # calculate icrf and trf=cond. mean of X given theta
 temp=irf( param, theta, weight, print=0, debug=0, plot=1, smallP=smallP )
 icrf=temp$ICRF
 TRF=temp$TRF
 fromP=temp$fromP
 toP=temp$toP
 vecv=temp$vecv
 rm(temp)

 if( debug > 0 ) Print(icrf, fromP,toP)

 # first derivative of icrf
 temp=dirf( param, theta, weight=weight, print=0, numderiv=numderiv
            , smallP=smallP )
 dicrf=temp$dICRF
 slope_trf=temp$dTRF
 dirf=temp$dIRF
 rm(temp)

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

 # information function with the locally optimal category weights
 dd=(dicrf^2)/icrf
 info_LO=rowSums(dd)
 info_item_LO=mapply( function(x,y) rowSums(dd[,x:y])
                      , fromP, toP, SIMPLIFY=FALSE )
 info_item_LO=matrix(unlist(info_item_LO),,nitems)
 colnames(info_item_LO)=iname
 rownames(info_item_LO)=thname
 names(info_LO)=thname

 # information function with the locally optimal item weights
 info_item_LOW=(dirf^2)/varuv_t
 info_LOW=rowSums(info_item_LOW)
 names(info_LOW)=thname

 # locally best item weight
 w_LO=dirf/varuv_t

 # locally best category weight: scaled so that v[0,j]=0.
 A=dicrf/icrf
 v_LO=mapply( function(x,y){A[,x:y]-A[,x]}, fromP, toP, SIMPLIFY=FALSE )
 v_LO=matrix(unlist(v_LO),npoints)
# v1=v_LO[,2]
# v_LO=mapply( function(x,y){v_LO[,x:y]/v1}, fromP, toP, SIMPLIFY=FALSE  )
# v_LO=matrix(unlist(v_LO),npoints)


 colnames(v_LO)=colnames(dicrf); rownames(v_LO)=thname
 colnames(A)=colnames(dicrf); rownames(A)=thname

 # trf with the locally best item category weight
 TRF_LO=rowSums( icrf*v_LO )

 # TRF_LO=TRF_LO/max(TRF_LO)*max(TRF)

 if( print > 0 ){
  cat( " info is the information function associated with given category and"
       , "item weights,\n"
       , "info_LO uses optimal category and item weights,\n"
       , "info_LOW uses optinal item weights given category weights.\n" )
  Print(theta,info,info_LO, info_LOW, fmt="6.3")
  if( print >= 2 ){
   cat(" Item information functions for various weights.\n")
   Print(info_item_LO, fmt="6.3")
   Print(info_item_LOW, fmt="6.3")
  }
  if( print >= 3 ){
   cat(" Locally Best Item Weights.\n")
   Print(w_LO, dirf, varuv_t, fmt="6.3")
   cat(" Locally Best Item Categorie Weights: scaled so that v[0,j]=0. \n")
   Print(v_LO, fmt="6.3")
  }
 }


 if( plot > 0 ){
  # trf

  matplot( theta, cbind(TRF, TRF_LO), type="l", col=1, lwd=2, lty=c(1,3)
           , main="TRF and TRF with the locally best item weights" )
  legend( "topleft", c("TRF","TRF with LO")
          , cex=0.8, col=1, lwd=3, lty=c(1,3)
          , pch=NA, title="legend" )


  # plot test info
  matplot( theta,cbind(info, info_LO, info_LOW), type="l"
           , xlim=c(min(theta),max(theta)), ylim=c(0,max(info_LO))
           , main="Information Functions", ylab="information")
  legend( range(theta)[1], max(info_LO)-.1, c("info","info_LO","info_LOW")
          , cex=0.8, col=c("black", "red", "green", "blue")
          , pch=NA, lty=c(1:3), title="legend" )

  # plot item info
  if( plot >= 2 ){
   matplot( theta,cbind(info_item_LO), type="l"
            , col=c("black", "red", "green", "blue"), lty=c(1:3)
            , xlim=c(min(theta),max(theta)), ylim=c(0,max(info_item_LO))
            , main="Item Information Functions", ylab="information")
   if( legend > 0 ){
    legend( range(theta)[1], max(info_item_LO)-.1, iname
            , cex=0.8, col=c("black", "red", "green", "blue")
            , pch=NA, lty=c(1:3), title="legend" )
   }

   for( j in 1:(nitems) ){
    title=paste("Item Information Function of "
                , iname[j]," ( type = ",type[j],", ncat=",ncat[j], " )")
    pp=param[j,4:(ncat[j]+3)]
    pp[is.na(pp)]=0
    sub=paste("param = ", paste(format(pp,digits=3)
                                , collapse=",  "), "  (with weights)")
    if( length(grep("^B[[:digit:]]*$", param$type[j])) > 0 ){
     pp=param[j,4:(ncat[j]+4)]
     pp[is.na(pp)]=0
     sub=paste("param = ", paste(format(pp,digits=3), collapse=",  "))
    }
    maxi=max(info_item_LO)
    if( maxinfo > 0 ) maxi=maxinfo
    plot(theta,info_item_LO[,j], main=title, sub=sub
         , xlim=c(min(theta),max(theta)), ylim=c(0,maxi)
         , type="l", ylab="information")
   }
  }
  if( plot >= 3 ){
   matplot( theta,w_LO, type="l"
            , col=c("black", "red", "green", "blue"), lty=c(1:3)
            , xlim=c(min(theta),max(theta)), ylim=c(0,max(w_LO))
            , main="Locally Best Item Weights", ylab="weight")
   if( legend > 0 ){
    legend( range(theta)[2]-1, range(w_LO)[2]/4*3, iname
            , cex=0.8, col=c("black", "red", "green", "blue")
            , pch=NA, lty=c(1:3), title="legend" )
   }
   for( j in 1:(nitems) ){
    title=paste("Locally Best Item Category Weights for"
                , iname[j]," ( type = ",type[j],", ncat=",ncat[j], " )")
    pp=param[j,4:(ncat[j]+3)]
    pp[is.na(pp)]=0
    sub=paste( "param = ", paste(format(pp,digits=3)
                                , collapse=",  ") )
    if( length(grep("^B[[:digit:]]*$", param$type[j])) > 0 ){
     pp=param[j,4:(ncat[j]+4)]
     pp[is.na(pp)]=0
     sub=paste("param = ", paste(format(pp,digits=3), collapse=",  "))
    }
    maxi=max(v_LO)
    #Print(j,theta,v_LO[,fromP[j]:toP[j]])
    matplot(theta,v_LO[,fromP[j]:toP[j]], main=title, sub=sub
         , xlim=c(min(theta),max(theta)), ylim=c(0,maxi)
         , type="l", ylab="weight")
   }


  }
 }

 res=named_list( theta, info, fromP, toP, TRF, TRF_LO, icrf, dicrf
               , info_LO, info_LOW, info_item_LO, info_item_LOW
               ,  w_LO, v_LO, A=A )


} # end of info_func





paramS0=paramS1[1:2,]
paramS0$type=c("Bn","Bn3")
paramS0$name=c("Q1n","Q2n")
paramS0=rbind(paramS1[1:2,],paramS0)
paramS0$p1=1
theta=seq(-4,4,length=121)
resInfo=info_func( paramS0, theta, plot=2, print=1 )





comments('

resInfo=info_func( paramS1, npoints=15, plot=1, print=1 )
resInfo=info_func( paramS1, weight=weightS11, npoints=15, plot=1, print=0 )
resInfo=info_func( paramS1, weight=weightS12, npoints=15, plot=1, print=0 )


resInfo=info_func( paramS2, plot=3, print=3, npoints=151 )



param=paramB2[c(1:3,7:9,13:15),]
param$p3=0.5
resInfo=info_func( param, plot=3, print=3, npoints=151 )


param=paramS2[c(1,3,5,7),]
param$p1=1
resInfo=info_func( param, plot=3, print=3, npoints=51 )


param=paramS2[c(1,3,5,7),]
param$p1=1
weight=weightS2[c(1,3,5,7),]
resInfo=info_func( param, weight, plot=3, print=3, npoints=51 )

param=paramS2[c(1,3,5,7),]
param$p1=1
weight=weightS2[c(1,3,5,7),]
weight[4,7]=4
resInfo=info_func( param, weight, plot=3, print=3, npoints=151 )


')








