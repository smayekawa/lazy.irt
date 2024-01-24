#' Calculation of Derivative of Item Response Function
#'
#' @param param Item Parameter Data Frame
#' @param theta Discrete theta values
#' @param weight Weight data frame
#' @param zero = 0 to exclude the xzero-th category from output
#' @param smallP Minimum value of probability
#' @param thmin Minimum value of discrete thata value
#' @param thmax Maximum value of discrete thata value
#' @param npoints # of discrete points for theta
#' @param DinP = 1 to include D=1.7 in logistic function
#' @param numderiv = 1 to use numerical derivatives
#' @param eps eps for JacobianMat
#' @param log = 1 to calculate log derivatives
#' @param print = 1 to print result
#' @param debug = 1 to print intemediate result
#' @param plot = 1 to plot result
#'
#' @return A list of \cr
#'   dICRF, dIRF, dTRF, fromP, toP=toP, vecv, minscore_i, maxscore_i,
#'   minscore_t, maxscore_tt, log
#'\cr
#'   where \cr
#'     dICRF  npoints x sum(ncat) \cr
#'     dIRF   npoints x nitems   weighted by item category weight \cr
#'     dTRF   npoints x 1        weighted by item category weight \cr
#'                                          and item weight \cr
#'     fromP, toP    location of each item category in dICRF \cr
#'     vectorize category weights \cr
#'     minscore_i mimimum score of each item \cr
#'     maxscore_i maximum score of each item \cr
#'     minscore_t mimimum score of test \cr
#'     maxscore_t maximum score of test
#' \cr
#' Note that when log=1, dICRF etc are the log derivatives, namely,
#' the derivative of log ICRF w.r.t. theta, etc.
#'
#' @examples
#' dirf( paramS1, plot=1 )
#'
#' # compare analytic and numeric derivative
#' res1=dirf( paramS1, print=0, plot=0 )$dICRF
#' res2=dirf( paramS1, print=0, numderiv=1, plot=0 )$dICRF
#' Print(max(abs(res1-res2)))
#'
#' @export
#'

dirf <- function( param, theta=NULL, weight=NULL, zero=1, smallP=1e-9
                  , thmin=-4, thmax=4, npoints=21, DinP=1
                  , numderiv=0, eps=1e-6, log=0
                  , print=1, debug=0, plot=0 ){
 # derivative of ICRF, IRF, and TRF
 # Shin-ichi Maykeawa
 # 120215
 # negative/non integer weight: 120223
 # iname as factor: 120223
 # exclude zero: 120223
 # checkparam: 120229
 # theta=NULL: 120304
 # print param value: 120305
 # 120306
 # type Bxx items: 120919,21
 # smallP: 121022
 # DinP: 121109(London)
 # type P items: 121110(Heathlow Airport)
 # plot: 121120
 # nominal response model: 121123
 # bugfix: 201506013
 # numderiv: 20180111
 # log: 20180128
 # smallP=1e-9: 20180129
 # digits -> fmt: 20180208
 # Bn: 20180208
 # Gn: 20180209
 #
 # Args:
 #
 #   param   data.frame containing item parameters as
 #     name type ncat  a  p1  p2  p3   .....   p[ncat-1]
 #    See the description in read.param.R
 #
 #   weight  data.frame containing item and category weight
 #     name type ncat  w  v0 V1  v2  v3   .....   p[ncat-1]
 #    See the description in read.weight.R
 #
 #
 # Value:  as list
 #
 #   list( dICRF, dIRF, dTRF, fromP, toP )
 #
 #   where
 #     ICRF  npoints x sum(ncat)
 #     IRF   npoints x nitems   weighted by item category weight
 #     TRF   npoints x 1        weighted by item category weight
 #                                          and item weight
 #     fromP, toP    location of each item category in ICRF
 #
 # Needs:
 #   icrfB, icrfBn, icrfG,  icrfPN
 #   dicrfB, dicrfBn, dicrfG,  dicrfPN
 #

 # argument name
 pdfname=as.character(substitute(param))
 wdfname=as.character(substitute(weight))

 # param and weight given as data frame
 param=checkparam( param, "ALL" , "irf" )
 if( is.null(param) ){
  cat("\n\nInput Item Parameter ERROR.\n\n")
  return()
 }
 #  if( !is.null(weight) ){
 #   weight=checkparam( weight, , "irf" )
 #   if( is.null(weight) ){
 #    cat("\n\nInput Item Weight ERROR.\n\n")
 #    return()
 #   }
 #  }

 nitems=nrow(param)
 ncp=ncol(param)
 iname=as.character(param$name)
 param$name=iname
 ncat=param$ncat
 type=param$type
 # Print(param,class(param),nitems, ncp,iname,ncat,type)

 if( is.null(theta) ){
  theta=seq(thmin,thmax,length.out=npoints)
 }
 else{
  if( is.matrix(theta) ) theta=as.vector(theta)
  npoints=length(theta)
  thmin=theta[1]; thmax=theta[npoints]
 }

 # item and category weight
 if( is.null(weight) ){
  # natural weight
  w=matrix(1,nitems); rownames(w)=iname; colnames(w)="w"
  v=matrix(0:(max(ncat)-1),nitems,max(ncat),byrow=1)
  rownames(v)=iname; colnames(v)=paste("v",0:(max(ncat)-1),sep="")
  for( i in 1:nitems ){
   if( ncat[i]+1 <= ncol(v) ) v[i,(ncat[i]+1):ncol(v)]=NA
  }
  weight=cbind(ncat,w,v); rownames(weight)=iname
  colnames(weight)=c("ncat","w",colnames(v))
 }
 else{
  w=matrix(as.numeric(weight$w),nitems,1 )
  rownames(w)=iname;
  colnames(w)="w"
  v=matrix( as.numeric( unlist( weight[,5:ncol(weight)] ) ), nitems )
 }

 maxscore_i=as.matrix(apply(v,1,max,na.rm=1),nitems)
 rownames(maxscore_i)=iname
 maxscore_t=sum(w*maxscore_i)
 minscore_i=as.matrix(apply(v,1,min,na.rm=1),nitems)
 rownames(minscore_i)=iname
 minscore_t=min(w*minscore_i)


 # locations in param of each type of items
 locB=grep("^B[[:digit:]]*$", type)
 locBn=grep("^Bn[[:digit:]]*$", type)
 locG=which(type=="G")
 locGn=which(type=="Gn")
 locP=which(type=="P")
 locPN=which(type=="PN")
 locN=which(type=="N")
 nB=length(locB)
 nBn=length(locBn)
 nG=length(locG)
 nGn=length(locGn)
 nN=length(locN)
 nP=length(locP)
 nPN=length(locPN)

 # locations in ircf-storage in the original order
 dicrfP=matrix(0,npoints,sum(ncat))
 cnP=character(max(ncat))
 toP=cumsum(ncat)
 fromP=toP-ncat+1
 ll=mapply(seq,fromP,toP)   # simplify removed: 120928
 if( is.matrix(ll) ){
  llB=as.vector(ll[,locB])
  llBn=as.vector(ll[,locBn])
  llG=as.vector(ll[,locG])
  llGn=as.vector(ll[,locGn])
  llN=as.vector(ll[,locN])
  llP=as.vector(ll[,locP])
  llPN=as.vector(ll[,locPN])
 }
 else if( is.list(ll) ){
  llB=unlist(ll[locB])
  llBn=unlist(ll[locBn])
  llG=unlist(ll[locG])
  llGn=unlist(ll[locGn])
  llN=unlist(ll[locN])
  llP=unlist(ll[locP])
  llPN=unlist(ll[locPN])
 }

 if( debug > 0 ) Print(locB,locG,locP,locPN, locN,ll,llB,llG,llP,llPN,llN)

 if( numderiv == 0 ){
  # for each item type, calculate icrf and store them in icrfP
  if( nB > 0 ){
   pB=as.matrix( param[locB,4:6,drop=0] )
   pB[is.na(pB)]=0                            # 20150613
   rownames(pB)=iname[locB]
   colnames(pB)=c("a","b","c")
   diccB=dicrfB( pB, theta, smallP=smallP, zero=1 )
   if( debug > 0 ) Print(pB,diccB,fmt="6.3")
   dicrfP[,llB]=diccB;
   cnP[llB]=colnames(diccB)
  }
  if( nBn > 0 ){
   pBn=as.matrix( param[locBn,4:6,drop=0] )
   pBn[is.na(pBn)]=0                            # 20150613
   rownames(pBn)=iname[locBn]
   colnames(pBn)=c("a","b","c")
   diccBn=dicrfBn( pBn, theta, smallP=smallP, zero=1 )
   if( debug > 0 ) Print(pBn,diccBn,fmt="6.3")
   dicrfP[,llBn]=diccBn;
   cnP[llBn]=colnames(diccBn)
  }
  if( nPN > 0 ){
   pP=as.matrix( param[locPN,3:ncp,drop=0] )
   rownames(pP)=iname[locPN]
   res=dicrfPN( pP, theta, smallP=smallP )
   diccP=res$dP
   if( debug > 0 ) Print(pP,diccP,fmt="6.3")
   dicrfP[,llPN]=diccP;
   cnP[llPN]=colnames(diccP)
  }
  if( nP > 0 ){
   pP=as.matrix( param[locP,3:ncp,drop=0] )
   rownames(pP)=iname[locP]
   res=dicrfP( pP, theta, DinP=DinP, smallP=smallP )
   diccP=res$dP
   if( debug > 0 ) Print(pP,diccP,fmt="6.3")
   dicrfP[,llP]=diccP;
   cnP[llP]=colnames(diccP)
  }
  if( nG > 0 ){
   pG=as.matrix( param[locG,3:ncp,drop=0] )
   rownames(pG)=iname[locG]
   res=dicrfG( pG, theta, smallP=smallP )
   diccG=res$dP
   if( debug > 0 ) Print(pG,diccG,fmt="6.3")
   dicrfP[,llG]=diccG;
   cnP[llG]=colnames(diccG)
  }
  if( nGn > 0 ){
   pGn=as.matrix( param[locGn,3:ncp,drop=0] )
   rownames(pGn)=iname[locGn]
   res=dicrfGn( pGn, theta, smallP=smallP )
   diccGn=res$dP
   if( debug > 0 ) Print(pGn,diccG,fmt="6.3")
   dicrfP[,llGn]=diccGn;
   cnP[llGn]=colnames(diccGn)
  }
  if( nN > 0 ){
   pN=as.matrix( param[locN,3:ncp,drop=0] )
   rownames(pN)=iname[locN]
   res=dicrfN( pN, theta, smallP=smallP )
   diccN=res$dP
   if( debug > 0 ) Print(pN,,diccN,fmt="6.3")
   dicrfP[,llN]=diccN;
   cnP[llN]=colnames(diccN)
  }
  colnames(dicrfP)=cnP; rownames(dicrfP)=format(theta,digits=3)
 } # end of analytic derivation
 else{
  # numeric derivatives
  dicrfP=dicrf_num( param=param, theta=theta, log=log, eps=eps )
  cnP=colnames(dicrfP)
 }




 # vector of category weights
 vecv=NULL
 for( i in 1:nitems ){
  vecv=c(vecv,v[i,1:ncat[i]])
 }
 vecv=as.matrix(vecv,,1)
 rownames(vecv)=cnP; colnames(vecv)=""
 # Print(vecv)

 # item response function and test response function
 dicrfPv=dicrfP*(matrix(1,npoints)%*%t(vecv))
 dirf=matrix(0,npoints,nitems)
 for( i in 1:nitems ){
  dirf[,i]=rowSums( dicrfPv[,fromP[i]:toP[i],drop=0] )
 }
 colnames(dirf)=iname; rownames(dirf)=format(theta,digits=3)
 dtrf=dirf%*%w; colnames(dtrf)=""
 # Print(irf,trf)


 if( numderiv !=1  && log ){
  logc="log"
  # cannot pass the weight d.f. when it is null.
  res=irf( param, theta, zero=1, smallP=smallP, print=0 )
  dicrfP=dicrfP/res$ICRF
  icrfPv=res$ICRF*(matrix(1,npoints)%*%t(vecv))
  dicrfPv=dicrfPv/icrfPv
  irf=matrix(0,npoints,nitems)
  for( i in 1:nitems ){
   irf[,i]=rowSums( icrfPv[,fromP[i]:toP[i],drop=0] )
  }
  dirf=dirf/irf
  dtrf=dtrf/res$TRF
 }
 else logc=""


 if( zero != 1 ){
  # remove the 0-th category
  dicrfP=dicrfP[,-fromP]
  fromP=fromP-0:(nitems-1)
  toP=toP-1:nitems
 }

 if( print ){
  cat("\n\nCalculation of the Derivative of ICRF \n")
  param=cbind(param,maxscore_i)
  weight=cbind(weight,maxscore_i)
  colnames(weight)[ncol(weight)]="maxscore_i"
  cat(" parameter data frame name =", pdfname
      , ",  item category weight data frame name =", wdfname,"\n")
  cat(" # of item parameters =",nitems
      , ", # of theta points =", npoints,"\n")
  cat("\n  D in P =", DinP, "\n")
  cat(" numderiv = ", numderiv, "( eps =", eps,")\n")
  cat(" item parameters and max score:  total score =", maxscore_t,"\n")
  print( param )
  cat("\n item and item category weight \n")
  print( weight )
  cat("\n")
  cat("\n dICRF\n")
  print(dicrfP,digits=3)
  cat("\n dIRF\n")
  print(dirf,digits=3)
  cat("\n dTRF\n")
  print(dtrf,digits=3)
  # Print(icrfP, irf, trf, digits=3)
 }

 if( plot ){
  for( j in 1:nitems ){

   # titles
   main=paste("Plot of the Derivative of", logc,"ICRF of "
               ,iname[j]," ( type = ",type[j],", ncat=",ncat[j]," )")
   sub=paste("param = ", paste(format(param[j,4:(ncat[j]+3)],digits=3)
                               , collapse=",  "))
   if( length(grep("^B", param$type[j])) > 0 )
    sub=paste("param = ", paste(format(param[j,4:(ncat[j]+4)],digits=3)
                                , collapse=",  "))

   # set up the plot
   plot(range(theta), c(min(dicrfP),max(dicrfP)), type="n"
        , xlab="theta", ylab="dicrf" )
   colors <- rainbow(ncat[j])
   linetype <- c(1:ncat[j])
   plotchar <- seq(18,18+ncat[j],1)

   # add lines
   for (k in 1:ncat[j]) {
    lines(theta, dicrfP[,fromP[j]+k-1], type="b", lwd=1.5,
          lty=linetype[k], col=colors[k], pch=plotchar[k])
   }

   # add a title and subtitle
   title(main, sub)

   # add a legend
   legend(range(theta)[1], max(dicrfP), (1:ncat[j])-1, cex=0.8, col=colors,
          pch=plotchar, lty=linetype, title="cat")
  }
  if( plot >= 2 ){
   for( j in 1:(nitems) ){
    title=paste("Plot of the Derivative of", logc, "IRF of "
                , iname[j]," ( type = ",type[j],", ncat=",ncat[j]
                , ", score range = [",minscore_i[j],",", maxscore_i[j],"] )")
    sub=paste("param = ", paste(format(param[j,4:(ncat[j]+3)],digits=3)
                                , collapse=",  "), "  (with weights)")
    if( length(grep("^B", param$type[j])) > 0 )
     sub=paste("param = ", paste(format(param[j,4:(ncat[j]+4)],digits=3)
                                 , collapse=",  "), "  (with weights)")
    minP=min(dirf); maxP=max(dirf)
    plot(theta,dirf[,j], main=title, sub=sub
         , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP)
         , type="l", ylab="dIRF")
   }
   title=paste("Plot of the Derivative of", logc
               , "TRF  (# of items =", nitems
               , ", score range = [",minscore_t,",", maxscore_t,"] )")
   minP=min(dtrf); maxP=max(dtrf)
   plot(theta,dtrf, main=title
        , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP)
        , type="l", ylab="dTRF")
  }
 }

 return( list(dICRF=dicrfP, dIRF=dirf, dTRF=dtrf, fromP=fromP, toP=toP
              , vecv=vecv, minscore_i=minscore_i, maxscore_i=maxscore_i
              , minscore_t=minscore_t, maxscore_t=maxscore_t, log=log) )


} # end of dirf




#' Numerical Derivatives of ICRF w.r.t. theta
#'
#' @param param Item Parameter Data Frame
#' @param theta Discrete theta values
#' @param log = 1 to calculate log derivatives.
#' @param second = 1 to calculate second derivatives
#' @param eps eps for JacobianMat
#'
#'
#' @details
#' This function uses lazy.mat::JacobianMat to calculate the first derivative.
#' Except for the second derivatives, the analytic functions are available
#' in \code{dirf}.
#'
#' @return
#' dICRF  length(theta) x sum(ncat) matrix
#'
#'
#' @examples
#' res1 <- dicrf_num( paramS1, -3:3 )
#' res2 <- dirf( paramS1, -3:3, print=0, plot=0 )$dICRF
#' Print(apply(abs(res1-res2),2,max))
#'
#' # second derivative of P
#' res3 <- dicrf_num( paramS1, -3:3, second=1 )
#'
#'
#' \dontrun{
#'
#' thmin <- -4
#' thmax <- 4
#' npoints <- 121
#' theta <- seq(thmin,thmax,length=npoints)
#'
#' param <- paramA1[1:7,]
#'
#'
#' log <- 1
#' res1 <- dicrf_num( param, theta, log=log, second=1 )
#'
#' # all
#' matplot(theta,res1, type="l", main="2nd deriv. of logP")
#'
#' # 0-th category: all negative
#' fromP <- get_range(param$ncat )[,1]
#' matplot(theta,res1[,fromP], type="l", main="2nd deriv. of logP: 0-th cat")
#'
#' # all negative except for the 1st category of 3PLM
#' matplot(theta,res1[,-c(4,8)], type="l"
#'         , main="2nd deriv. of logP: all but 1st of 3PLM")
#'
#' # 1st category of 3PLM
#' matplot(theta,res1[,c(4,8)], type="l", main="2nd deriv. of logP: 1st of 3PLM")
#'
#' }
#'
#' @export
#'
#'

dicrf_num <- function( param, theta, log=0, second=0, eps=1e-6 ){
 # numerical version of dirf
 # Shin-ichi Mayekawa
 # 20180111
 # log: 20180128
 # second: 20180129
 # roxgen2 example: 20180203, 20221206,08
 #

 icrf <- function( x, param, log ){
  res=irf( param, x, plot=0, print=0 )$ICRF
  if( log ) res=log(res)
  return( res )
 } # end of icrf


 dicrf <- function( x, param, log ){
  res=dirf( param, x, log=log, plot=0, print=0 )$dICRF
  return( res )
 } # end of dicrf


 nitem=nrow(param)
 ncolP=sum(param$ncat)
 npoint=length(theta)
 res=matrix(0,npoint,ncolP)
 rownames(res)=format(theta,digits=3)

 if( second ){
  for( i in 1:npoint ){
   # res[i,]=jacobian( icrf, theta[i], param=param )
   res[i,]=JacobianMat( theta[i], dicrf, ..eps..=eps, param=param, log=log )
  }
 }
 else{
  for( i in 1:npoint ){
   # res[i,]=jacobian( icrf, theta[i], param=param )
   res[i,]=JacobianMat( theta[i], icrf, ..eps..=eps, param=param, log=log )
  }
 }
 colnames(res)=colnames(icrf(theta[1],param,log))
 return( res )

} # end of dicrf_num









#' Calculation of the Derivative of Item Response Function of Binary
#' Logistic Items
#'
#' @param param Item Parameter Data Frame
#' @param theta Discrete theta values
#' @param maxZ Maximum value of the argument to logistic function
#' @param zero = 0 to exclude the xzero-th category from output
#' @param smallP Minimum value of probability
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#'
#' @return
#'   dP       length(theta) x nitems   dICRF matrix when zero=0 \cr
#'            length(theta) x 2*nitems dICRF matrix when zero=1
#'
#'
#' @export
#'
dicrfB <- function( param, theta, maxZ=700, zero=0, smallP=0
                    , print=0, plot=0 ){
 # derivative of Binary Logistic ICRF
 # Shin-ichi Mayekawa
 # 120214,15
 # bugfix: 120217
 # name as factor: 120224
 # checkparam: 120224,29
 # smallP: 121022
 # when c=NA: 20180208
 # digits -> fmt: 20180208
 #
 # Args:
 #  param    nitems x 3 item parameters matrix or data frame
 #  theta    npoints x 1 theta values
 #  zero     = 1 to include 1-P as the 0-th category probability
 #
 # Values:
 #   dP       length(theta) x nitems   ICRF matrix when zero=0
 #            length(theta) x 2*nitems ICRF matrix when zero=1
 #
 # Needs:
 #   icrfB, checkparam
 #

 # argument name
 paramname=as.character(substitute(param))

 isdf=0
 if( is.data.frame(param) ){
  # param and weight given as data frame
  isdf=1
  param=checkparam( param, "B", "dicrfB" )
  if( is.null(param) ){
   cat("\n\nInput Item Parameter ERROR.\n\n")
   return()
  }
 }

 nitems=nrow(param)
 iname=rownames(param)
 if( length(which(colnames(param) == "ncat")) > 0 )
  param=param[,-which(colnames(param) == "ncat")]
 a=param[,1]; c=param[,3]
 c[is.na(c)]=0

 if( is.matrix(theta) ) theta=as.vector(theta)
 npoints=length(theta)

 # calculate icrf
 P=icrfB( param, theta, smallP=smallP, maxZ=maxZ, zero=0, print=0, plot=0 )
 cc=matrix(c,npoints,nitems,byrow=1)
 dP=1.7*matrix(a,npoints,nitems,byrow=1)*(1-P)*(P-cc)/(1-cc)
 colnames(dP)=rownames(param)
 if( zero == 1 ){
  P0=matrix(0,npoints,2*nitems)
  P0[,2*(1:nitems)]=dP
  P0[,2*(1:nitems)-1]=-dP
  catname=outer(paste(iname,"_",sep=""),0:1,paste,sep="")
  colnames(P0)=matrix(t(catname),1)
  dP=P0
 }
 rownames(dP)=format(theta,fmt="6.3")
 if( print > 0 ){
  cat("\nDerivative of IRF of Binary Logistic Items \n ")
  if( isdf )
   cat(" item parameter data frame name =", paramname,"\n")
  else
   cat(" item parameter matrix name =", paramname,"\n")
  cat("  # of item parameters =",nitems,
      ",  # of theta points =", npoints,
      ",  zero category =",zero, "\n")
  Print(param)
  Print(dP, fmt="6.3")
 }
 if( plot > 0 ){
  matplot(theta,dP, type = "l", ylab="IRF"
          , xlim=c(min(theta),max(theta)), ylim=
           , main="Plot of the Derivatives of ICRF of Binary Items")
 }
 return( dP )

} # end of dicrfB



#' Calculation of the Derivative of Item Response Function of Binary
#' Normal CDF Items
#'
#' @param param Item Parameter Data Frame
#' @param theta Discrete theta values
#' @param maxZ Maximum value of the argument to logistic function
#' @param zero = 0 to exclude the xzero-th category from output
#' @param smallP Minimum value of probability
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#'
#' @return
#'   dP       length(theta) x nitems   dICRF matrix when zero=0 \cr
#'            length(theta) x 2*nitems dICRF matrix when zero=1
#'
#'
#' @export
#'
dicrfBn <- function( param, theta, maxZ=700, zero=0, smallP=0
                    , print=0, plot=0 ){
 # derivative of Binary Normal CDF ICRF
 # Shin-ichi Mayekawa
 # dicrfB: 20120214-20180209
 # dicrfB modified: 20180200
 #
 # Args:
 #  param    nitems x 3 item parameters matrix or data frame
 #  theta    npoints x 1 theta values
 #  zero     = 1 to include 1-P as the 0-th category probability
 #
 # Values:
 #   dP       length(theta) x nitems   ICRF matrix when zero=0
 #            length(theta) x 2*nitems ICRF matrix when zero=1
 #
 # Needs:
 #   icrfBn, checkparam
 #

 # argument name
 paramname=as.character(substitute(param))

 isdf=0
 if( is.data.frame(param) ){
  # param and weight given as data frame
  isdf=1
  param=checkparam( param, "Bn", "dicrfBn" )
  if( is.null(param) ){
   cat("\n\nInput Item Parameter ERROR.\n\n")
   return()
  }
 }

 nitems=nrow(param)
 iname=rownames(param)
 if( length(which(colnames(param) == "ncat")) > 0 )
  param=param[,-which(colnames(param) == "ncat")]
 a=param[,1]; b=param[,2]; c=param[,3]
 c[is.na(c)]=0

 if( is.matrix(theta) ) theta=as.vector(theta)
 npoints=length(theta)

 # calculate icrf
 Z=-(a*outer(b,theta,"-"))
 Z[which(Z > maxZ)]=maxZ
 Z[which(Z < -maxZ)]=-maxZ
 dP=t( a*(1-c)*dnorm(Z) )
 colnames(dP)=rownames(param)
 if( zero == 1 ){
  P0=matrix(0,npoints,2*nitems)
  P0[,2*(1:nitems)]=dP
  P0[,2*(1:nitems)-1]=-dP
  catname=outer(paste(iname,"_",sep=""),0:1,paste,sep="")
  colnames(P0)=matrix(t(catname),1)
  dP=P0
 }
 rownames(dP)=format(theta,digits=3)
 if( print > 0 ){
  cat("\nDerivative of IRF of Binary Normal CDF Items \n ")
  if( isdf )
   cat(" item parameter data frame name =", paramname,"\n")
  else
   cat(" item parameter matrix name =", paramname,"\n")
  cat("  # of item parameters =",nitems,
      ",  # of theta points =", npoints,
      ",  zero category =",zero, "\n")
  Print(param)
  Print(dP, fmt=".3")
 }
 if( plot > 0 ){
  matplot(theta,dP, type = "l", ylab="IRF"
          , xlim=c(min(theta),max(theta)), ylim=
           , main="Plot of the Derivatives of ICRF of Binary Items")
 }
 return( dP )

} # end of dicrfBn





#' Calculation of the Derivative of Item Response Function of Graded
#' Response Items
#'
#' @param param Item Parameter Data Frame
#' @param theta Discrete theta values
#' @param smallP Minimum value of probability
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#' @param debug = 1 to print intemediate result
#'
#' @return A list of \cr
#'   dP       length(theta) x sum(ncat) dICRF matrix \cr
#'   dPP      length(theta) x sum(ncat) dICBRF matrix \cr
#'   dPt      length(theta) x nitems derivative of item response function \cr
#'   fromP, toP    indexes of each ICRF in P \cr
#'   fromPP, toPP  indexes of each ICBRF in PP
#'
#' @export
#'
dicrfG <- function( param, theta, smallP=0, print=0, plot=0, debug=0 ){
 # derivative of GRM ICRF
 # Shin-ichi Mayekawa
 # 120214
 # checkparam: 120224,29
 # smallP: 121022
 #
 # Args:
 #  param    nitems x max(ncat) matrix of item parameters or data frame:
 #    An Example of Item Parameter datastet:  maxncat = 5
 #    Note that the subscript of the b-parameter ranges between
 #         1 and ncat-1 (= maxscore)  (We set c0=0.)
 #
 #          ncat     A            b1     b2     b3     b4
 #
 #    Item1  2     1.0           0.0     .      .      .
 #    Item2  2     1.0           0.0     .      .      .
 #    Item3  3     1.5          -1.0    0.0     .      .
 #    Item4  4     1.5          -1.0    0.0    1.0     .
 #    Item5  5     1.5          -1.0    0.0    1.0    2.0
 #
 #  theta    npoints x 1 theta values
 #
 # Values: as list
 #   P             length(theta) x sum(ncat) ICRF matrix
 #   PP            length(theta) x sum(ncat) ICBRF matrix
 #   Pt            length(theta) x nitems item response function
 #   fromP, toP    indexes of each ICRF in P
 #   fromPP, toPP  indexes of each ICBRF in PP
 #
 #
 # Needs:
 #   icrfB, checkparam
 #
 #


 # argument name
 paramname=as.character(substitute(param))

 isdf=0
 if( is.data.frame(param) ){
  # param and weight given as data frame
  isdf=1
  param=checkparam( param, "G", "dicrfG" )
  if( is.null(param) ){
   cat("\n\nInput Item Parameter ERROR.\n\n")
   return()
  }
 }


 # const
 if( is.matrix(theta) ) theta=as.vector(theta)
 npoints=length(theta)
 thetac=format(theta,digits=3)

 nitems=nrow(param)
 iname=rownames(param)
 ncat=param[,1]
 ncat1=ncat-1
 a=param[,2]
 b=param[,3:(1+max(ncat)),drop=0]

 # locations of dicrf and icbrf of each item in P and PP
 toP=cumsum(ncat)
 fromP=toP-ncat+1
 toPP=cumsum(ncat+1)
 fromPP=toPP-ncat
 if( debug ) Print(fromPP,toPP,fromP,toP)

 # names of icrf and icbrf
 PPname=character(sum(ncat+1))
 Pname=character(sum(ncat))
 catname=outer(paste(iname,"_",sep=""),0:max(ncat),paste,sep="")
 for( j in 1:nitems ){
  PPname[fromPP[j]:toPP[j]]=catname[j,1:(ncat[j]+1)]
  Pname[fromP[j]:toP[j]]=catname[j,1:ncat[j]]
 }
 if( debug ) Print(nitems,catname)
 rm(catname)
 if( debug ) Print(PPname, Pname)


 # main body
  if( debug ) Print(a,b,ncat,ncat1)
 dPP=matrix( 0,npoints, sum(ncat+1))
 dP=matrix( 0,npoints, sum(ncat))
 dPt=matrix(0,npoints,nitems)
 rownames(dPP)=thetac; colnames(dPP)=PPname
 rownames(dP)=thetac;  colnames(dP)=Pname
 rownames(dPt)=thetac;  colnames(dPt)=iname

 if( debug ) print1=1
 else print1=0

 # for each item
 for( j in 1:nitems ){

  # derivative of binary logistic icrf
  paramj=cbind( rep(a[j],ncat1[j]), b[j,1:ncat1[j]],rep(0,ncat1[j]) )
  dPPj=dicrfB( paramj, theta, smallP=smallP, print=print1 )

  # icbrf
  dPP[,fromPP[j]]=0;
  dPP[,(fromPP[j]+1):(toPP[j]-1)]=dPPj
  dPP[,toPP[j]]=0;

  # dicrf as difference of icbrf
  #    pj=PP[,fromPP[j]:toPP[j]-1]-PP[,fromPP[j]+1:toPP[j]];
  dP[,fromP[j]:toP[j]]=
      dPP[,fromPP[j]:(toPP[j]-1)]-dPP[,(fromPP[j]+1):toPP[j]];

  # dicrf
  dPt[,j]=rowSums(dP[,fromP[j]:toP[j]]%*%diag(c(0:(ncat[j]-1))))

 } # end of j loop

 if( print > 0 ){
  cat("\nICRF of Graded Response Items \n ")
  if( isdf )
   cat(" item parameter data frame name =", paramname,"\n")
  else
   cat(" item parameter matrix name =", paramname,"\n")
  cat("  # of item parameters =",nitems,
      ",  # of theta points =", npoints,"\n")
  Print(param)
  cat("\nDerivative of Item Category Response Functions (ICRF)\n")
  Print(dP, digits=2)
  cat("\nDerivative ofItem Response Functions (IRF)"
       ," with natural category weights\n")
  Print(dPt, digits=2)
  if( print >= 2 ){
   cat("\nDerivative of Item Category Boundary Response Functions (ICBRF)\n")
   Print(dPP,digits=2)
  }
 }
 if( plot ){
  for( j in 1:(nitems) ){
   title=paste("Plot of the Derivative of ICRF of "
             , iname[j]," (ncat=",ncat[j],")")
   minP=min(dP); maxP=max(dP)
   for( k in 1:ncat1[j] ){
    plot(theta,dP[,fromP[j]+k-1]
        , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP), type="l", ylab="")
    par(new=1)
   }
   plot(theta,dP[,toP[j]], main=title, sub="Graded Responce Model"
       , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP)
       , type="l", ylab="dICRF")
   par(new=0)
  }
 }
 if( plot >= 2 ){
  for( j in 1:(nitems) ){
   title=paste("Plot of the Derivative of IRF"
               ,"(with natural category weight) of "
               , iname[j]," (ncat=",ncat[j],")")
   plot(theta,dPt[,j], main=title, sub="Graded Responce Model"
        , xlim=c(min(theta),max(theta)), ylim=, type="l", ylab="dIRF")
   par(new=0)
  }
 }
 if( plot >= 3 ){
  for( j in 1:(nitems) ){
   minP=min(dPP); maxP=max(dPP)
   title=paste("Plot of the Derivative of ICBRF of ",iname[j]
             , " (ncat=",ncat[j],")")
   for( k in 1:ncat[j] ){
    plot(theta,dPP[,fromPP[j]+k-1]
         , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP), type="l", ylab="")
    par(new=1)
   }
   plot(theta,dPP[,toPP[j]], main=title, sub="Graded Response Model"
        , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP)
        , type="l", ylab="dICBRF")
   par(new=0)
  }
 }

 return( list(dP=dP, dPP=dPP, dPt=dPt
              , fromP=fromP, toP=toP, fromPP=fromPP, toPP=toPP) )

} # end of dicrfG



#' Calculation of the Derivative of Item Response Function of Graded
#' Response Items with Normal Ogive MOdel
#'
#' @param param Item Parameter Data Frame
#' @param theta Discrete theta values
#' @param smallP Minimum value of probability
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#' @param debug = 1 to print intemediate result
#'
#' @return A list of \cr
#'   dP       length(theta) x sum(ncat) dICRF matrix \cr
#'   dPP      length(theta) x sum(ncat) dICBRF matrix \cr
#'   dPt      length(theta) x nitems derivative of item response function \cr
#'   fromP, toP    indexes of each ICRF in P \cr
#'   fromPP, toPP  indexes of each ICBRF in PP
#'
#' @export
#'

dicrfGn <- function( param, theta, smallP=smallP, print=0, plot=0, debug=0 ){
 # derivative of GRM ICRF
 # Shin-ichi Mayekawa
 # dicrfG: 20120214-20120229
 # dicrfG modified: 20180209
 #
 # Args:
 #  param    nitems x max(ncat) matrix of item parameters or data frame:
 #    An Example of Item Parameter datastet:  maxncat = 5
 #    Note that the subscript of the b-parameter ranges between
 #         1 and ncat-1 (= maxscore)  (We set c0=0.)
 #
 #          ncat     A            b1     b2     b3     b4
 #
 #    Item1  2     1.0           0.0     .      .      .
 #    Item2  2     1.0           0.0     .      .      .
 #    Item3  3     1.5          -1.0    0.0     .      .
 #    Item4  4     1.5          -1.0    0.0    1.0     .
 #    Item5  5     1.5          -1.0    0.0    1.0    2.0
 #
 #  theta    npoints x 1 theta values
 #
 # Values: as list
 #   P             length(theta) x sum(ncat) ICRF matrix
 #   PP            length(theta) x sum(ncat) ICBRF matrix
 #   Pt            length(theta) x nitems item response function
 #   fromP, toP    indexes of each ICRF in P
 #   fromPP, toPP  indexes of each ICBRF in PP
 #
 #
 # Needs:
 #   icrfB, checkparam
 #
 #


 # argument name
 paramname=as.character(substitute(param))

 isdf=0
 if( is.data.frame(param) ){
  # param and weight given as data frame
  isdf=1
  param=checkparam( param, "Gn", "dicrfGn" )
  if( is.null(param) ){
   cat("\n\nInput Item Parameter ERROR.\n\n")
   return()
  }
 }


 # const
 if( is.matrix(theta) ) theta=as.vector(theta)
 npoints=length(theta)
 thetac=format(theta,digits=3)

 nitems=nrow(param)
 iname=rownames(param)
 ncat=param[,1]
 ncat1=ncat-1
 a=param[,2]
 b=param[,3:(1+max(ncat)),drop=0]

 # locations of dicrf and icbrf of each item in P and PP
 toP=cumsum(ncat)
 fromP=toP-ncat+1
 toPP=cumsum(ncat+1)
 fromPP=toPP-ncat
 if( debug ) Print(fromPP,toPP,fromP,toP)

 # names of icrf and icbrf
 PPname=character(sum(ncat+1))
 Pname=character(sum(ncat))
 catname=outer(paste(iname,"_",sep=""),0:max(ncat),paste,sep="")
 for( j in 1:nitems ){
  PPname[fromPP[j]:toPP[j]]=catname[j,1:(ncat[j]+1)]
  Pname[fromP[j]:toP[j]]=catname[j,1:ncat[j]]
 }
 if( debug ) Print(nitems,catname)
 rm(catname)
 if( debug ) Print(PPname, Pname)


 # main body
 if( debug ) Print(a,b,ncat,ncat1)
 dPP=matrix( 0,npoints, sum(ncat+1))
 dP=matrix( 0,npoints, sum(ncat))
 dPt=matrix(0,npoints,nitems)
 rownames(dPP)=thetac; colnames(dPP)=PPname
 rownames(dP)=thetac;  colnames(dP)=Pname
 rownames(dPt)=thetac;  colnames(dPt)=iname

 if( debug ) print1=1
 else print1=0

 # for each item
 for( j in 1:nitems ){

  # derivative of binary logistic icrf
  paramj=cbind( rep(a[j],ncat1[j]), b[j,1:ncat1[j]],rep(0,ncat1[j]) )
  dPPj=dicrfBn( paramj, theta, smallP=smallP, print=print1 )

  # icbrf
  dPP[,fromPP[j]]=0;
  dPP[,(fromPP[j]+1):(toPP[j]-1)]=dPPj
  dPP[,toPP[j]]=0;

  # dicrf as difference of icbrf
  #    pj=PP[,fromPP[j]:toPP[j]-1]-PP[,fromPP[j]+1:toPP[j]];
  dP[,fromP[j]:toP[j]]=
   dPP[,fromPP[j]:(toPP[j]-1)]-dPP[,(fromPP[j]+1):toPP[j]]

  # dicrf
  dPt[,j]=rowSums(dP[,fromP[j]:toP[j]]%*%diag(c(0:(ncat[j]-1))))

 } # end of j loop

 if( print > 0 ){
  cat("\nDerivative of ICRF of Graded Response (2PNM) Items \n ")
  if( isdf )
   cat(" item parameter data frame name =", paramname,"\n")
  else
   cat(" item parameter matrix name =", paramname,"\n")
  cat("  # of item parameters =",nitems,
      ",  # of theta points =", npoints,"\n")
  Print(param)
  cat("\nDerivative of Item Category Response Functions (ICRF)\n")
  Print(dP, fmt="6.2")
  cat("\nDerivative ofItem Response Functions (IRF)"
      ," with natural category weights\n")
  Print(dPt, fmt="6.2")
  if( print >= 2 ){
   cat("\nDerivative of Item Category Boundary Response Functions (ICBRF)\n")
   Print(dPP, fmt="6.2")
  }
 }
 if( plot ){
  for( j in 1:(nitems) ){
   title=paste("Plot of the Derivative of ICRF of "
               , iname[j]," (ncat=",ncat[j],")")
   minP=min(dP); maxP=max(dP)
   for( k in 1:ncat1[j] ){
    plot(theta,dP[,fromP[j]+k-1], type="l", ylab=""
        , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP) )
    par(new=1)
   }
   plot(theta,dP[,toP[j]], main=title, sub="Graded Responce Model (2PNM)"
        , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP)
        , type="l", ylab="dICRF")
   par(new=0)
  }
 }
 if( plot >= 2 ){
  for( j in 1:(nitems) ){
   title=paste("Plot of the Derivative of IRF"
               ,"(with natural category weight) of "
               , iname[j]," (ncat=",ncat[j],")")
   plot(theta,dPt[,j], main=title
        , sub="Graded Responce Model (2PNM)", ylab="dIRF"
        , xlim=c(min(theta),max(theta)), ylim=, type="l" )
   par(new=0)
  }
 }
 if( plot >= 3 ){
  for( j in 1:(nitems) ){
   minP=min(dPP); maxP=max(dPP)
   title=paste("Plot of the Derivative of ICBRF of ",iname[j]
               , " (ncat=",ncat[j],")")
   for( k in 1:ncat[j] ){
    plot(theta,dPP[,fromPP[j]+k-1], type="l", ylab=""
        , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP) )
    par(new=1)
   }
   plot(theta,dPP[,toPP[j]], main=title, sub="Graded Response Model (2PNM)"
        , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP)
        , type="l", ylab="dICBRF")
   par(new=0)
  }
 }

 return( list(dP=dP, dPP=dPP, dPt=dPt
              , fromP=fromP, toP=toP, fromPP=fromPP, toPP=toPP) )

} # end of dicrfGn






#' Calculation of the Derivative of Item Response Function of Partial
#' Credit Items in Nominal Format
#'
#' @param paramPN Item Parameter Data Frame for Partial Credit Model
#' in Nominal Format
#' @param theta Discrete theta values
#' @param smallP Minimum value of probability
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#' @param debug = 1 to print intemediate result
#'
#' @return A list of \cr
#'   dP       length(theta) x sum(ncat) dICRF matrix \cr
#'   dPt      length(theta) x nitems derivative of item response function \cr
#'   fromP, toP    indexes of each ICRF in P \cr
#'   fromPP, toPP  indexes of each ICBRF in PP
#'
#'
#' @export
#'
dicrfPN <- function( paramPN, theta, smallP=0, print=0, plot=0, debug=0 ){
 # derivative of Generalized Partial Credit ICRF in Nominal format
 # Shin-ichi Mayekawa
 # 120215
 # checkparam: 120224
 # converted to b-type paramters input: 120301
 # b-type paramters corrected: 120307
 # smallP: 121022
 #
 #
 # Args:
 #  paramPN  nitems x max(ncat) matrix of item parameters or data frame
 #    An Example of Item Parameter datastet:  maxncat = 5
 #    Note that the subscript of the b-parameter ranges between
 #         1 and ncat-1 (= maxscore)  (We set c0=0.)
 #
 #          ncat     A            b1     b2     b3     b4
 #
 #    Item1  2     1.0           0.0     .      .      .
 #    Item2  2     1.0           0.0     .      .      .
 #    Item3  3     1.5          -1.0    0.0     .      .
 #    Item4  4     1.5          -1.0    0.0    1.0     .
 #    Item5  5     1.5          -1.0    0.0    1.0    2.0
 #
 #    ICRF of category k of item j at theta, namely, p_{jk}(theta)
 #    is defined as
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0}^{ncat[j]-1} Ez_{jk}(theta)
 #      where
 #     Ez_{jk}(theta) = exp( z_{jk}(theta) )
 #       and
 #     z_{jk}(theta) = a_j k (theta - b_{k}), k=0,1, ..., ncat[j]-1
 #    with b_{j0}=0.
 #
 #
 #  theta    npoints x 1 theta values
 #
 # Values: as list
 #   P       length(theta) x sum(ncat) icrf matrix
 #   PP      length(theta) x sum(ncat) icbf  matrix
 #   Pt      length(theta) x nitems item response function
 #
 #

 # argument name
 paramname=as.character(substitute(paramPN))

 isdf=0
 if( is.data.frame(paramPN) ){
  # param and weight given as data frame
  isdf=1
  paramPN=checkparam( paramPN, "PN", "dicrfPN" )
  if( is.null(paramPN) ){
   cat("\n\nInput Item Parameter ERROR.\n\n")
   return()
  }
 }

 # convert to c-type params (120307)
 paramPN0=paramPN
 paramPN0[,2]=1.7*paramPN[,2]
 colnames(paramPN0)=c("ncat","a",paste("c",1:(ncol(paramPN0)-2),sep=""))
 for( k in 3:ncol(paramPN0) ){
  paramPN0[,k]=-1.7*paramPN[,k]*paramPN[,2]*(k-2)
 }

 # const
 if( is.matrix(theta) ) theta=as.vector(theta)
 npoints=length(theta)
 thetac=format(theta,digits=3)
 nitems=nrow(paramPN0)
 iname=rownames(paramPN0)
 ncat=paramPN0[,1]
 ncat1=ncat-1
 a=paramPN0[,2,drop=0]
 c=paramPN0[,3:(1+max(ncat)),drop=0]

 # icrf
 temp=icrfPN0( paramPN0, theta, smallP=smallP )
 P=temp$P
 fromP=temp$fromP
 toP=temp$toP
 Pname=colnames(P)
 if( debug ) Print(P,fromP,toP)

 # sum_{k=0}^m-j  k P_{kj}
 vecv=matrix( as.vector( unlist( mapply(seq,0,ncat1) ) ),,1 )
 sumP=matrix(0,npoints,nitems)
 for( j in 1:nitems ){
  sumP[,j]=P[,fromP[j]:toP[j]]%*%vecv[fromP[j]:toP[j],1]
 }
 if( debug ) Print(P,vecv,sumP)
 # d
 dP=P; dPt=matrix(0,npoints,nitems)
 for( j in 1:nitems ){
  dP[,fromP[j]:toP[j]]=a[j]*P[,fromP[j]:toP[j]] *
          ( matrix(1,npoints)%*%t(vecv[fromP[j]:toP[j],1]) - sumP[,j] )
  # dicrf
  dPt[,j]=dP[,fromP[j]:toP[j]]%*%vecv[fromP[j]:toP[j],1]
 }


 if( print > 0 ){
  cat("\nDerivative of ICRF of Generalized Partial Credit Items \n ")
  if( isdf )
   cat(" item parameter data frame name =", paramname,"\n")
  else
   cat(" item parameter matrix name =", paramname,"\n")
  cat("  # of item parameters =",nitems,
      ",  # of theta points =", npoints,"\n")
  cat(" Item Parameters (type b)\n")
  Print(paramPN, digits=3)
  cat("\nDerivative of Item Category Response Functions (ICRF)\n")
  Print(dP, digits=2)
  cat("\nDerivative of Item Response Functions (IRF)"
     ," with natural category weights\n")
  Print(dPt, digits=2)
 }

 if( plot ){
  for( j in 1:(nitems) ){
   title=paste("Plot of the Derivative of ICRF of "
             , iname[j]," (ncat=",ncat[j],")")
   minP=min(dP); maxP=max(dP)
   for( k in 1:ncat1[j] ){
    plot(theta,dP[,fromP[j]+k-1]
         , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP), type="l", ylab="")
    par(new=1)
   }
   plot(theta,dP[,toP[j]], main=title, sub="Generalized Partial Credit Model"
        , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP)
        , type="l", ylab="dICRF")
   par(new=0)
  }
  if( plot >= 2 ){
   for( j in 1:(nitems) ){
    title=paste("Plot of the Derivative of IRF"
               , " (with natural category weight) of "
                , iname[j]," (ncat=",ncat[j],")")
    plot(theta,dPt[,j], main=title, sub="Generalized Partial Credit Model"
         , xlim=c(min(theta),max(theta)), ylim=, type="l", ylab="dIRF")
    par(new=0)
   }
  }
 }


 return( list(dP=dP, dPt=dPt
              , fromP=fromP, toP=toP) )

} # end of dicrfPN




#' Calculation of the Derivative of Item Response Function of Partial
#' Credit Items
#'
#' @param paramP Item Parameter Data Frame of Partial Credit Items
#' @param theta Discrete theta values
#' @param DinP = 1 to include D=1.7 in logistic function
#' @param smallP Minimum value of probability
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#' @param debug = 1 to print intemediate result
#'
#' @return A list of \cr
#'   dP       length(theta) x sum(ncat) dICRF matrix \cr
#'   dPt      length(theta) x nitems derivative of item response function \cr
#'   fromP, toP    indexes of each ICRF in P \cr
#'   fromPP, toPP  indexes of each ICBRF in PP
#'
#' @export
#'
dicrfP <- function( paramP, theta, DinP=1, smallP=0, print=0
                  , plot=0, debug=0 ){
 # derivative of Generalized Partial Credit ICRF in Nominal format
 # Shin-ichi Mayekawa
 # converted from dicrfPN: 121016
 # smallP: 121022
 # DinP: 121109(London), 20
 # bugfix: 121120
 #
 #
 # Args:
 #  paramP  nitems x max(ncat) matrix of item parameters or data frame
 #    An Example of Item Parameter datastet:  maxncat = 5
 #    Note that the subscript of the b-parameter ranges between
 #         1 and ncat-1 (= maxscore)  (We set c0=0.)
 #
 #          ncat     A            b1     b2     b3     b4
 #
 #    Item3  3     1.5          -1.0    0.0     .      .
 #    Item4  4     1.5          -1.0    0.0    1.0     .
 #    Item5  5     1.5          -1.0    0.0    1.0    2.0
 #
 #    ICRF of category k of item j at theta, namely, p_{jk}(theta)
 #    is defined as
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0}^{ncat[j]-1} Ez_{jk}(theta)
 #      where
 #     Ez_{jk}(theta) = exp( z_{jk}(theta) )
 #       and
 #     z_{jk}(theta) = a_j k (theta - b_{k}), k=0,1, ..., ncat[j]-1
 #    with b_{j0}=0.
 #
 #
 #  theta    npoints x 1 theta values
 #
 # Values: as list
 #   P       length(theta) x sum(ncat) icrf matrix
 #   PP      length(theta) x sum(ncat) icbf  matrix
 #   Pt      length(theta) x nitems item response function
 #
 #
 # Needs:
 #  icrfPN
 #

 # argument name
 paramname=as.character(substitute(paramP))

 isdf=0
 if( is.data.frame(paramP) ){
  # param and weight given as data frame
  isdf=1
  paramP=checkparam( paramP, "P", "dicrfP" )
  if( is.null(paramP) ){
   cat("\n\nInput Item Parameter ERROR.\n\n")
   return()
  }
 }

 # convert to c-type params (121016)
 paramP0=paramP
 if( DinP == 1 ) paramP0[,2]=1.7*paramP0[,2]
 bs=paramP[,3:ncol(paramP),drop=0]
 c=t( apply(bs,1,cumsum) )
 paramP0[,3:ncol(paramP)]=-paramP0[,2]*c
 colnames(paramP0)=c("ncat","a",paste("c",1:(ncol(paramP0)-2),sep=""))

 # convert to type="PN" form: we need this because we use icrfPN to get ICRF.
 paramPN=paramP
 # b=t( apply( bs,1,cumsum) )%*%diag(1/(1:ncol(bs))) # does not work if NA.
 b=t( apply( bs,1,cumsum)*(1/(1:ncol(bs))) )
 paramPN[,2]=paramPN[,2]
 if( DinP == 0 ) paramPN[,2]=paramPN[,2]/1.7
 paramPN[,3:ncol(paramPN)]=b

 # const
 if( is.matrix(theta) ) theta=as.vector(theta)
 npoints=length(theta)
 thetac=format(theta,digits=3)

 nitems=nrow(paramP0)
 iname=rownames(paramP0)
 ncat=paramP0[,1]
 ncat1=ncat-1
 a=paramP0[,2,drop=0]
 c=paramP0[,3:(1+max(ncat)),drop=0]

 # icrf
 temp=icrfPN( paramPN, theta, smallP=smallP )
 P=temp$P
 fromP=temp$fromP
 toP=temp$toP
 Pname=colnames(P)
 if( debug ) Print(P,fromP,toP)

 # sum_{k=0}^m-j  k P_{kj}
 vecv=matrix( as.vector( unlist( mapply(seq,0,ncat1) ) ),,1 )
 sumP=matrix(0,npoints,nitems)
 for( j in 1:nitems ){
  sumP[,j]=P[,fromP[j]:toP[j]]%*%vecv[fromP[j]:toP[j],1]
 }
 if( debug ) Print(P,vecv,sumP)
 # d
 dP=P; dPt=matrix(0,npoints,nitems)
 for( j in 1:nitems ){
  dP[,fromP[j]:toP[j]]=a[j]*P[,fromP[j]:toP[j]] *
   ( matrix(1,npoints)%*%t(vecv[fromP[j]:toP[j],1]) - sumP[,j] )
  # dicrf
  dPt[,j]=dP[,fromP[j]:toP[j]]%*%vecv[fromP[j]:toP[j],1]
 }


 if( print > 0 ){
  cat("\nDerivative of ICRF of Generalized Partial Credit Items \n ")
  if( isdf )
   cat(" item parameter data frame name =", paramname,"\n")
  else
   cat(" item parameter matrix name =", paramname,"\n")
  cat("  # of item parameters =",nitems,
      ",  # of theta points =", npoints,"\n")
  cat("\n item parameters with D in P =", DinP, "\n")
  Print(paramP, digits=3)
  cat("\nDerivative of Item Category Response Functions (ICRF)\n")
  Print(dP, digits=2)
  cat("\nDerivative of Item Response Functions (IRF)"
      ," with natural category weights\n")
  Print(dPt, digits=2)
 }

 if( plot ){
  for( j in 1:(nitems) ){
   title=paste("Plot of the Derivative of ICRF of "
               , iname[j]," (ncat=",ncat[j],")")
   minP=min(dP); maxP=max(dP)
   for( k in 1:ncat1[j] ){
    plot(theta,dP[,fromP[j]+k-1]
         , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP), type="l", ylab="")
    par(new=1)
   }
   plot(theta,dP[,toP[j]], main=title, sub="Generalized Partial Credit Model"
        , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP)
        , type="l", ylab="dICRF")
   par(new=0)
  }
  if( plot >= 2 ){
   for( j in 1:(nitems) ){
    title=paste("Plot of the Derivative of IRF"
                , " (with natural category weight) of "
                , iname[j]," (ncat=",ncat[j],")")
    plot(theta,dPt[,j], main=title, sub="Generalized Partial Credit Model"
         , xlim=c(min(theta),max(theta)), ylim=, type="l", ylab="dIRF")
    par(new=0)
   }
  }
 }


 return( list(dP=dP, dPt=dPt
              , fromP=fromP, toP=toP) )

} # end of dicrfP


#' Calculation of the Derivative of Item Response Function of Nominal
#'  Items
#'
#' @param param Item Parameter Data Frame
#' @param theta Discrete theta values
#' @param smallP Minimum value of probability
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#' @param debug = 1 to print intemediate result
#' @param DinP = 1 to include D=1.7 in logistic function
#'
#'
#'
#'
#' @export
#'
dicrfN <- function( param, theta, smallP=0, print=0, plot=0, debug=0, DinP=1 ){
 # derivative of Generalized Partial Credit ICRF in Nominal format
 # Shin-ichi Mayekawa
 # 121123
 #
 #
 #  theta    npoints x 1 theta values
 #
 # Values: as list
 #   P       length(theta) x sum(ncat) icrf matrix
 #   PP      length(theta) x sum(ncat) icbf  matrix
 #   Pt      length(theta) x nitems item response function
 #
 #

 # argument name
 paramname=as.character(substitute(param))

 isdf=0
 if( is.data.frame(param) ){
  # param and weight given as data frame
  isdf=1
  param=checkparam( param, "N", "dicrfN" )
  if( is.null(param) ){
   cat("\n\nInput Item Parameter ERROR.\n\n")
   return()
  }
 }


 # argument name
 paramname=as.character(substitute(param))
 isdf=0
 if( is.data.frame(param) ){
  # param and weight given as data frame
  isdf=1
  param=checkparam( param, "N", "icrfN" )
  if( is.null(param) ){
   cat("\n\nInput Item Parameter ERROR.\n\n")
   return()
  }
 }

 # const
 if( is.matrix(theta) ) theta=as.vector(theta)
 npoints=length(theta)
 thetac=format(theta,digits=3)

 nitems=nrow(param)
 iname=rownames(param)
 ncat=param[,1]
 ncat1=ncat-1
 ncp=ncol(param)
 a=matrix(NA,nitems,max(ncat1))
 c=matrix(NA,nitems,max(ncat1))
 for( j in 1:nitems ){
  a[j,1:ncat1[j]]=param[j,2:(2+ncat1[j]-1)]
  c[j,1:ncat1[j]]=param[j,(ncat[j]+1):(2*ncat1[j]+1)]
 }

 # icrf
 temp=icrfN( param, theta, smallP=smallP )
 P=temp$P
 fromP=temp$fromP
 toP=temp$toP
 Pname=colnames(P)
 if( debug ) Print(P,fromP,toP)


 # sum_{k=0}^m-j  k P_{kj}
 vecv=matrix( as.vector( unlist( mapply(seq,0,ncat1) ) ),,1 )
 sumP=matrix(0,npoints,nitems)
 for( j in 1:nitems ){
  sumP[,j]=P[,fromP[j]:toP[j]]%*%matrix(c(0,a[j,1:ncat1[j]]),,1)
 }
 if( debug ) Print(P,vecv,sumP)
 # d
 dP=P; dPt=matrix(0,npoints,nitems)
 for( j in 1:nitems ){
  Pj=P[,fromP[j]:toP[j]]
  dP[,fromP[j]:toP[j]]=Pj%*%diag(c(0,a[j,1:ncat1[j]]))-Pj*sumP[,j]
  # dirf
  dPt[,j]=dP[,fromP[j]:toP[j]]%*%matrix((0:ncat1[j]),,1)
 }


 if( print > 0 ){
  cat("\nDerivative of ICRF of Generalized Partial Credit Items \n ")
  if( isdf )
   cat(" item parameter data frame name =", paramname,"\n")
  else
   cat(" item parameter matrix name =", paramname,"\n")
  cat("  # of item parameters =",nitems,
      ",  # of theta points =", npoints,"\n")
  Print(param, digits=3)
  cat("\nDerivative of Item Category Response Functions (ICRF)\n")
  Print(dP, digits=2)
  cat("\nDerivative of Item Response Functions (IRF)"
      ," with natural category weights\n")
  Print(dPt, digits=2)
 }

 if( plot ){
  for( j in 1:(nitems) ){
   title=paste("Plot of the Derivative of ICRF of "
               , iname[j]," (ncat=",ncat[j],")")
   minP=min(dP); maxP=max(dP)
   for( k in 1:ncat1[j] ){
    plot(theta,dP[,fromP[j]+k-1]
         , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP), type="l", ylab="")
    par(new=1)
   }
   plot(theta,dP[,toP[j]], main=title, sub="Generalized Partial Credit Model"
        , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP)
        , type="l", ylab="dICRF")
   par(new=0)
  }
  if( plot >= 2 ){
   for( j in 1:(nitems) ){
    title=paste("Plot of the Derivative of IRF"
                , " (with natural category weight) of "
                , iname[j]," (ncat=",ncat[j],")")
    plot(theta,dPt[,j], main=title, sub="Generalized Partial Credit Model"
         , xlim=c(min(theta),max(theta)), ylim=, type="l", ylab="dIRF")
    par(new=0)
   }
  }
 }


 return( list(dP=dP, dPt=dPt
              , fromP=fromP, toP=toP) )

} # end of dicrfN





#' Calculation of the Derivative of Partial Credit ICRF in Nominal format
#'
#' @param param Item Parameter Data Frame
#' @param theta Discrete theta values
#' @param smallP Minimum value of probability
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#' @param debug = 1 to print intemediate result
#'
#' @return A list of \cr
#'   dP             length(theta) x sum(ncat) dICRF matrix \cr
#'   dPt            length(theta) x nitems derivative of
#'                   item response function \cr
#'   fromP, toP    indexes of each ICRF in P \cr
#'   fromPP, toPP  indexes of each ICBRF in PP
#'
#' @export
#'
dicrfPN0 <- function( param, theta, smallP=0, print=0, plot=0, debug=0 ){
 # derivative of Generalized Partial Credit ICRF in Nominal format
 # Shin-ichi Mayekawa
 # 120215
 # checkparam: 120224
 # renamed as dicrfPN0: 120301
 # smallP: 121022
 #
 # Args:
 #  param    nitems x max(ncat) matrix of item parameters or data frame
 #    An Example of Item Parameter datastet:  maxncat = 5
 #    Note that the subscript of the c-parameter ranges between
 #         1 and ncat-1 (= maxscore)  (We set c0=0.)
 #
 #          ncat     A            c1     c2     c3     c4
 #
 #    Item1  2     1.0           0.0     .      .      .
 #    Item2  2     1.0           0.0     .      .      .
 #    Item3  3     1.5          -1.0    0.0     .      .
 #    Item4  4     1.5          -1.0    0.0    1.0     .
 #    Item5  5     1.5          -1.0    0.0    1.0    2.0
 #
 #    ICRF of category k of item j at theta, namely, p_{jk}(theta)
 #    is defined as
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0}^{ncat[j]-1} Ez_{jk}(theta)
 #      where
 #     Ez_{jk}(theta) = exp( z_{jk}(theta) )
 #       and
 #     z_{jk}(theta) = a_j k theta + c_{k}, k=0,1, ..., ncat[j]-1
 #    with c_{j0}=0.
 #
 #
 #  theta    npoints x 1 theta values
 #
 # Values: as list
 #   P       length(theta) x sum(ncat) icrf matrix
 #   PP      length(theta) x sum(ncat) icbf  matrix
 #   Pt      length(theta) x nitems item response function
 #
 #

 # argument name
 paramname=as.character(substitute(param))

 isdf=0
 if( is.data.frame(param) ){
  # param and weight given as data frame
  isdf=1
  param=checkparam( param, "PN", "dicrfPN0" )
  if( is.null(param) ){
   cat("\n\nInput Item Parameter ERROR.\n\n")
   return()
  }
 }


 # const
 if( is.matrix(theta) ) theta=as.vector(theta)
 npoints=length(theta)
 thetac=format(theta,digits=3)

 nitems=nrow(param)
 iname=rownames(param)
 ncat=param[,1]
 ncat1=ncat-1
 a=param[,2]
 c=param[,3:(1+max(ncat)),drop=0]

 # icrf
 temp=icrfPN0( param, theta, smallP=smallP )
 P=temp$P
 fromP=temp$fromP
 toP=temp$toP
 Pname=colnames(P)
 if( debug ) Print(P,fromP,toP)


 # sum_{k=0}^m-j  k P_{kj}
 vecv=matrix( as.vector( unlist( mapply(seq,0,ncat1) ) ),,1 )
 sumP=matrix(0,npoints,nitems)
 for( j in 1:nitems ){
  sumP[,j]=P[,fromP[j]:toP[j]]%*%vecv[fromP[j]:toP[j],1]
 }
 if( debug ) Print(P,vecv,sumP)
 # d
 dP=P; dPt=matrix(0,npoints,nitems)
 for( j in 1:nitems ){
  dP[,fromP[j]:toP[j]]=a[j]*P[,fromP[j]:toP[j]] *
          ( matrix(1,npoints)%*%t(vecv[fromP[j]:toP[j],1]) - sumP[,j] )
  # dicrf
  dPt[,j]=dP[,fromP[j]:toP[j]]%*%vecv[fromP[j]:toP[j],1]
 }


 if( print > 0 ){
  cat("\nDerivative of ICRF of Generalized Partial Credit Items \n ")
  if( isdf )
   cat(" item parameter data frame name =", paramname,"\n")
  else
   cat(" item parameter matrix name =", paramname,"\n")
  cat("  # of item parameters =",nitems,
      ",  # of theta points =", npoints,"\n")
  Print(param, digits=3)
  cat("\nDerivative of Item Category Response Functions (ICRF)\n")
  Print(dP, digits=2)
  cat("\nDerivative of Item Response Functions (IRF)"
     ," with natural category weights\n")
  Print(dPt, digits=2)
 }

 if( plot ){
  for( j in 1:(nitems) ){
   title=paste("Plot of the Derivative of ICRF of "
             , iname[j]," (ncat=",ncat[j],")")
   minP=min(dP); maxP=max(dP)
   for( k in 1:ncat1[j] ){
    plot(theta,dP[,fromP[j]+k-1]
         , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP), type="l", ylab="")
    par(new=1)
   }
   plot(theta,dP[,toP[j]], main=title, sub="Generalized Partial Credit Model"
        , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP)
        , type="l", ylab="dICRF")
   par(new=0)
  }
  if( plot >= 2 ){
   for( j in 1:(nitems) ){
    title=paste("Plot of the Derivative of IRF"
               , " (with natural category weight) of "
                , iname[j]," (ncat=",ncat[j],")")
    plot(theta,dPt[,j], main=title, sub="Generalized Partial Credit Model"
         , xlim=c(min(theta),max(theta)), ylim=, type="l", ylab="dIRF")
    par(new=0)
   }
  }
 }


 return( list(dP=dP, dPt=dPt
              , fromP=fromP, toP=toP) )

} # end of dicrfPN0








#' Calculation of the Derivatives of the Item Category Response Function
#' with respect to Item Parameters
#'
#' @param paramj item parameters data frame for ONE item
#' @param theta Discrete theta values
#' @param weight Weight data frame:  NOT used.
#' @param smallP Minimum value of probability
#' @param DinP = 1 to include D=1.7 in logistic function
#' @param thmin Minimum value of discrete thata value
#' @param thmax Maximum value of discrete thata value
#' @param npoints # of discrete points for theta
#' @param Pj icrf:   npoints x (ncatj-1)  (no zero category)  or NULL
#' @param PPj icbrf of the Graded Response Model   or NULL
#' @param zero = 1 to include the zero-th category in output
#' @param cat.first  = 1 to chage the category fist in the rows of Jack.
#' @param log = 1 to obtain the log Jacobian: d log(ICRF) d param
#' @param print = 1 to print result
#'
#' @return
#'   list of (Jack, Pj, PPj) \cr
#'   Jack     (length(theta) x (ncatj-1)) x ncatj \cr
#'             derivative of vec(Pj) with respect to (a, b1, b2, ...) \cr
#'           If cat.first = 0 \cr
#'             theta changes first, then k changes from 1 to ncatj \cr
#'           If cat.first = 1, category(k) changes first, then theta. \cr
#'   Pj Item Category Response Function (icrf)  (length(theta) x (ncatj-1))\cr
#'   PPj will be output when type="G" or "Gn".
#'
#'
#' @examples
#' res=dirf_p( paramS1[3,], npoints=5, print=1 )
#' res=dirf_p( paramS1[3,], npoints=5, print=1, zero=1 )
#' res=dirf_p( paramS1[3,], npoints=5, print=1, cat.first=1 )
#'
#'
#' @export
#'

dirf_p <- function( paramj, theta=NULL, weight=NULL, smallP=0, DinP=1
                  , thmin=-4, thmax=4, npoints=21
                  , Pj=NULL, PPj=NULL
                  , zero=0, cat.first=0, log=0, print=0 ){
 # derivative of ICRF with respect to item paramters
 # Shin-ichi Mayekawa
 # 121013
 # zero=1: 121018
 # cat.first: 121022
 # log: 121022,24
 # smallP: 121022
 # calculation of PPj from Pj for Graded: 121105
 # type = "P": 121106,08(London)
 # DinP: 121109(London)
 # Check Pj in denominator by smallP: 121115
 # Check Pj in denominator by 1e-7: 121119
 # P with ncat=2: 121122
 # nominal response model: 121123
 # bugfix for type "B": 20150614
 # error check: 201506014
 # rownames for Jack: 20150615
 # check if c is na or null: 20150615
 # Bn and Gn: 20180211
 # roxgen2 comment: 20180212
 #
 #
 #
 # Args:
 #  paramj   item parameters data frame for one item
 #  theta    lenth x 1 theta values
 #  zero     = 1 to include the 0-th category
 #  Pj       icrf:   npoints x (ncatj-1)  (no zero category)
 #  PPj      icbrf for the Graded Response Model
 #  cat.first = 1 to chage the category fist in the rows of Jack.
 #           ( 10000 theta points or less.)
 #  log      = 1 to obtain the log Jacobian: d log(ICRF) d param
 #  smallp   = small value for icrf etc
 #
 # Values:
 #   list of (Jack, Pj, PPj)
 #   Jack     (length(theta) x (ncatj-1)) x ncatj
 #             derivative of vec(Pj) with respect to (a, b1, b2, ...)
 #           If cat.first = 0
 #             theta changes first, then k changes from 1 to ncatj
 #           If cat.first = 1, category(k) changes first, then theta.
 #
 # Needs:
 #   icrfB, icrfG, icrfPN
 #


 JacobianMat <- function( parvec, infcn, eps = 1e-06 ){
  # Function to calculate the difference-quotient approx gradient matrix
  # of an arbitrary input function infcn.
  # Now recoded to use central differences !
  #
  # Original version was Gradmat in Sec4NotF09.pdf
  #

  dd = length(parvec)
  aa = length(infcn(parvec))
  epsmat = (diag(dd) * eps)/2
  gmat = array(0, dim = c(aa, dd))
  for(i in 1:dd)
   gmat[, i]=( infcn(parvec + epsmat[, i] )
               - infcn(parvec - epsmat[, i] ) )/eps

  # return a vector if infcn is unidimensional
  if( aa == 1 ) gmat=c(gmat)
  return( gmat )

 } # end of JacobianMat


 icrfP00 <- function( paramnum ){
  # This is to be used in conjunction with JacobianMat
  paramj0=paramj
  paramj0[,4:(4+length(paramnum)-1)]=paramnum
  res=irf( paramj0, theta, zero=0, print=0 )$ICRF
  res=matrix(res,,1)
  # vertically stack all categories
  return(  matrix( res, , 1 )  )
 } # end of icrfP00


 # error check
 if( nrow(paramj) > 1 ){
  cat("\n\nerror1:(dirf_p)  paramj data frame contains more than 2 items.\n ")
  return()
 }

 # const
 ncatj=paramj$ncat
 ncatj1=ncatj-1
 typej=paramj$type

 # # of parameters to be estimated
 if( typej %in% c("B3","Bn3") ) np=ncatj+1
 else if( typej == "N" ) np=2*ncatj1
 else np=ncatj

 if( is.null(theta) ){
  theta=seq(thmin,thmax,length.out=npoints)
 }
 else{
  if( is.matrix(theta) ) theta=as.vector(theta)
  thmin=theta[1]; thmax=theta[length(theta)]
 }
 lenth=length(theta)

 if( DinP == 1 ) DD=1.7
 else DD=1

 # item param for type Bx, G, P, PN
 locparam=grep("^p[[:digit:]]$",colnames(paramj))
 if( length(locparam) >= 1 ) locparam=locparam[1]
 else locparam=grep("^a$",colnames(paramj))
 a=paramj$p1
 b=as.matrix(paramj[1,(locparam+1):(locparam+ncatj1)],1)
 if( typej %in% c("G","Gn") ) b=cbind(b,9999)
 c=paramj$p3
 # 20150615
 if( is.null(c) ) c=0; if( is.na(c) ) c=0

 # storage
 Jack=matrix(0,ncatj1*lenth,np)
 colnames(Jack)=paste("p",1:np,sep="")
 # thname=paste( "t", 1:lenth, sep="" )
 thname=format(theta,digits=3)
 catname=paste( "c", 1:ncatj1, sep="" )
 rnJ=apply( expand.grid( thname, catname ), 1, paste, collapse="-" )
 rownames(Jack)=rnJ

  # icrf at current paramnum
 if( length(grep("^B[[:digit:]]*$", typej)) > 0 ){

  if( is.null(Pj) ){
   # icrf for k=1
   Pj=icrfB( paramj, theta, smallP=smallP, zero=0 )
  }
  # Jacobian
  Jack[,1]=1.7/(1-c)*(theta-b[1])*(Pj-c)*(1-Pj)
  Jack[,2]=-1.7/(1-c)*a*(Pj-c)*(1-Pj)
  if( typej == "B3" ) Jack[,3]=(1-Pj)/(1-c)

 } # end of type = "B"
 else if( length(grep("^Bn[[:digit:]]*$", typej)) > 0 ){

  dd=dnorm(a*(theta-b[1]))
  if( is.null(Pj) ){
   # icrf for k=1
   Pj=icrfBn( paramj, theta, smallP=smallP, zero=0 )
  }
  # Jacobian
  Jack[,1]=dd*(theta-b[1])*(1-c)
  Jack[,2]=-dd*a*(1-c)
  if( typej == "Bn3" ) Jack[,3]=(1-Pj)/(1-c)

 } # end of type = "Bn"
 else if( typej == "G" ){

  if( is.null(Pj) ){
   # icrf : no zero-th category
   temp=icrfG( paramj, theta, smallP=smallP )
   Pj=temp$P[,2:ncatj,drop=0];
   PPj=temp$PP[,2:ncatj,drop=0]; PPj=cbind(PPj,0)
  }
  else if( is.null(PPj) ){
   # calculation of icbrf from icrf
   PPj=matrix( 0,lenth, ncatj)
   for( k in ncatj1:1 ){
    PPj[,k]=Pj[,k] + PPj[,k+1]
   }
  }
  # Jacobian
  for( k in 1:ncatj1 ){
   temp=1.7*(
    (theta-b[k])*PPj[,k]*(1-PPj[,k])-(theta-b[k+1])*PPj[,k+1]*(1-PPj[,k+1]) )
   Jack[((k-1)*lenth+1):(k*lenth),1]=temp
   for( z in 1:ncatj1 ){
    temp=0
    if( k == z ) temp=-1.7*a*PPj[,z]*(1-PPj[,z])
    else if( (k+1) == z ) temp=1.7*a*PPj[,z]*(1-PPj[,z])
    Jack[((k-1)*lenth+1):(k*lenth),1+z]=temp
   } # z
  } # k

 } # end of type = "G"
 else if( typej == "Gn" ){

  maxZ=700
  Z=-(a*outer(c(b),theta,"-"))
  Z[which(Z > maxZ)]=maxZ
  Z[which(Z < -maxZ)]=-maxZ
  dPj=t( dnorm(Z) )

  if( is.null(Pj) ){
   # icrf : no zero-th category
   temp=icrfGn( paramj, theta, smallP=smallP )
   Pj=temp$P[,2:ncatj,drop=0];
   PPj=temp$PP[,2:ncatj,drop=0]; PPj=cbind(PPj,0)
  }
  else if( is.null(PPj) ){
   # calculation of icbrf from icrf
   PPj=matrix( 0,lenth, ncatj)
   for( k in ncatj1:1 ){
    PPj[,k]=Pj[,k] + PPj[,k+1]
   }
  }
  # Jacobian
  for( k in 1:ncatj1 ){
   temp=( (theta-b[k])*dPj[,k]-(theta-b[k+1])*dPj[,k+1] )
   Jack[((k-1)*lenth+1):(k*lenth),1]=temp
   for( z in 1:ncatj1 ){
    temp=0
    if( k == z ) temp=-a*dPj[,z]
    else if( (k+1) == z ) temp=a*dPj[,z]
    Jack[((k-1)*lenth+1):(k*lenth),1+z]=temp
   } # z
  } # k

 } # end of type = "Gn"
 else if( typej == "PN" ){

  if( is.null(Pj) ){
   # icrf : no zero-th category
   Pj=icrfPN( paramj, theta, smallP=smallP )$P[,2:ncatj,drop=0]
  }
  # Jacobian
  thetab=matrix(theta,lenth,ncatj1)-matrix(1,lenth,1)%*%b
  Pqthetab=rowSums( (Pj%*%diag(1:ncatj1))*thetab )
  for( k in 1:ncatj1 ){
   temp=1.7*Pj[,k]*(k*thetab[,k]-Pqthetab)
   Jack[((k-1)*lenth+1):(k*lenth),1]=temp
   for( z in 1:ncatj1 ){
    if( k == z ) temp=-1.7*a*z*Pj[,z]*(1-Pj[,z])
    else  temp=1.7*a*z*Pj[,k]*Pj[,z]
    Jack[((k-1)*lenth+1):(k*lenth),1+z]=temp
   } # z
  } # k

 } # end of type = "PN"
 else if( typej == "P" ){

  if( is.null(Pj) ){
   # icrf : no zero-th category
   Pj=icrfP( paramj, theta, smallP=smallP, DinP=DinP )$P[,2:ncatj,drop=0]
  }
  # Jacobian
  thetab=matrix(theta,lenth,ncatj1)-matrix(1,lenth,1)%*%b
  thetabc=matrix( t( apply(thetab,1,cumsum) ), lenth )
  Pqthetabc=rowSums( Pj*thetabc )
  # reverse cum sum
  Pqac=DD*a*matrix( t(apply(Pj[,ncatj-(1:ncatj1),drop=0],1,cumsum)) ,lenth)
  Pqac=Pqac[,ncatj-(1:ncatj1),drop=0] # reverse again
  for( k in 1:ncatj1 ){
   temp=DD*Pj[,k]*(thetabc[,k]-Pqthetabc)
   Jack[((k-1)*lenth+1):(k*lenth),1]=temp
   for( z in 1:ncatj1 ){
    if( k >= z ) temp=-DD*a + Pqac[,z]
    else  temp=Pqac[,z]
    Jack[((k-1)*lenth+1):(k*lenth),1+z]=Pj[,k]*temp
   } # z
  } # k
  #  Jack=JacobianMat( c(a,b), icrfP00, eps = 1e-06 )
  #  Print(Jack)

 } # end of type = "P"
 else if( typej == "N" ){

  if( is.null(Pj) ){
   # icrf : no zero-th category
   Pj=icrfN( paramj, theta, smallP=smallP, DinP=DinP )$P[,2:ncatj,drop=0]
  }

  # Jacobian
  for( k in 1:ncatj1 ){
   for( z in 1:ncatj1 ){
    if( k == z ){
     Jack[((k-1)*lenth+1):(k*lenth),k]=theta*Pj[,k]*(1-Pj[,k])
     Jack[((k-1)*lenth+1):(k*lenth),ncatj1+k]=Pj[,k]*(1-Pj[,k])
    }
    else{
     Jack[((k-1)*lenth+1):(k*lenth),z]=-theta*Pj[,k]*Pj[,z]
     Jack[((k-1)*lenth+1):(k*lenth),ncatj1+z]=-Pj[,k]*Pj[,z]
    }
   }
  }

 } # end of type = "N"



 if( zero == 1 ){
  # dP0dz = d(1-sum_k Pj)dz = 0 - sum_k dPjdz
  temp=0
  rntemp=paste(thname,"-c0",sep="")
  for( k in 1:ncatj1 ){
   temp=temp+Jack[((k-1)*lenth+1):(k*lenth),]
  } # k
  Jack=rbind(-temp, Jack)
  rownames(Jack)=c(rntemp,rnJ)
  Pj=cbind(1-apply(Pj,1,sum),Pj)
  colnames(Pj)[1]="0"
 }

 if( log == 1 ){
  # Pj[which(Pj < smallP)]=smallP
  # Pj[which(Pj < 1e-307)]=1e-307
  # Pj[which(Pj < 1e-30)]=1e-30
  # This 1e-7 seems to be the best.
  Pj[which(Pj < 1e-7)]=1e-7
  Jack=Jack/as.vector(matrix(Pj,,1))
 }

 if( cat.first ){
  # sort the rows of Jack to make the category change first.
  ncat1=ncatj-1; if( zero == 1 ) ncat1=ncatj
  index=10000*rep(1:lenth,ncat1) + as.vector(
                  sapply(1:ncat1, function(x) rep(x,lenth)) )
  Jack=Jack[order(index),,drop=0]
 }

 if( print > 0 ){
  Print(lenth, zero, log, cat.first, smallP,paramj)
  Print(Pj,Jack, fmt="5.2")
 }

 return( list(Jack=Jack, Pj=Pj, PPj=PPj) )

} # end of dirf_p



