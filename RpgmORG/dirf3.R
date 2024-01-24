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
#' @param print = 1 to print result
#' @param debug = 1 to print intemediate result
#' @param plot = 1 to plot result
#'
#' @return A list of \cr
#'   dICRF, dIRF, dTRF, fromP, toP=toP, vecv, minscore_i, maxscore_i,
#'   minscore_t, maxscore_tt
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
#'
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

dirf <- function( param, theta=NULL, weight=NULL, zero=1, smallP=0
                  , thmin=-4, thmax=4, npoints=21, DinP=1
                  , numderiv=0, eps=1e-6
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
 #   icrfB, icrfG,  icrfPN
 #   dicrfB, dicrfG,  dicrfPN
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
 locG=which(type=="G")
 locP=which(type=="P")
 locPN=which(type=="PN")
 locN=which(type=="N")
 nB=length(locB)
 nN=length(locN)
 nP=length(locP)
 nPN=length(locPN)
 nG=length(locG)

 # locations in ircf-storage in the original order
 dicrfP=matrix(0,npoints,sum(ncat))
 cnP=character(max(ncat))
 toP=cumsum(ncat)
 fromP=toP-ncat+1
 ll=mapply(seq,fromP,toP)   # simplify removed: 120928
 if( is.matrix(ll) ){
  llB=as.vector(ll[,locB])
  llG=as.vector(ll[,locG])
  llN=as.vector(ll[,locN])
  llP=as.vector(ll[,locP])
  llPN=as.vector(ll[,locPN])
 }
 else if( is.list(ll) ){
  llB=unlist(ll[locB])
  llG=unlist(ll[locG])
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
   if( debug > 0 ) Print(pB,diccB, digits=2)
   dicrfP[,llB]=diccB;
   cnP[llB]=colnames(diccB)
  }
  if( nPN > 0 ){
   pP=as.matrix( param[locPN,3:ncp,drop=0] )
   rownames(pP)=iname[locPN]
   res=dicrfPN( pP, theta, smallP=smallP )
   diccP=res$dP
   if( debug > 0 ) Print(pP,diccP, digits=2)
   dicrfP[,llPN]=diccP;
   cnP[llPN]=colnames(diccP)
  }
  if( nP > 0 ){
   pP=as.matrix( param[locP,3:ncp,drop=0] )
   rownames(pP)=iname[locP]
   res=dicrfP( pP, theta, DinP=DinP, smallP=smallP )
   diccP=res$dP
   if( debug > 0 ) Print(pP,diccP, digits=2)
   dicrfP[,llP]=diccP;
   cnP[llP]=colnames(diccP)
  }
  if( nG > 0 ){
   pG=as.matrix( param[locG,3:ncp,drop=0] )
   rownames(pG)=iname[locG]
   res=dicrfG( pG, theta, smallP=smallP )
   diccG=res$dP
   if( debug > 0 ) Print(pG,diccG, digits=2)
   dicrfP[,llG]=diccG;
   cnP[llG]=colnames(diccG)
  }
  if( nN > 0 ){
   pN=as.matrix( param[locN,3:ncp,drop=0] )
   rownames(pN)=iname[locN]
   res=dicrfN( pN, theta, smallP=smallP )
   diccN=res$dP
   if( debug > 0 ) Print(pN,,diccN, digits=2)
   dicrfP[,llN]=diccN;
   cnP[llN]=colnames(diccN)
  }
  colnames(dicrfP)=cnP; rownames(dicrfP)=format(theta,digits=3)
 } # end of analytic derivation
 else{
  # numeric derivatives
  dicrfP=dicrf_num( param=param, theta=theta, eps=eps )
  Print(dicrfP)
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
   main=paste("Plot of the Derivative of ICRF of "
              ,iname[j]," ( type = ",type[j],", ncat=",ncat[j]," )")
   sub=paste("param = ", paste(format(param[j,4:(ncat[j]+3)],digits=3)
                               , collapse=",  "))
   if( length(grep("^B[[:digit:]]*$", param$type[j])) > 0 )
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
    title=paste("Plot of the Derivative of IRF of "
                , iname[j]," ( type = ",type[j],", ncat=",ncat[j]
                , ", score range = [",minscore_i[j],",", maxscore_i[j],"] )")
    sub=paste("param = ", paste(format(param[j,4:(ncat[j]+3)],digits=3)
                                , collapse=",  "), "  (with weights)")
    if( length(grep("^B[[:digit:]]*$", param$type[j])) > 0 )
     sub=paste("param = ", paste(format(param[j,4:(ncat[j]+4)],digits=3)
                                 , collapse=",  "), "  (with weights)")
    minP=min(dirf); maxP=max(dirf)
    plot(theta,dirf[,j], main=title, sub=sub
         , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP)
         , type="l", ylab="dIRF")
   }
   title=paste("Plot of the Derivative of TRF  (# of items =", nitems
               , ", score range = [",minscore_t,",", maxscore_t,"] )")
   minP=min(dtrf); maxP=max(dtrf)
   plot(theta,dtrf, main=title
        , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP)
        , type="l", ylab="dTRF")
  }
 }

 return( list(dICRF=dicrfP, dIRF=dirf, dTRF=dtrf, fromP=fromP, toP=toP
              , vecv=vecv, minscore_i=minscore_i, maxscore_i=maxscore_i
              , minscore_t=minscore_t, maxscore_t=maxscore_t) )


} # end of dirf





res1=dirf( paramS1, plot=0 )$dICRF
res2=dirf( paramS1, numderiv=1, plot=0 )$dICRF
Print(max(abs(res1-res2)))

