#' Calculation of Item Response Function
#'
#' This function calculates icrf, irf and trf. \cr
#' Japanese help file: (\link[lazy.irt]{irf_JPH})
#'
#' @param param Item Parameter Data Frame
#' @param theta Discrete theta values
#' @param weight Weight data frame
#' @param zero = 0 to exclude the zero-th category from output
#' @param smallP Minimum value of probability
#' @param thmin Minimum value of discrete thata valuepackage
#' @param thmax Maximum value of discrete thata value
#' @param npoints # of discrete points for theta
#' @param DinP = 1 to include D=1.7 in logistic function
#' @param print = 1 to print result
#' @param debug = 1 to print intemediate result
#' @param plot = 1 to plot result
#' @param colors = color of the lines in plot
#' @param plotchar = pch vector to be used in the plot
#' @param linetype = lty vector to be used in the plot
#'
#' @details
#' It is assumed that the data frame has the following columns
#'   in the following order.\cr
#'\cr
#'   name, type, ncat, p1, p2, ..., p[ncat+1]\cr
#'\cr
#'\cr
#'   When item_type == "B"  or "B2"  or  "B3"   or   "Bn", ( ncat == 2 )\cr
#'     p1 is the discrimination parameter\cr
#'     p2 is the difficulty parameter\cr
#'     p3 is the gussing parameter\cr
#'\cr\cr
#'   When item_type == "G" or "Gn",\cr
#'     p1 is the discrimination parameter\cr
#'     p2 is the category threshold parameter b_{j1}\cr
#'     p3 is the category threshold parameter b_{j2}\cr
#'\cr
#'     p_ncat is the category threshold parameter b_{j,[ncat-1]}\cr
#'\cr\cr
#'   When item_type == "P",\cr
#'     p1 is the discrimination parameter\cr
#'     p2 is the step parameter b_{j1}\cr
#'     p3 is the step parameter b_{j2}\cr
#'\cr
#'     p_ncat is the step parameter b_{j,[ncat-1]}\cr
#'\cr\cr
#'   When item_type == "PN,\cr
#'     p1 is the slope parameter\cr
#'     p2 is the intercept parameter b_{j1}\cr
#'     p3 is the intercept parameter b_{j2}\cr
#'\cr
#'     p_ncat is the intercept parameter b_{j,[ncat-1]}\cr
#'\cr\cr
#'   When item_type == "N,\cr
#'     p1 - p_(ncat[j]-1) are the slope parameters\cr
#'     p_ncat[j] - p_2*(ncat[j]-1) are the intercept parameters\cr
#'\cr
#'   Regardless of the item_type, there will be ncat[j]\cr
#'   item paramters for item j, except for the Binary Items\cr
#'   which have the gusseing param as the ncat[j]+1st,\cr
#'   and the Nominal Items which have 2*(ncat[j]-1) item parameters.\cr
#'
#' @return
#'    list( ICRF, IRF, TRF, fromP, toP )\cr
#' \cr
#'    where\cr
#'      ICRF  npoints x sum(ncat)\cr
#'      IRF   npoints x nitems   weighted by item category weight\cr
#'      TRF   npoints x 1        weighted by item category weight\cr
#'                                           and item weight\cr
#'      fromP, toP    location of each item category in ICRF\cr
#' \cr
#'  Needs:\cr
#'    icrfB, icrfBN, icrfG, idrfGn, icrfPN, icrfP, icrfN, checkparam\cr
#' \cr
#'
#' @examples
#' Cards="
#' name type ncat   p1   p2   p3   p4
#'   Q1    B    2  0.9  0.0,,
#'   Q2    B3   2  1.0  0.0  0.2,
#'   Q3    G    4  1.0 -2.0  0.0  2.0
#'   Q4    P    4  0.8 -2.0  0.0  2.0
#' "; param <- cards( Cards, header=1 )
#' irf( param, plot=1 )
#'
#' @export
#'

irf <- function( param, theta=NULL, weight=NULL, zero=1, smallP=1e-9
               , thmin=-4, thmax=4, npoints=21, DinP=1
               , print=1, debug=0
               , plot=0, colors="black", plotchar=NULL, linetype=NULL ){
 # calculation of ICRF, IRF, and TRF
 # Shin-ichi Maykeawa
 # 120208,09,10,11
 # case when mapply returns matrix ll: 120212
 # 120213,15
 # negative/non integer weight: 120221
 # iname as factor: 120223
 # exclude zero: 120223
 # clean input data frame with NA: 120224
 # checkparam: 120224,29
 # theta=NULL: 120304
 # print param value: 120305
 # 120306
 # type Bxx items: 120919,21
 # simplify removed: 120928
 # type P items: 121016
 # smallP: 121022
 # DinP: 121109(London)
 # plot legend: 121119
 # nominal response model: 121123
 # bug fix for colors=: 20150315
 # bug fix for plot and zero=0: 20150321 IC
 # new plotchar and linetype: 20161123
 # bugfix for nitems=1 and zero=0: 20171130
 # smallP=1e-9: 20180129
 # digits -> fmt: 20180208
 # Bn: 20180208
 # GN: 20180209
 # remove excessive cat weight: 20230718
 # JPH: 20231028cot
 #

 # Args:
 #
 #   param   data.frame containing item parameters as
 #     name type ncat  a  p1  p2  p3   .....   p[ncat-1]
 #    See the description in read.param.R
 #
 #   theta   vector of length npoints of theta points
 #
 #   weight  data.frame containing item and category weight
 #     name type ncat  w  v0 V1  v2  v3   .....   p[ncat-1]
 #    See the description in read.weight.R
 #
 #    Weight will be used to calculate item response functions and
 #    test response function.
 #    If weight == NULL, a set of natural category wieight will be used
 #
 #   zero  = 0 to exclude the probability of the zero-th category
 #
 #
 # Value
 #
 #   list( ICRF, IRF, TRF, fromP, toP )
 #
 #   where
 #     ICRF  npoints x sum(ncat)
 #     IRF   npoints x nitems   weighted by item category weight
 #     TRF   npoints x 1        weighted by item category weight
 #                                          and item weight
 #     fromP, toP    location of each item category in ICRF
 #
 # Needs:
 #   icrfB, icrfG, icrfPN, icrfP, icrfN, checkparam
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
 #Print(param,nitems,ncp,iname,ncat,type)

 # clearn up
 loca=which(colnames(param) == "p1")
 if( length(loca) == 0 ) loca=which(colnames(param) == "a")
 for( i in 1:nitems ){
  if( ( type[i] == "Bn" || length(grep("^B[[:digit:]]*$", type[i])) > 0 )
      & ncat[i]+loca+1 <= ncp )
   param[i,(ncat[i]+loca+1):ncp]=NA
  else if( (type[i] == "G"  |  type[i] == "PN") & ncat[i]+loca <= ncp )
   param[i,(ncat[i]+loca):ncp]=NA
 }
 locna=is.na(param)

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
  # remove excessive category weights: 20230718
  for( i in 1:nitems ){
   if( ncat[i] < ncol(v)) v[i,(ncat[i]+1):ncol(v)]=NA
  }
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
 ICRFP=matrix(0,npoints,sum(ncat))
 cnP=character(max(ncat))
 toP=cumsum(ncat)
 fromP=toP-ncat+1
 ll=mapply(seq,fromP,toP)  # simplify removed: 120928
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

 if( debug > 0 ) Print(locB, locBn,locG, locGn,locP,locN
                       ,ll,llB, llBn,llG,llGn,llP,llN)

 # for each item type, calculate icrf and store them in ICRFP
 if( nB > 0 ){
  # 3PLM
  pB=as.matrix( param[locB,4:6,drop=0] )
  rownames(pB)=iname[locB]
  colnames(pB)=c("a","b","c")
  iccB=icrfB( pB, theta, zero=1, smallP=smallP )
  if( debug > 0 ) Print(pB,iccB,fmt="6.3")
  ICRFP[,llB]=iccB;
  cnP[llB]=colnames(iccB)
 }
 if( nBn > 0 ){
  # 3PLM
  pBn=as.matrix( param[locBn,4:6,drop=0] )
  rownames(pBn)=iname[locBn]
  colnames(pBn)=c("a","b","c")
  iccBn=icrfBn( pBn, theta, zero=1, smallP=smallP )
  if( debug > 0 ) Print(pBn,iccBn,fmt="6.3")
  ICRFP[,llBn]=iccBn;
  cnP[llBn]=colnames(iccBn)
 }
 if( nP > 0 ){
  # Generalized Partial Credit with step parameters and DinP
  pP=as.matrix( param[locP,3:ncp,drop=0] )
  rownames(pP)=iname[locP]
  res=icrfP( pP, theta, smallP=smallP, DinP=DinP )
  iccP=res$P
  if( debug > 0 ) Print(pP,iccP,fmt="6.3")
  ICRFP[,llP]=iccP;
  cnP[llP]=colnames(iccP)
 }
 if( nPN > 0 ){
  # Generalized Partial Credit in Nominal form wih 1.7
  pPN=as.matrix( param[locPN,3:ncp,drop=0] )
  rownames(pPN)=iname[locPN]
  res=icrfPN( pPN, theta, smallP=smallP )
  iccPN=res$P
  if( debug > 0 ) Print(pPN,iccPN,fmt="6.3")
  ICRFP[,llPN]=iccPN;
  cnP[llPN]=colnames(iccPN)
 }
 if( nG > 0 ){
  # Graded Response Model (2PLM)
  pG=as.matrix( param[locG,3:ncp,drop=0] )
  rownames(pG)=iname[locG]
  res=icrfG( pG, theta, smallP=smallP )
  iccG=res$P
  if( debug > 0 )Print(pG,iccG,fmt="6.3")
  ICRFP[,llG]=iccG
  cnP[llG]=colnames(iccG)
 }
 if( nGn > 0 ){
  # Graded Response Model (2PNM)
  pGn=as.matrix( param[locGn,3:ncp,drop=0] )
  rownames(pGn)=iname[locGn]
  res=icrfGn( pGn, theta, smallP=smallP )
  iccGn=res$P
  if( debug > 0 )Print(pGn,iccG,fmt="6.3")
  ICRFP[,llGn]=iccGn
  cnP[llGn]=colnames(iccGn)
 }
 if( nN > 0 ){
  # Nominal Response Model: *** not available yet ***
  pN=as.matrix( param[locN,3:ncp,drop=0] )
  rownames(pN)=iname[locN]
  res=icrfN( pN, theta, smallP=smallP )
  iccN=res$P
  if( debug > 0 )Print(pN,,iccN,fmt="6.3")
  ICRFP[,llN]=iccN;
  cnP[llN]=colnames(iccN)
 }
 colnames(ICRFP)=cnP; rownames(ICRFP)=format(theta,digits=3)


 # vector of category weights
 vecv=NULL
 for( i in 1:nitems ){
  vecv=c(vecv,v[i,1:ncat[i]])
 }
 vecv=as.matrix(vecv,,1)
 rownames(vecv)=cnP; colnames(vecv)=""
 # Print(vecv)

 # item response function and test response function
 ICRFPv=ICRFP*matrix(1,npoints)%*%t(vecv)
 # Print(ICRFPv)
 irf=matrix(0,npoints,nitems)
 for( i in 1:nitems ){
  irf[,i]=rowSums( ICRFPv[,fromP[i]:toP[i],drop=0] )
 }
 colnames(irf)=iname; rownames(irf)=format(theta,digits=3)
 trf=irf%*%w; colnames(trf)=""
 # Print(irf,trf)

 ncat2=ncat
 if( zero != 1 ){
  # remove the 0-th category
  ICRFP=ICRFP[,-fromP,drop=0] # 2017130
  fromP=fromP-0:(nitems-1)
  toP=toP-1:nitems
  ncat2=ncat-1
 }


 if( print ){
  cat("\n\ncalculation of item response functions etc \n")
  param=cbind(param,maxscore_i)
  weight=cbind(weight,maxscore_i)
  colnames(weight)[ncol(weight)]="maxscore_i"
  cat(" parameter data frame name =", pdfname
    , ",  item category weight data frame name =", wdfname,"\n")
  cat(" # of item parameters =",nitems
      , ", # of theta points =", npoints,"\n")
  cat(" include zero-th category =",zero,"\n")
  cat(" range of observed score = [",minscore_t,",", maxscore_t,"]\n")
  cat("\n item parameters with D in P =", DinP, "\n")
  print( param )
  cat("\n item and item category weight \n")
  print( weight )
  cat("\n")
  cat("\n ICRF\n")
  print(ICRFP,digits=3)
  cat("\n IRF\n")
  print(irf,digits=3)
  cat("\n TRF\n")
  print(trf,digits=3)
  # Print(ICRFP, irf, trf, digits=3)
 }

 if( plot ){

  if( is.null(plotchar) )
   plotchar=rep( c( 16,17,18,15, 21,24,23,22, 1,2,5,0 ),max(ncat) )
  else if( length(plotchar) == 1 ) plotchar=rep(plotchar,max(ncat))
  if( is.null(linetype) ) linetype=rep(1:3,max(ncat))
  else if( length(linetype) == 1 ) linetype=rep(linetype,max(ncat))
  if( is.null(colors) ) colors=rainbow(max(ncat2))
  else if( length(colors) == 1 ) colors=rep(colors,max(ncat2))

  for( j in 1:nitems ){

   # titles
   main=paste("Plot of ICRF of "
              ,iname[j]," ( type = ",type[j],", ncat=",ncat[j]," )")
   sub=paste("param = ", paste(format(param[j,4:(ncat[j]+3)],digits=3)
                               , collapse=",  "))
   if( length(grep("^B", param$type[j])) > 0 )
    sub=paste("param = ", paste(format(param[j,4:(ncat[j]+4)],digits=3)
                                , collapse=",  "))
   if( param$type[j] == "N" )
    sub=paste("param = ", paste(format(param[j,4:(2*(ncat[j]-1)+3)],digits=3)
                                , collapse=",  "))

   # set up the plot
   plot(range(theta), c(0,1), type="n", xlab="theta", ylab="icrf" )

   # add lines
   for (k in 1:ncat2[j]) {
    lines(theta, ICRFP[,fromP[j]+k-1], type="b", lwd=1.5,
          lty=linetype[k], col=colors[k], pch=plotchar[k])
   }

   # add a title and subtitle
   title(main, sub)

   # add a legend
   if( zero == 1 )
    legend( range(theta)[1], 0.8, (1:ncat2[j])-1, cex=0.8, col=colors
         ,  pch=plotchar, lty=linetype, title="cat" )
   else
    legend( range(theta)[1], 0.8, (1:ncat2[j]), cex=0.8, col=colors
            ,  pch=plotchar, lty=linetype, title="cat" )

  } # end of j


  if( plot >= 2 ){
   for( j in 1:(nitems) ){
    title=paste("Plot of IRF of "
                , iname[j]," ( type = ",type[j],", ncat=",ncat[j]
                , ", score range = [",minscore_i[j],",", maxscore_i[j],"] )")
    sub=paste("param = ", paste(format(param[j,4:(ncat[j]+3)],digits=3)
                                , collapse=",  "), "  (with weights)")
    if( length(grep("^B", param$type[j])) > 0 )
     sub=paste("param = ", paste(format(param[j,4:(ncat[j]+4)],digits=3)
                                 , collapse=",  "), "  (with weights)")
    plot(theta,irf[,j], main=title,sub=sub
         , xlim=c(min(theta),max(theta)), ylim=c(minscore_i[j],maxscore_i[j])
         , type="l", ylab="IRF")
   }
   title=paste("Plot of TRF  ( # of items =", nitems
               , ", score range = [",minscore_t,",", maxscore_t,"] )")
   plot(theta,trf, main=title
        , xlim=c(min(theta),max(theta)), ylim=c(minscore_t,maxscore_t)
        , type="l", ylab="TRF")
  }
 } # end of plot


 return( list(theta=theta, ICRF=ICRFP, IRF=irf, TRF=trf, fromP=fromP, toP=toP
            , vecv=vecv, minscore_i=minscore_i, maxscore_i=maxscore_i
            , minscore_t=minscore_t, maxscore_t=maxscore_t, zero=zero) )


} # end of irf










#' Calculation of Item Response Function of Binary Logistic Items
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
#'   P       length(theta) x nitems   ICRF matrix when zero=0 \cr
#'           length(theta) x 2*nitems ICRF matrix when zero=1
#'
#'
#' @export
#'
icrfB <- function( param, theta, maxZ=700, zero=0, smallP=0
                   , print=0, plot=0 ){
 # calculation of Binary Logistic ICRF
 # Shin-ichi Mayekawa
 # 120201
 # maxZ: 120202
 # zero: 120209
 # renamed as irfB -> icrfB: 120210
 # dataframe input: 120213,14
 # checkparam: 120224,29
 # when nrow(param) == 1: 120308
 # bugfix: 120908
 # when c-parameter is NA: 120921
 # smallP: 121022
 # diag(a) avoided: 20180208
 # digits -> fmt: 20180208
 #
 # Args:
 #  param    nitems x 3 item parameters matrix or data frame
 #  theta    npoints x 1 theta values
 #  zero     = 1 to include 1-P as the 0-th category probability
 #
 # Values:
 #   P       length(theta) x nitems   ICRF matrix when zero=0
 #           length(theta) x 2*nitems ICRF matrix when zero=1
 #
 #
 # Needs:
 #   checkparam
 #

 # argument name
 paramname=as.character(substitute(param))
 isdf=0
 if( is.data.frame(param) ){
  # param and weight given as data frame
  isdf=1
  param=checkparam( param, "B", "icrfB" )
  if( is.null(param) ){
   cat("\n\nInput Item Parameter ERROR.\n\n")
   return()
  }
 }

 #Print("after checkparam in ICRFB", param)

 # const
 nitems=nrow(param)
 iname=rownames(param)
 if( length(which(colnames(param) == "ncat")) > 0 )
  param=param[,-which(colnames(param) == "ncat"),drop=F]
 a=param[,1]; b=param[,2]; c=param[,3]

 # if the c-parameter is NA   120921
 c[is.na(c)]=0

 if( is.matrix(theta) ) theta=as.vector(theta)
 npoints=length(theta)

 Z=1.7*t(a*outer(b,theta,"-"))
 Z[which(Z > maxZ)]=maxZ
 Z[which(Z < -maxZ)]=-maxZ
 cc=matrix(c,npoints,nitems,byrow=1)
 P=cc+(1-cc)/(1+exp(Z))
 colnames(P)=rownames(param)
 if( zero == 1 ){
  P0=matrix(0,npoints,2*nitems)
  P0[,2*(1:nitems)]=P
  P0[,2*(1:nitems)-1]=1-P
  catname=outer(paste(iname,"_",sep=""),0:1,paste,sep="")
  colnames(P0)=matrix(t(catname),1)
  P=P0
 }
 P[which(P < smallP)]=smallP
 P[which(P > 1-smallP)]=1-smallP
 rownames(P)=format(theta,fmt="6.3")
 if( print > 0 ){
  cat("\nIRF of Binary Logistic Items \n ")
  if( isdf )
   cat(" item parameter data frame name =", paramname,"\n")
  else
   cat(" item parameter matrix name =", paramname,"\n")
  cat("  # of item parameters =",nitems,
      ",  # of theta points =", npoints,
      ",  zero category =",zero, "\n")
  Print(param)
  Print(P, fmt="6.3")
 }
 if( plot > 0 ){
  matplot(theta,P, type = "l", ylab="IRF"
          , xlim=c(min(theta),max(theta)), ylim=c(0,1)
          , main="Plot of ICRF of Binary Items")
 }
 return( P )

} # end of icrfB



#' Calculation of Item Response Function of Binary Norma CDF Items
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
#'   P       length(theta) x nitems   ICRF matrix when zero=0 \cr
#'           length(theta) x 2*nitems ICRF matrix when zero=1
#'
#'
#' @export
#'
icrfBn <- function( param, theta, maxZ=700, zero=0, smallP=0
                   , print=0, plot=0 ){
 # calculation of Binary Normal CDF ICRF
 # Shin-ichi Mayekawa
 # icrfB: 20120201-20180208
 # icrfB modified: 20180208
 #
 # Args:
 #  param    nitems x 3 item parameters matrix or data frame
 #  theta    npoints x 1 theta values
 #  zero     = 1 to include 1-P as the 0-th category probability
 #
 # Values:
 #   P       length(theta) x nitems   ICRF matrix when zero=0
 #           length(theta) x 2*nitems ICRF matrix when zero=1
 #
 #
 # Needs:
 #   checkparam
 #

 # argument name
 paramname=as.character(substitute(param))
 isdf=0
 if( is.data.frame(param) ){
  # param and weight given as data frame
  isdf=1
  param=checkparam( param, "Bn", "icrfBn" )
  if( is.null(param) ){
   cat("\n\nInput Item Parameter ERROR.\n\n")
   return()
  }
 }

 #Print("after checkparam in ICRFBn", param)

 # const
 nitems=nrow(param)
 iname=rownames(param)
 if( length(which(colnames(param) == "ncat")) > 0 )
  param=param[,-which(colnames(param) == "ncat"),drop=F]
 a=param[,1]; b=param[,2]; c=param[,3]

 # if the c-parameter is NA   120921
 c[is.na(c)]=0

 if( is.matrix(theta) ) theta=as.vector(theta)
 npoints=length(theta)

 Z=-t(a*outer(b,theta,"-"))
 Z[which(Z > maxZ)]=maxZ
 Z[which(Z < -maxZ)]=-maxZ
 cc=matrix(c,npoints,nitems,byrow=1)
 P=cc+(1-cc)*pnorm(Z)
 colnames(P)=rownames(param)
 if( zero == 1 ){
  P0=matrix(0,npoints,2*nitems)
  P0[,2*(1:nitems)]=P
  P0[,2*(1:nitems)-1]=1-P
  catname=outer(paste(iname,"_",sep=""),0:1,paste,sep="")
  colnames(P0)=matrix(t(catname),1)
  P=P0
 }
 P[which(P < smallP)]=smallP
 P[which(P > 1-smallP)]=1-smallP
 rownames(P)=format(theta,digits=3)
 if( print > 0 ){
  cat("\nIRF of Binary Normal CDF Items \n ")
  if( isdf )
   cat(" item parameter data frame name =", paramname,"\n")
  else
   cat(" item parameter matrix name =", paramname,"\n")
  cat("  # of item parameters =",nitems,
      ",  # of theta points =", npoints,
      ",  zero category =",zero, "\n")
  Print(param)
  Print(P, digits=3)
 }
 if( plot > 0 ){
  matplot(theta,P, type = "l", ylab="IRF"
          , xlim=c(min(theta),max(theta)), ylim=c(0,1)
          , main="Plot of ICRF of Binary Items")
 }
 return( P )

} # end of icrfBn



#' Calculation of Item Response Function of Graded Response Items
#'
#' @param param Item Parameter Data Frame
#' @param theta Discrete theta values
#' @param smallP Minimum value of probability
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#' @param debug = 1 to print intemediate result
#' @details
#'     An Example of Item Parameter datastet:  maxncat = 5\cr
#'     Note that the subscript of the b-parameter ranges between\cr
#'          1 and ncat-1 (= maxscore)  (We set c0=0.)\cr
#' \cr
#'     name  ncat     A            b1     b2     b3     b4\cr
#' \cr
#'     Item1  2     1.0           0.0     .      .      .\cr
#'     Item2  2     1.0           0.0     .      .      .\cr
#'     Item3  3     1.5          -1.0    0.0     .      .\cr
#'     Item4  4     1.5          -1.0    0.0    1.0     .\cr
#'     Item5  5     1.5          -1.0    0.0    1.0    2.0\cr
#'
#' @return
#' list \cr
#'   P             length(theta) x sum(ncat) ICRF matrix \cr
#'   PP            length(theta) x sum(ncat) ICBRF matrix \cr
#'  Pt            length(theta) x nitems item response function \cr
#'   fromP, toP    indexes of each ICRF in P \cr
#'   fromPP, toPP  indexes of each ICBRF in PP \cr
#'
#' @export
#'
icrfG <- function( param, theta, smallP=0, print=0, plot=0, debug=0 ){
 # calculation of GRM ICRF
 # Shin-ichi Mayekawa
 # 120201,02
 # bugfix: 120208
 # renamed: irfG -> icrfG : 120210
 # dataframe input: 120213,14,15
 # checkparam: 120224,29
 #
 # Args:
 #  param    nitems x max(ncat) matrix of item parameters or data frame
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
  param=checkparam( param, "G", "icrfG" )
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

 # locations of icrf and icbrf of each item in P and PP
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
 PP=matrix( 0,npoints, sum(ncat+1))
 P=matrix( 0,npoints, sum(ncat))
 Pt=matrix(0,npoints,nitems)
 rownames(PP)=thetac; colnames(PP)=PPname
 rownames(P)=thetac;  colnames(P)=Pname
 rownames(Pt)=thetac;  colnames(Pt)=iname

 if( debug ) print1=1
 else print1=0

 # for each item
 for( j in 1:nitems ){

  # binary logistic icrf
  paramj=cbind( rep(a[j],ncat1[j]), b[j,1:ncat1[j]],rep(0,ncat1[j]) )
  PPj=icrfB( paramj, theta, smallP=smallP, print=print1 )

  # icbrf
  PP[,fromPP[j]]=1;
  PP[,(fromPP[j]+1):(toPP[j]-1)]=PPj
  PP[,toPP[j]]=0;

  # icrf as difference of icbrf
  #    pj=PP[,fromPP[j]:toPP[j]-1]-PP[,fromPP[j]+1:toPP[j]];
  temp=PP[,fromPP[j]:(toPP[j]-1)]-PP[,(fromPP[j]+1):toPP[j]]

  temp[which(temp < smallP)]=smallP
  temp[which(temp > 1-smallP)]=1-smallP

  P[,fromP[j]:toP[j]]=temp

  # icrf
  Pt[,j]=rowSums(P[,fromP[j]:toP[j]]%*%diag(c(0:(ncat[j]-1))))

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
  cat("\nItem Category Response Functions (ICRF)\n")
  Print(P, fmt="5.3")
  cat("\nItem Response Functions (IRF) with natural category weights\n")
  Print(Pt, fmt="5.3")
  if( print >= 2 ){
   cat("\nItem Category Boundary Response Functions (ICBRF)\n")
   Print(PP, fmt="5.3")
  }
 }

 if( plot ){
  for( j in 1:(nitems) ){
   title=paste("Plot of ICRF of ",iname[j]," (ncat=",ncat[j],")")
   sub=paste("param = ", paste(param[j,1:ncat[j]+1], collapse=",  "))
   for( k in 1:ncat1[j] ){
    plot(theta,P[,fromP[j]+k-1]
        , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="")
    par(new=1)
   }
   plot(theta,P[,toP[j]], main=title, sub=sub
       , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="ICRF")
   par(new=0)
  }
  if( plot >= 2 ){
   for( j in 1:(nitems) ){
    title=paste("Plot of IRF (with natural category weight) of "
                , iname[j]," (ncat=",ncat[j],")")
    sub=paste("param = ", paste(param[j,1:ncat[j]+1], collapse=",  "))
    plot(theta,Pt[,j], main=title, sub=sub
         , xlim=c(min(theta),max(theta)), ylim=c(0,ncat1[j])
         , type="l", ylab="IRF")
    par(new=0)
   }
  }
 }
 if( plot >= 3 ){
  for( j in 1:(nitems) ){
   title=paste("Plot of ICBRF of ",iname[j]," (ncat=",ncat[j],")")
   sub=paste("param = ", paste(param[j,1:ncat[j]+1], collapse=",  "))
   for( k in 1:ncat[j] ){
    plot(theta,PP[,fromPP[j]+k-1]
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="")
    par(new=1)
   }
   plot(theta,PP[,toPP[j]], main=title, sub=sub
        , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="ICBRF")
   par(new=0)
  }
 }

 return( list(P=P, PP=PP, Pt=Pt
              , fromP=fromP, toP=toP, fromPP=fromPP, toPP=toPP) )

} # end of icrfG



#' Calculation of Item Response Function of Graded Response Items
#'  with Normal Ogive Model (2PNM)
#'
#' @param param Item Parameter Data Frame
#' @param theta Discrete theta values
#' @param smallP Minimum value of probability
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#' @param debug = 1 to print intemediate result
#' @details
#'     An Example of Item Parameter datastet:  maxncat = 5\cr
#'     Note that the subscript of the b-parameter ranges between\cr
#'          1 and ncat-1 (= maxscore)  (We set c0=0.)\cr
#' \cr
#'     name  ncat     A            b1     b2     b3     b4\cr
#' \cr
#'     Item1  2     1.0           0.0     .      .      .\cr
#'     Item2  2     1.0           0.0     .      .      .\cr
#'     Item3  3     1.5          -1.0    0.0     .      .\cr
#'     Item4  4     1.5          -1.0    0.0    1.0     .\cr
#'     Item5  5     1.5          -1.0    0.0    1.0    2.0\cr
#'
#' @return
#' list \cr
#'   P             length(theta) x sum(ncat) ICRF matrix \cr
#'   PP            length(theta) x sum(ncat) ICBRF matrix \cr
#'  Pt            length(theta) x nitems item response function \cr
#'   fromP, toP    indexes of each ICRF in P \cr
#'   fromPP, toPP  indexes of each ICBRF in PP \cr
#'
#' @export
#'
icrfGn <- function( param, theta, smallP=0, print=0, plot=0, debug=0 ){
 # calculation of GRM ICRF
 # Shin-ichi Mayekawa
 # icrfG: 20120201-20120229
 # icrfG modified: 20180209
 #
 # Args:
 #  param    nitems x max(ncat) matrix of item parameters or data frame
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
  param=checkparam( param, "Gn", "icrfGn" )
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

 # locations of icrf and icbrf of each item in P and PP
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
 PP=matrix( 0,npoints, sum(ncat+1))
 P=matrix( 0,npoints, sum(ncat))
 Pt=matrix(0,npoints,nitems)
 rownames(PP)=thetac; colnames(PP)=PPname
 rownames(P)=thetac;  colnames(P)=Pname
 rownames(Pt)=thetac;  colnames(Pt)=iname

 if( debug ) print1=1
 else print1=0

 # for each item
 for( j in 1:nitems ){

  # binary logistic icrf
  paramj=cbind( rep(a[j],ncat1[j]), b[j,1:ncat1[j]],rep(0,ncat1[j]) )
  PPj=icrfBn( paramj, theta, smallP=smallP, print=print1 )

  # icbrf
  PP[,fromPP[j]]=1;
  PP[,(fromPP[j]+1):(toPP[j]-1)]=PPj
  PP[,toPP[j]]=0;

  # icrf as difference of icbrf
  #    pj=PP[,fromPP[j]:toPP[j]-1]-PP[,fromPP[j]+1:toPP[j]];
  temp=PP[,fromPP[j]:(toPP[j]-1)]-PP[,(fromPP[j]+1):toPP[j]]

  temp[which(temp < smallP)]=smallP
  temp[which(temp > 1-smallP)]=1-smallP

  P[,fromP[j]:toP[j]]=temp

  # icrf
  Pt[,j]=rowSums(P[,fromP[j]:toP[j]]%*%diag(c(0:(ncat[j]-1))))

 } # end of j loop


 if( print > 0 ){
  cat("\nICRF of Graded Response Items (2PNM) \n ")
  if( isdf )
   cat(" item parameter data frame name =", paramname,"\n")
  else
   cat(" item parameter matrix name =", paramname,"\n")
  cat("  # of item parameters =",nitems,
      ",  # of theta points =", npoints,"\n")
  Print(param)
  cat("\nItem Category Response Functions (ICRF)\n")
  Print(P, fmt="5.3")
  cat("\nItem Response Functions (IRF) with natural category weights\n")
  Print(Pt, fmt="5.3")
  if( print >= 2 ){
   cat("\nItem Category Boundary Response Functions (ICBRF)\n")
   Print(PP, fmt="5.3")
  }
 }

 if( plot ){
  for( j in 1:(nitems) ){
   title=paste("Plot of ICRF of ",iname[j]," (ncat=",ncat[j],")")
   sub=paste("param = ", paste(param[j,1:ncat[j]+1], collapse=",  "))
   for( k in 1:ncat1[j] ){
    plot(theta,P[,fromP[j]+k-1]
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="")
    par(new=1)
   }
   plot(theta,P[,toP[j]], main=title, sub=sub
        , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="ICRF")
   par(new=0)
  }
  if( plot >= 2 ){
   for( j in 1:(nitems) ){
    title=paste("Plot of IRF (with natural category weight) of "
                , iname[j]," (ncat=",ncat[j],")")
    sub=paste("param = ", paste(param[j,1:ncat[j]+1], collapse=",  "))
    plot(theta,Pt[,j], main=title, sub=sub
         , xlim=c(min(theta),max(theta)), ylim=c(0,ncat1[j])
         , type="l", ylab="IRF")
    par(new=0)
   }
  }
 }
 if( plot >= 3 ){
  for( j in 1:(nitems) ){
   title=paste("Plot of ICBRF of ",iname[j]," (ncat=",ncat[j],")")
   sub=paste("param = ", paste(param[j,1:ncat[j]+1], collapse=",  "))
   for( k in 1:ncat[j] ){
    plot(theta,PP[,fromPP[j]+k-1]
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="")
    par(new=1)
   }
   plot(theta,PP[,toPP[j]], main=title, sub=sub
        , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="ICBRF")
   par(new=0)
  }
 }

 return( list(P=P, PP=PP, Pt=Pt
              , fromP=fromP, toP=toP, fromPP=fromPP, toPP=toPP) )

} # end of icrfGn





#' Calculation of Item Response Function of Partial Credit Items \cr
#' in Nominal Model Format
#'
#' @param paramPN Item Parameter Data Frame of Partial Credit Items
#'  in Nominal format
#' @param theta Discrete theta values
#' @param smallP Minimum value of probability
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#' @param debug = 1 to print intemediate result
#'
#' @details
#'     An Example of Item Parameter datastet:  maxncat = 5\cr
#'     Note that the subscript of the b-parameter ranges between\cr
#'          1 and ncat-1 (= maxscore)  (We set c0=0.)\cr
#' \cr
#'           ncat     A            b1     b2     b3     b4\cr
#' \cr
#'     Item1  2     1.0           0.0     .      .      .\cr
#'     Item2  2     1.0           0.0     .      .      .\cr
#'     Item3  3     1.5          -1.0    0.0     .      .\cr
#'     Item4  4     1.5          -1.0    0.0    1.0     .\cr
#'     Item5  5     1.5          -1.0    0.0    1.0    2.0\cr
#' \cr
#'     ICRF of category k of item j at theta, namely, \eqn{p_{jk}(theta)}\cr
#'     is defined as\cr
#'      \eqn{p_{jk}(theta) = Ez_{jk}(theta)
#'      / sum_{k=0}^{ncat[j]-1} Ez_{jk}(theta)}\cr
#'       where\cr
#'      \eqn{Ez_{jk}(theta) = exp( 1.7 z_{jk}(theta) )}\cr
#'        and\cr
#'      \eqn{z_{jk}(theta) = a_j k (theta - b_{k}), k=0,1, ..., ncat[j]-1}\cr
#'     with \eqn{b_{j0}=0}.\cr
#' \cr
#'
#' @return A list of: \cr
#'   P       length(theta) x sum(ncat) icrf matrix \cr
#'   PP      length(theta) x sum(ncat) icbf  matrix \cr
#'   Pt      length(theta) x nitems item response function \cr

#'
#' @export
#'
icrfPN <- function( paramPN, theta, smallP=0, print=0, plot=0, debug=0 ){
 # calculation of Generalized Partial Credit ICRF in Nominal format
 # Shin-ichi Mayekawa
 # 120201,02
 # renamed: irfP -> irfPN : 120208
 # renamed: irfPN -> icrfPN : 120210
 # dataframe input: 120213,14,15
 # checkparam: 120224,29
 # converted to b-type paramters input: 120301
 # b-type paramters corrected: 120307
 # smallP: 121022
 #
 #
 # Args:
 #  paramPN    nitems x max(ncat) matrix of item parameters or data frame
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
 #     Ez_{jk}(theta) = exp( 1.7 z_{jk}(theta) )
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
 #    checkparam
 #

 # argument name
 paramname=as.character(substitute(paramPN))
 isdf=0
 if( is.data.frame(paramPN) ){
  # param and weight given as data frame
  isdf=1
  paramPN=checkparam( paramPN, "PN", "icrfPN" )
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

 # locations of icrf and icbrf of each item in P and PP
 toP=cumsum(ncat)
 fromP=toP-ncat+1
 toPP=cumsum(ncat+1)
 fromPP=toPP-ncat

 # names of icrf and icbrf
 PPname=character(sum(ncat+1))
 Pname=character(sum(ncat))
 catname=outer(paste(iname,"_",sep=""),0:max(ncat),paste,sep="")

 for( j in 1:nitems ){
  PPname[fromPP[j]:toPP[j]]=catname[j,1:(ncat[j]+1)]
  Pname[fromP[j]:toP[j]]=catname[j,1:ncat[j]]
 }
 rm(catname)


 # main body
 if( debug ) Print(a,c,ncat,ncat1)
 P=matrix( 0,npoints, sum(ncat))
 Pt=matrix(0,npoints,nitems)
 PP=matrix( 0,npoints, sum(ncat+1))
 rownames(PP)=thetac; colnames(PP)=PPname
 rownames(P)=thetac;  colnames(P)=Pname
 rownames(Pt)=thetac;  colnames(Pt)=iname

 if( debug ) print1=1
 else print1=0

 ########## testing new type-b item params ########
 # P1=matrix( 0,npoints, sum(ncat))
 # a1=paramPN[,2,drop=1]
 # b1=paramPN[,3:(1+max(ncat)),drop=0]
 ##################################################
 # calculation of icrf
 # for each item
 for( j in 1:nitems ){

  # exp( linear components ) adjusted to sum to unity
  cj=c[j,1:ncat1[j],drop=0];
  beta=a[j]*((1:ncat[j])-1)
  alpha=cbind(0 , cj);
  if( debug ) Print(cj,beta,alpha)
  # z=theta*beta+j(npoints,1,1)*alpha;
  Z=as.matrix( outer(theta,beta,"*")+matrix(1,npoints)%*%alpha )
  Pj=exp(Z);
  if( debug ) Print(Z, Pj, rowSums(Pj))
  if( npoints > 1 ) Pj=diag(1/rowSums(Pj))%*%Pj
  else Pj=Pj/rowSums(Pj)

  Pj[which(Pj < smallP)]=smallP
  Pj[which(Pj > 1-smallP)]=1-smallP

  P[,fromP[j]:toP[j]]=Pj;

  # icrf
  Pt[,j]=rowSums(P[,fromP[j]:toP[j]]%*%diag(c(0:(ncat[j]-1))))

  ########## testing new type-b item params ########
  #   ajk=matrix(a1[j]*(1:(ncat[j]-1)),npoints,ncat[j]-1,byrow=1)
  #   Pj=exp( 1.7*ajk*(theta-matrix(1,npoints)%*%b1[j,1:(ncat[j]-1)]) )
  #   Pj=cbind(1,Pj)
  #   if( npoints > 1 ) Pj=diag(1/rowSums(Pj))%*%Pj
  #   else Pj=Pj/rowSums(Pj)
  #   #  Print(j,Pj,fromP[j],toP[j])
  #   P1[,fromP[j]:toP[j]]=Pj;
  ##################################################

 } # end of j loop


 # calculation of icbrf
 for( j in 1:nitems ){
  for( k in ncat[j]:1 ){
   #   PP[,fromtoPP[j,1]-1 +k]=
   #    P[,fromtoP[j,1]-1 +k]+PP[,fromtoPP[j,1]-1 +k+1];
   PP[,fromPP[j]-1 + k]=P[,fromP[j]-1 + k] + PP[,fromPP[j]-1 + k+1]
  }
 }


 if( print > 0 ){
  cat("\nICRF of Generalized Partial Credit Items \n ")
  if( isdf )
   cat(" item parameter data frame name =", paramname,"\n")
  else
   cat(" item parameter matrix name =", paramname,"\n")
  cat("  # of item parameters =",nitems,
      ",  # of theta points =", npoints,"\n")
  cat(" Item Parameters (type b)\n")
  Print(paramPN)
  cat("\nItem Category Response Functions (ICRF)\n")
  Print(P, digits=2)
  cat("\nItem Response Functions (IRF) with natural category weights\n")
  Print(Pt, digits=2)
  if( print >= 2 ){
   cat("\nItem Category Boundary Response Functions (ICBRF)\n")
   Print(PP,digits=2)
  }
 }

 if( plot ){
  for( j in 1:(nitems) ){
   title=paste("Plot of ICRF of ",iname[j]," (ncat=",ncat[j],")")
   sub=paste("param = ", paste(paramPN[j,1:ncat[j]+1], collapse=",  "))
   for( k in 1:ncat1[j] ){
    plot(theta,P[,fromP[j]+k-1]
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="")
    par(new=1)
   }
   plot(theta,P[,toP[j]], main=title, sub=sub
        , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="ICRF")
   par(new=0)
  }
  if( plot >= 2 ){
   for( j in 1:(nitems) ){
    title=paste("Plot of IRF (with natural category weight) of "
                , iname[j]," (ncat=",ncat[j],")")
    sub=paste("param = ", paste(paramPN[j,1:ncat[j]+1], collapse=",  "))
    plot(theta,Pt[,j], main=title, sub=sub
         , xlim=c(min(theta),max(theta)), ylim=c(0,ncat1[j])
         , type="l", ylab="IRF")
    par(new=0)
   }
  }
 }
 if( plot >= 3 ){
  for( j in 1:(nitems) ){
   title=paste("Plot of ICBRF of ",iname[j]," (ncat=",ncat[j],")")
   sub=paste("param = ", paste(paramPN[j,1:ncat[j]+1], collapse=",  "))
   for( k in 1:ncat[j] ){
    plot(theta,PP[,fromPP[j]+k-1]
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="")
    par(new=1)
   }
   plot(theta,PP[,toPP[j]], main=title, sub=sub
        , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="ICBRF")
   par(new=0)
  }
 }

 return( list(P=P, PP=PP, Pt=Pt
              , fromP=fromP, toP=toP, fromPP=fromPP, toPP=toPP) )

} # end of icrfPN


#' Calculation of Item Response Function of Partial Credit Items
#'
#' @param paramP Item Parameter Data Frame of Partial Credit Items
#' @param theta Discrete theta values
#' @param DinP = 1 to include D=1.7 in logistic function
#' @param smallP Minimum value of probability
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#' @param debug = 1 to print intemediate result
#'
#' @details
#'     An Example of Item Parameter datastet:  maxncat = 5\cr
#'     Note that the subscript of the b-parameter ranges between\cr
#'          1 and ncat-1 (= maxscore)  (We set c0=0.)\cr
#' \cr
#'           ncat     A            b1     b2     b3     b4\cr
#' \cr
#'     Item3  3     1.5          -1.0    0.0     .      .\cr
#'     Item4  4     1.5          -1.0    0.0    1.0     .\cr
#'     Item5  5     1.5          -1.0    0.0    1.0    2.0\cr
#' \cr
#'     ICRF of category k of item j at theta, namely, \eqn{p_{jk}(theta)}\cr
#'     is defined as\cr
#'      \eqn{p_{jk}(theta) = Ez_{jk}(theta)
#'      / sum_{k=0}^{ncat[j]-1} Ez_{jk}(theta)}\cr
#'       where\cr
#'  \eqn{Ez_{jk}(theta) = exp( 1.7^DinP a*_j k (theta - sum_{h=0}^k b*_{jh}) )}
#'      \cr
#'               \eqn{= exp( 1.7^DinP a*_j k (theta - sum_{h=1}^k b*_{jh}) )}
#'                     \cr
#'                                      \eqn{, k=0,1, ..., ncat[j]-1}\cr
#'     with \eqn{b*_{j0}=0}.\cr
#' \cr
#'     This \eqn{b*_{jk}} is the original step parameter and it is\cr
#'     the value of theta where \eqn{P_{jk-1}} and \eqn{P_{jk}} intersect.\cr
#'
#' @return A list of: \cr
#'   P       length(theta) x sum(ncat) icrf matrix \cr
#'   PP      length(theta) x sum(ncat) icbf  matrix \cr
#'   Pt      length(theta) x nitems item response function \cr
#'
#' @export
#'
icrfP <- function( paramP, theta, DinP=1
                 , smallP=0, print=0, plot=0, debug=0 ){
 # calculation of Generalized Partial Credit ICRF with original parametrization
 # Shin-ichi Mayekawa
 # converted from icrfPN: 121016
 # smallP: 121022
 # DinP option:121109(London)
 #
 #
 # Args:
 #  paramP    nitems x max(ncat) matrix of item parameters or data frame
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
 #     Ez_{jk}(theta) = exp( 1.7^DinP a*_j k (theta - sum_{h=0}^k b*_{jh}) )
 #                    = exp( 1.7^DinP a*_j k (theta - sum_{h=1}^k b*_{jh}) )
 #                                              , k=0,1, ..., ncat[j]-1
 #    with b*_{j0}=0.
 #
 #    This b*_{jk} is the original step parameter and it is
 #    the value of theta where P_{jk-1} and P_{jk} intersect.
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
 #    checkparam
 #

 # argument name
 paramname=as.character(substitute(paramP))
 isdf=0
 if( is.data.frame(paramP) ){
  # param and weight given as data frame
  isdf=1
  paramP=checkparam( paramP, "P", "icrfP" )
  if( is.null(paramP) ){
   cat("\n\nInput Item Parameter ERROR.\n\n")
   return()
  }
 }

 # convert to c-type params (121016)
 paramP0=paramP
 if( DinP == 1) paramP0[,2]=1.7*paramP0[,2]
 bs=paramP[,3:ncol(paramP),drop=0]
 c=t( apply(bs,1,cumsum) )
 paramP0[,3:ncol(paramP)]=-paramP0[,2]*c
 colnames(paramP0)=c("ncat","a",paste("c",1:(ncol(paramP0)-2),sep=""))

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

 # locations of icrf and icbrf of each item in P and PP
 toP=cumsum(ncat)
 fromP=toP-ncat+1
 toPP=cumsum(ncat+1)
 fromPP=toPP-ncat

 # names of icrf and icbrf
 PPname=character(sum(ncat+1))
 Pname=character(sum(ncat))
 catname=outer(paste(iname,"_",sep=""),0:max(ncat),paste,sep="")

 for( j in 1:nitems ){
  PPname[fromPP[j]:toPP[j]]=catname[j,1:(ncat[j]+1)]
  Pname[fromP[j]:toP[j]]=catname[j,1:ncat[j]]
 }
 rm(catname)


 # main body
 if( debug ) Print(a,c,ncat,ncat1)
 P=matrix( 0,npoints, sum(ncat))
 Pt=matrix(0,npoints,nitems)
 PP=matrix( 0,npoints, sum(ncat+1))
 rownames(PP)=thetac; colnames(PP)=PPname
 rownames(P)=thetac;  colnames(P)=Pname
 rownames(Pt)=thetac;  colnames(Pt)=iname

 if( debug ) print1=1
 else print1=0

 # calculation of icrf
 # for each item
 for( j in 1:nitems ){

  # exp( linear components ) adjusted to sum to unity
  cj=c[j,1:ncat1[j],drop=0];
  beta=a[j]*((1:ncat[j])-1)
  alpha=cbind(0 , cj);
  if( debug ) Print(cj,beta,alpha)
  # z=theta*beta+j(npoints,1,1)*alpha;
  Z=as.matrix( outer(theta,beta,"*")+matrix(1,npoints)%*%alpha )
  Pj=exp(Z);
  if( debug ) Print(Z, Pj, rowSums(Pj))
  if( npoints > 1 ) Pj=diag(1/rowSums(Pj))%*%Pj
  else Pj=Pj/rowSums(Pj)

  Pj[which(Pj < smallP)]=smallP
  Pj[which(Pj > 1-smallP)]=1-smallP

  P[,fromP[j]:toP[j]]=Pj;

  # icrf
  Pt[,j]=rowSums(P[,fromP[j]:toP[j]]%*%diag(c(0:(ncat[j]-1))))

 } # end of j loop


 # calculation of icbrf
 for( j in 1:nitems ){
  for( k in ncat[j]:1 ){
   #   PP[,fromtoPP[j,1]-1 +k]=
   #    P[,fromtoP[j,1]-1 +k]+PP[,fromtoPP[j,1]-1 +k+1];
   PP[,fromPP[j]-1 + k]=P[,fromP[j]-1 + k] + PP[,fromPP[j]-1 + k+1]
  }
 }


 if( print > 0 ){
  cat("\nICRF of Generalized Partial Credit Items (original and DinP) \n ")
  if( isdf )
   cat(" item parameter data frame name =", paramname,"\n")
  else
   cat(" item parameter matrix name =", paramname,"\n")
  cat("  # of item parameters =",nitems,
      ",  # of theta points =", npoints,"\n")
  cat("\n item parameters with D in P =", DinP, "\n")
  Print(paramP)
  cat("\nItem Category Response Functions (ICRF)\n")
  Print(P, digits=2)
  cat("\nItem Response Functions (IRF) with natural category weights\n")
  Print(Pt, digits=2)
  if( print >= 2 ){
   cat("\nItem Category Boundary Response Functions (ICBRF)\n")
   Print(PP,digits=2)
  }
 }

 if( plot ){
  for( j in 1:(nitems) ){
   title=paste("Plot of ICRF of ",iname[j]," (ncat=",ncat[j],")")
   sub=paste("param = ", paste(paramP[j,1:ncat[j]+1], collapse=",  "))
   for( k in 1:ncat1[j] ){
    plot(theta,P[,fromP[j]+k-1]
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="")
    par(new=1)
   }
   plot(theta,P[,toP[j]], main=title, sub=sub
        , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="ICRF")
   par(new=0)
  }
  if( plot >= 2 ){
   for( j in 1:(nitems) ){
    title=paste("Plot of IRF (with natural category weight) of "
                , iname[j]," (ncat=",ncat[j],")")
    sub=paste("param = ", paste(paramP[j,1:ncat[j]+1], collapse=",  "))
    plot(theta,Pt[,j], main=title, sub=sub
         , xlim=c(min(theta),max(theta)), ylim=c(0,ncat1[j])
         , type="l", ylab="IRF")
    par(new=0)
   }
  }
 }
 if( plot >= 3 ){
  for( j in 1:(nitems) ){
   title=paste("Plot of ICBRF of ",iname[j]," (ncat=",ncat[j],")")
   sub=paste("param = ", paste(paramP[j,1:ncat[j]+1], collapse=",  "))
   for( k in 1:ncat[j] ){
    plot(theta,PP[,fromPP[j]+k-1]
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="")
    par(new=1)
   }
   plot(theta,PP[,toPP[j]], main=title, sub=sub
        , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="ICBRF")
   par(new=0)
  }
 }

 return( list(P=P, PP=PP, Pt=Pt
              , fromP=fromP, toP=toP, fromPP=fromPP, toPP=toPP) )

} # end of icrfP








#' Calculation of  Partial Credit ICRF in Nominal format
#'
#' @param param Item Parameter Data Frame
#' @param theta Discrete theta values
#' @param smallP Minimum value of probability
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#' @param debug = 1 to print intemediate result
#'
#' @details
#'     An Example of Item Parameter datastet:  maxncat = 5\cr
#'     Note that the subscript of the c-parameter ranges between\cr
#'          1 and ncat-1 (= maxscore)  (We set c0=0.)\cr
#' \cr
#'           ncat     A            c1     c2     c3     c4\cr
#' \cr
#'     Item1  2     1.0           0.0     .      .      .\cr
#'     Item2  2     1.0           0.0     .      .      .\cr
#'     Item3  3     1.5          -1.0    0.0     .      .\cr
#'     Item4  4     1.5          -1.0    0.0    1.0     .\cr
#'     Item5  5     1.5          -1.0    0.0    1.0    2.0\cr
#' \cr
#'     ICRF of category k of item j at theta, namely, \eqn{p_{jk}(theta)}\cr
#'     is defined as\cr
#'      \eqn{p_{jk}(theta) = Ez_{jk}(theta)
#'      / sum_{k=0}^{ncat[j]-1} Ez_{jk}(theta)}\cr
#'       where\cr
#'      \eqn{Ez_{jk}(theta) = exp( z_{jk}(theta) )}\cr
#'        and\cr
#'      \eqn{z_{jk}(theta) = a_j k theta + c_{k}, k=0,1, ..., ncat[j]-1}\cr
#'     with \eqn{c_{j0}=0}.\cr
#'
#' @return A list of \cr
#'   P             length(theta) x sum(ncat) ICRF matrix \cr
#'   PP            length(theta) x sum(ncat) ICBRF matrix \cr
#'   Pt            length(theta) x nitems item response function \cr
#'   fromP, toP    indexes of each ICRF in P \cr
#'   fromPP, toPP  indexes of each ICBRF in PP
#'
#'
#'
#' @export
#'
icrfPN0 <- function( param, theta, smallP=0, print=0, plot=0, debug=0 ){
 # calculation of Generalized Partial Credit ICRF in Nominal format
 # Shin-ichi Mayekawa
 # 120201,02
 # renamed: irfP -> irfPN : 120208
 # renamed: irfPN -> icrfPN : 120210
 # dataframe input: 120213,14,15
 # checkparam: 120224,29
 # renamed as icrfPN0: 120301
 # plot: 121123
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
 # Needs:
 #    checkparam
 #

 # argument name
 paramname=as.character(substitute(param))
 isdf=0
 if( is.data.frame(param) ){
  # param and weight given as data frame
  isdf=1
  param=checkparam( param, "PN", "icrfPN" )
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
 a=param[,2,drop=0]
 c=param[,3:(1+max(ncat)),drop=0]

 # locations of icrf and icbrf of each item in P and PP
 toP=cumsum(ncat)
 fromP=toP-ncat+1
 toPP=cumsum(ncat+1)
 fromPP=toPP-ncat

 # names of icrf and icbrf
 PPname=character(sum(ncat+1))
 Pname=character(sum(ncat))
 catname=outer(paste(iname,"_",sep=""),0:max(ncat),paste,sep="")
 for( j in 1:nitems ){
  PPname[fromPP[j]:toPP[j]]=catname[j,1:(ncat[j]+1)]
  Pname[fromP[j]:toP[j]]=catname[j,1:ncat[j]]
 }
 rm(catname)


 # main body
 if( debug ) Print(a,c,ncat,ncat1)
 P=matrix( 0,npoints, sum(ncat))
 Pt=matrix(0,npoints,nitems)
 PP=matrix( 0,npoints, sum(ncat+1))
 rownames(PP)=thetac; colnames(PP)=PPname
 rownames(P)=thetac;  colnames(P)=Pname
 rownames(Pt)=thetac;  colnames(Pt)=iname

 if( debug ) print1=1
 else print1=0

 # calculation of icrf
 # for each item
 for( j in 1:nitems ){

  # exp( linear components ) adjusted to sum to unity
  cj=c[j,1:ncat1[j],drop=0];
  beta=a[j]*((1:ncat[j])-1)
  alpha=cbind(0 , cj);
  if( debug ) Print(cj,beta,alpha)
  # z=theta*beta+j(npoints,1,1)*alpha;
  Z=as.matrix( outer(theta,beta,"*")+matrix(1,npoints)%*%alpha )
  Pj=exp(Z);
  if( debug ) Print(Z, Pj, rowSums(Pj))
  if( npoints > 1 ) Pj=diag(1/rowSums(Pj))%*%Pj
  else Pj=Pj/rowSums(Pj)

  Pj[which(Pj < smallP)]=smallP
  Pj[which(Pj > 1-smallP)]=1-smallP

  P[,fromP[j]:toP[j]]=Pj;

  # icrf
  Pt[,j]=rowSums(P[,fromP[j]:toP[j]]%*%diag(c(0:(ncat[j]-1))))

 } # end of j loop

 # calculation of icbrf
 for( j in 1:nitems ){
  for( k in ncat[j]:1 ){
   #   PP[,fromtoPP[j,1]-1 +k]=
   #    P[,fromtoP[j,1]-1 +k]+PP[,fromtoPP[j,1]-1 +k+1];
   PP[,fromPP[j]-1 + k]=P[,fromP[j]-1 + k] + PP[,fromPP[j]-1 + k+1]
  }
 }


 if( print > 0 ){
  cat("\nICRF of Generalized Partial Credit Items \n ")
  if( isdf )
   cat(" item parameter data frame name =", paramname,"\n")
  else
   cat(" item parameter matrix name =", paramname,"\n")
  cat("  # of item parameters =",nitems,
      ",  # of theta points =", npoints,"\n")
  Print(param)
  cat("\nItem Category Response Functions (ICRF)\n")
  Print(P, digits=2)
  cat("\nItem Response Functions (IRF) with natural category weights\n")
  Print(Pt, digits=2)
  if( print >= 2 ){
   cat("\nItem Category Boundary Response Functions (ICBRF)\n")
   Print(PP,digits=2)
  }
 }

 if( plot ){
  for( j in 1:(nitems) ){
   title=paste("Plot of ICRF of ",iname[j]," (ncat=",ncat[j],")")
   sub=paste("param = ", paste(param[j,1:ncat[j]+1], collapse=",  "))
   for( k in 1:ncat1[j] ){
    plot(theta,P[,fromP[j]+k-1]
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="")
    par(new=1)
   }
   plot(theta,P[,toP[j]], main=title, sub=sub
        , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="ICRF")
   par(new=0)
  }
  if( plot >= 2 ){
   for( j in 1:(nitems) ){
    title=paste("Plot of IRF (with natural category weight) of "
                , iname[j]," (ncat=",ncat[j],")")
    sub=paste("param = ", paste(param[j,1:ncat[j]+1], collapse=",  "))
    plot(theta,Pt[,j], main=title, sub=sub
         , xlim=c(min(theta),max(theta)), ylim=c(0,ncat1[j])
         , type="l", ylab="IRF")
    par(new=0)
   }
  }
 }
 if( plot >= 2 ){
  for( j in 1:(nitems) ){
   title=paste("Plot of ICBRF of ",iname[j]," (ncat=",ncat[j],")")
   sub=paste("param = ", paste(param[j,1:ncat[j]+1], collapse=",  "))
   for( k in 1:ncat[j] ){
    plot(theta,PP[,fromPP[j]+k-1]
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="")
    par(new=1)
   }
   plot(theta,PP[,toPP[j]], main=title, sub=sub
        , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="ICBRF")
   par(new=0)
  }
 }

 return( list(P=P, PP=PP, Pt=Pt
              , fromP=fromP, toP=toP, fromPP=fromPP, toPP=toPP) )

} # end of icrfPN0





#' Generation of Simulated Item Response Data
#'
#' @param Ntot # of total observations to be generated
#'            If Ntot == 1, one observation per each theta point
#'            will be generated, \cr
#'            otherwize, the # of observations per each theta point is \cr
#'            proportional to round(Ntot*thd) where thd is \cr
#'            the npoint x 1 vector of theta distribution. \cr
#'            Therefore, Ntot must be large enough. \cr
#' @param param Item Parameter Data Frame
#' @param DinP = 1 to include D=1.7 in logistic function
#' @param Nmat npoints x nitems matrix consisting of # of trials, or NULL.
#' This will be used only when Ntot=1 and compress=0.
#' @param zero = 0 to exclude zero-th category when not compressed.
#' @param theta Discrete theta values
#' @param npoints # of discrete points for theta, or # of random variables
#' @param thmin Minimum value of discrete thata value
#' @param thmax Maximum value of discrete thata value
#' @param thd theta distribution probability vector
#' @param thdist = "NORMAL" or "UNIFORM" or "rnorm"  or  "runif" \cr
#'            When "NORMAL" or "rnorm", \cr
#'                thmean and thstd will be used to generate thd. \cr
#'            When "UNIFORM" or "runif", \cr
#'                thmin and thmax will be used to generate thd. \cr
#'            When "rnorm" or "runif", theta will be generated using \cr
#'                npoints random numbers and thd  and Ntot is set equal to 1.
#' @param thmean mean of theta distribution
#' @param thstd std of theta distribution
#' @param nomiss = 1 to exclude those theta points with N = 0.
#' @param compress = 1 to compress the data when Ntot = 1.
#' @param sort = 1 to sort the result according to theta.
#' @param thetaname = Name of the variable which will contain theta value
#' or NULL.
#'
#' @return A list of \cr
#'  U, N, npoints, theta, thd, Ntot, fromP, toP, type, ncat, thmean, thstd
#'  , compress, sort
#'
#' \cr
#' where
#' \cr
#' \cr
#'            If compress = 0 and Ntot is large, \cr
#'              U    is   npoints x sum of ncat[j] \cr
#'               U[,fromP[j]:toP[j]]  is npoints x ncat[j] \cr
#'               U[,fromP[j]:toP[j]][i,k] contains # of responses to \cr
#'               the k-th category (k=0,1,...,ncat[j]) at theta[i] \cr
#'               to which N[i]=round(Ntot*thd[i]) observations belong. \cr
#' \cr
#' \cr
#'            If compress = 1 and Ntot = 1, \cr
#'              U    is   npoints x nitems \cr
#'               U[i,j] contains the response to the j-th item at theta[i] \cr
#'               where 0 <= U[i,j] <= ncat[j]-1. \cr
#'              In this case, fromP and toP are not compressed. \cr
#' \cr
#'            (theta, thd, N)   are   npoints x 1 \cr
#'            N[i]  is the # of observations at theta[i] \cr
#'
#' @examples
#' # 10 observations with theta from N( 0, 1 ): not compressed
#' set.seed(1702)
#' res1 <- gendataIRT( 1, paramS1, npoints=10, thdist="rnorm" )
#' Print(res1$N,res1$U,res1$theta)
#' #
#' # 10 observations with theta from N( 0, 1 ): compressed
#' set.seed(1702)
#' res1 <- gendataIRT( 1, paramS1, npoints=10, thdist="rnorm", compress=1 )
#' Print(res1$N,res1$U,res1$theta)
#' #
#' # 10 observations with theta from N( 0, 1 ): no zero category
#' set.seed(1702)
#' res1 <- gendataIRT( 1, paramS1, npoints=10, thdist="rnorm", zero=0 )
#' Print(res1$N,res1$U,res1$theta)
#' #
#' #
#' #
#' # seven theta points in [-3,3]
#' set.seed(1701)
#' res1 <- gendataIRT( 1, paramS1, npoints=7, thmin=-3, thmax=3 )
#' Print(res1$N,res1$U,res1$theta)
#' #
#' # seven theta points in [-3,3] and compressed
#' set.seed(1701)
#' res1 <- gendataIRT( 1, paramS1, npoints=7, thmin=-3, thmax=3, compress=1 )
#' Print(res1$N,res1$U,res1$theta)
#' #
#' set.seed(1701)
#' # 100 observations tabulated at seven theta points in [-3,3]
#' res1 <- gendataIRT( 1000, paramS1, npoints=7, thmin=-3, thmax=3 )
#' Print(res1$N,res1$U,res1$theta)
#' #
#' #
#' # Using Nmat argument
#' set.seed(1701)
#' npoints=20
#' param=paramS1
#' nitems=nrow(param)
#' Nmat=matrix(100,npoints,nitems)
#' res1 <- gendataIRT( 1, param, Nmat=Nmat, npoints=npoints
#'                   , thdist="rnorm", compress=1 )
#' Print(res1$N,res1$U,res1$theta)
#'
#'
#' @export
#'

gendataIRT <- function( Ntot, param, DinP=1, Nmat=NULL, zero=1
                      , theta=NULL, npoints=31, thmin=-3, thmax=3
                      , thd=NULL, thdist="NORMAL", thmean=0, thstd=1
                      , nomiss=0, compress=0, sort=1, thetaname=NULL ){
 # generation of Item response data
 # Shin-ichi Mayekawa
 # 121017,18,24
 # DinP: 121110(London)
 # renamed as gendataIRT: 121110(London)
 # make the data integer when possible: 121126
 # bugfix: 20150614
 # nomiss: 20150615
 # theta as rownames abandoned: 20161103
 # Nmat: 20161103
 # zero: 20170303
 # bugfix: 20170307
 #
 #
 # Args:
 #
 #   npoints  # of theta points:  or the length of theta
 #
 #   Ntot     # of total observations to be generated
 #            If Ntot == 1, one observation per each theta point,
 #            else the # of observations per each theta point is
 #            proportional to round(Ntot*thd*) where thd is
 #            the npoint x 1 vector or theta distribution.
 #            Therefore, Ntot must be large enough.
 #   param    parameter data frame
 #   theta    theta points
 #            If null, npoints, thmin and thmax will be used to generate it
 #   thd      theta distribution probability vector
 #   thdist   = "NORMAL" or "UNIFORM" or "rnorm"  or  "runif"
 #            When "NORMAL" or "rnorm",
 #                thmean and thstd will be used to generate thd.
 #            When "UNIFORM" or "runif",
 #                thmin and thmax will be used to generate thd.
 #            When "rnorm" or "runif", theta will be generated using
 #                npoints random numbers and thd and Ntot is set equal to 1.
 #   nomiss   = 1 to omit those theta points with N=0.
 #   compress = 1 to compress the data when Ntot = 1.
 #            See the Values section.
 #
 #   zero     = 0 to exclude zero-th category
 #
 #
 # Values:
 #            list of U, fromP, toP, N, theta, thd
 #
 #            If compress = 0 and Ntot > 1,
 #              U    is   npoints x sum of ncat[j]
 #               U[,fromP[j]:toP[j]]  is npoints x ncat[j]
 #               U[,fromP[j]:toP[j]][i,k] contains # of responses to
 #               the k-th category (k=0,1,...,ncat[j]) at theta[i]
 #               to which N[i]=round(Ntot*thd[i]) observations belong.
 #
 #
 #            If compress = 1 and Ntot = 1,
 #              U    is   npoints x nitems
 #               U[i,j] contains the response to the j-th item at theta[i]
 #               where 0 <= U[i,j] <= ncat[j]-1.
 #              fromP and toP are not compressed.
 #
 #            (theta, thd, N)   are   npoints x 1
 #            N[i]  is the # of observations at theta[i]
 #
 #            Response probability out of Ntot can be obtained as
 #              P = U/N
 #
 #            U will be of integer type.
 #

 # const
 nitems=nrow(param)
 type=param$type
 ncat=param$ncat

 # generate theta
 if( is.null(theta) ){
  theta=seq(thmin,thmax,length.out=npoints)
 }
 else{
  npoints=length(as.vector(theta))
 }

 # error check
 if( toupper(thdist) == "NORMAL"  ||  toupper(thdist) == "UNIFORM" )
  if( Ntot > 1  &&  Ntot < 3*npoints ){
   cat("\n\nwarning1(gendataIRT): Ntot may be too small. \n\n")
  }

 if( Ntot > 1  &  compress == 1 ){
  cat("\n\nwarning1(gendataIRT): compress=1 cannot be used when
      Ntot > 1 or Nmat is given. \n\n")
 }

 if( !is.null(Nmat) ){
  Nmatin=1
  if( Ntot != 1 ){
   cat("\n\nerror1(gendataIRT): Ntot must be 1 when Nmat is given.\n\n")
   return()
  }
  if( npoints != nrow(Nmat) || nrow(param) != ncol(Nmat) ){
   cat("\n\nerror1(gendataIRT): Nmat does not comform to param.\n\n")
   return()
  }
 }
 else Nmatin=0

 # theta distribution
 if( is.null(thd) ){
  if( toupper(thdist) == "NORMAL" ){
   thd=exp(-0.5*( (theta-thmean)/thstd )^2); thd=thd/sum(thd)
  }
  else if( toupper(thdist) == "UNIFORM" ){
   thd=1/npoints
  }
  else if( toupper(thdist) == "RNORM" ){
   thd=1; Ntot=1;
   theta=rnorm(npoints,thmean,thstd)
  }
  else if( toupper(thdist) == "RUNIF" ){
   thd=1; Ntot=1;
   theta=runif(npoints,thmin,thmax)
  }
  else{
   thd=1;
   Ntot=1
  }
 }

 if( sort == 1 ) theta=sort(theta)

 # # of obs. per theta point
 if( Ntot > 1 ){
  N=round(Ntot*thd)
  if( sum(N) < Ntot ){
   N[round(npoints/2):(round(npoints/2)-(Ntot-sum(N))+1)]=
                N[round(npoints/2):(round(npoints/2)-(Ntot-sum(N))+1)]+1
  }
 }
 else N=rep(1,npoints)

 if( !( Ntot == 1  && !is.null(Nmat) ) ){
  Nmat=matrix(N,length(N),nitems)
 }

 # generate irf
 res=irf( param, theta, DinP=DinP, zero=1, print=0, plot=0 )
 icrf=res$ICRF; fromP=res$fromP; toP=res$toP; ncat=res$maxscore_i+1
 thname=rownames(icrf)
 rm(res)
 # Print(icrf, fromP, toP, ncat)

 # generate dummy expanded response according to icrf
 Y=matrix(NA,nrow(icrf),ncol(icrf))
 colnames(Y)=colnames(icrf)
 rownames(Y)=1:nrow(Y)
 for( j in 1:nitems ){
  domj=c(0:(ncat[j]-1))
  Pj=icrf[,fromP[j]:toP[j]]
  # Print(ncat[j],Pj,digits=3)
  for( i in 1:npoints ){
   Y[i,fromP[j]:toP[j]]=t( rmultinom( 1, Nmat[i,j], Pj[i,] ) )
  }
 }

 if( nomiss == 1 ){
  # delete rows with N=0
  locobs=which( rowSums(Nmat) > 0 )
  N=N[locobs]; theta=theta[locobs]; thd=thd[locobs]
  icrf=icrf[locobs,]; Y=Y[locobs,]
  npoints=length(locobs)
  Nmat=Nmat[locobs,]
 }

 # compress if requested
 if( Ntot == 1 && compress == 1 && all( Nmat==1 ) ){
  U=matrix(NA,npoints,nitems)
  for( j in 1:nitems ){
   loc=apply( Y[,fromP[j]:toP[j]], 1, function(y) which( y == 1) )
   U[,j]=loc-1
  }
  U=matrix(as.integer(U),nrow(U),ncol(U))
  colnames(U)=param$name
 }
 else U=Y
 rm(Y)

 rownames(U)=1:nrow(U)

 if( !is.null(thetaname) ){
  U=cbind(U, theta)
  colnames(U)[ncol(U)]=thetaname
 }

 # zero
 if( zero == 0 ){
  loc=setdiff(1:ncol(U),grep("_0",colnames(U)))
  U=U[,loc]
 }

 res=list( U=U, N=N, npoints=npoints, theta=theta, thd=thd, icrf=icrf
           , Ntot=Ntot, fromP=fromP, toP=toP, type=type, ncat=ncat
           , thmean=thmean, thstd=thstd
           , compress=compress, sort=sort, zero=zero )
 if( Nmatin ) res=append(res,list(Nmat=Nmat))
 return( res )

} # end of gendataIRT








#' Reading Parameter File
#'
#' @param infile Input file name
#' @param skip # of lies to skip
#' @param nrows NOT used
#' @param na.strings A character indicating missing value
#' @param sep Delimite character
#' @param print = 1 to print result
#'
#' @details
#' It is assumed that the data file has the following quantities
#'   in the following order.\cr
#'\cr
#'   item_name, item_type, ncat, p1, p2, ..., p[ncat+1]\cr
#'\cr
#' The first line must consist of the variable names shown above. \cr
#' If the data file has header lines berore variable names, those lines
#' must be skipped by specifying skip= parameter.
#'\cr
#'   When item_type == "B", ( ncat == 2 )\cr
#'     p1 is the discrimination parameter\cr
#'     p2 is the difficulty parameter\cr
#'     p3 is the gussing parameter\cr
#'\cr
#'   When item_type == "G",\cr
#'     p1 is the discrimination parameter\cr
#'     p2 is the category threshold parameter b_{j1}\cr
#'     p3 is the category threshold parameter b_{j2}\cr
#'\cr
#'     p_ncat is the category threshold parameter b_{j,[ncat-1]}\cr
#'\cr
#'   When item_type == "P",\cr
#'     p1 is the discrimination parameter\cr
#'     p2 is the step parameter b_{j1}\cr
#'     p3 is the step parameter b_{j2}\cr
#'\cr
#'     p_ncat is the step parameter b_{j,[ncat-1]}\cr
#'\cr
#'   When item_type == "PN,\cr
#'     p1 is the slope parameter\cr
#'     p2 is the intercept parameter b_{j1}\cr
#'     p3 is the intercept parameter b_{j2}\cr
#'\cr
#'     p_ncat is the intercept parameter b_{j,[ncat-1]}\cr
#'\cr
#'   When item_type == "N,\cr
#'     p1 - p_(ncat[j]-1) are the slope parameters\cr
#'     p_ncat[j] - p_2*(ncat[j]-1) are the intercept parameters\cr
#'\cr
#'   Regardless of the item_type, there will be ncat[j]\cr
#'   item paramters for item j, except for the Binary Items\cr
#'   which have the gusseing param as the ncat[j]+1st,\cr
#'   and the Nominal Items which have 2*(ncat[j]-1) item parameters.\cr
#'\cr
#' @export
#'
read.param <- function( infile, skip=0, nrows=-1, na.strings="NA"
                        , sep="", print=0 ){
 # reading IRT parameter files
 # Shin-ichi Mayekawa
 # 120208,09,12
 # clean input data frame with NA: 120224
 # use of col.names: 120224
 # a dropped: 120225
 # sep: 120228,29
 # type Bxx items: 120919,21
 # simplified version: 121130
 #
 #
 # It is assumed that the data file has the following quantities
 #   in the following order.
 #
 #   item_name, item_type, ncat, p1, p2, ..., p[ncat+1]
 #
 #
 #   When item_type == "B", ( ncat == 2 )
 #     p1 is the discrimination parameter
 #     p2 is the difficulty parameter
 #     p3 is the gussing parameter
 #
 #   When item_type == "G",
 #     p1 is the discrimination parameter
 #     p2 is the category threshold parameter b_{j1}
 #     p3 is the category threshold parameter b_{j2}
 #
 #     p_ncat is the category threshold parameter b_{j,[ncat-1]}
 #
 #   When item_type == "P",
 #     p1 is the discrimination parameter
 #     p2 is the step parameter b_{j1}
 #     p3 is the step parameter b_{j2}
 #
 #     p_ncat is the step parameter b_{j,[ncat-1]}
 #
 #   When item_type == "PN,
 #     p1 is the slope parameter
 #     p2 is the intercept parameter b_{j1}
 #     p3 is the intercept parameter b_{j2}
 #
 #     p_ncat is the intercept parameter b_{j,[ncat-1]}
 #
 #   When item_type == "N,
 #     p1 - p_(ncat[j]-1) are the slope parameters
 #     p_ncat[j] - p_2*(ncat[j]-1) are the intercept parameters
 #
 #   Regardless of the item_type, there will be ncat[j]
 #   item paramters for item j, except for the Binary Items
 #   which have the gusseing param as the ncat[j]+1st,
 #   and the Nominal Items which have 2*(ncat[j]-1) item parameters.
 #
 #
 # Value: as data.frame
 #
 #   name type ncat  p1  p2  p3   .....
 #

 param=read.table( infile, header=1, fill=1
                 , na.strings=na.strings, skip=skip, flush=1, sep=sep
                 , stringsAsFactors=0 )

 locnotna=which( !apply( param, 2, function(x) all(is.na(x)) ) )
 param=param[,locnotna]
 cname=colnames(param)

 # in case the header has a variable 'name'.
 if( cname[1] == "name" ){
  rownames(param)=param[,1]
  param=param[,-1]
 }

 # param names from the row names
 param=cbind(name=rownames(param),param)
 cname=colnames(param)

 nitems=nrow(param)
 ncol=ncol(param)


 # structural error
 if( !(cname[1] %in% c("name","name")) ){
  cat("\n\n** error1(read.param) Must have a variable named 'name'"
      , " as the first variable.\n\n")
  return()
 }
 if( !(cname[2] == "type") ){
  cat("\n\n** error1(read.param) Must have a variable named 'type'"
      , " as the second variable.\n\n")
  return()
 }
 if( !(cname[3] == "ncat") ){
  cat("\n\n** error1(read.param) Must have a variable named 'ncat'"
      , " as the third variable.\n\n")
  return()
 }

 if( length(grep("^p[[:digit:]]$",cname)) == 0 ){
  cat("\n\n** error1(read.param) Must have a variable named 'p_num_'"
      , " at colums 4, 5, ....\n\n")
  return()
 }
 ncat=param$ncat
 locNO=which( param$type  %in% c("N","NO","O") )
 if( !any(cname %in% paste("p",max(ncat),sep="")) ){
  cat("\n\n** error1(read.param) Too few item paramters given."
      , " Must have 'p1' through ", paste("p",max(ncat),sep=""), ".\n\n")
  return()
 }
 if( length(locNO) > 0 )
  if( !any(cname %in% paste("p",max(ncat[locNO]-1)*2,sep="")) ){
   cat("\n\n** error1(read.param) Too few item paramters given."
       , " Must have 'p1' through "
       , paste("p",max(ncat[locNO]-1)*2,sep=""), ".\n\n")
   return()
  }



 # item names and types
 type=param$type
 locB=grep("^B[[:digit:]]*$", type)

 # use 0 for the c-parameter
 if( length(locB) > 0 ) param[locB,6][is.na(param[locB,6])]=0

 if( print ){
  cat("\n\n The following item parameter data frame was created from ",
      infile,".\n\n")
  Print(param)
 }

 return( param )

} # end of read.param







#' Reading Weight File
#'
#' @param infile Input file name
#' @param skip # of lies to skip
#' @param nrows NOT used
#' @param na.strings A character indicating missing value
#' @param sep Delimite character
#' @param print = 1 to print result
#'
#' @details
#' It is assumed that the data file has the following quantities
#'   in the following order.
#' \cr
#'   item_name item_type  ncat  w  v0  v1  v2  ..., v[ncat-1] \cr
#' \cr
#'  where w is the item weight and v0,v1..., are the category weight. \cr
#'
#'
#'
#' @export
#'
read.weight <- function( infile, skip=0, nrows=-1, na.strings="."
                       , sep="", print=0 ){
 # reading IRT item weight files
 # Shin-ichi Mayekawa
 # 120209,12,25
 # sep: 120229
 #
 #
 # It is assumed that the data file has the following quantities
 #   in the following order.
 #
 #   item_name item_type  ncat  w  v0  v1  v2  ..., v[ncat-1]
 #
 #  where w is the item weight and v0,v1..., are the category weight.
 #
 #
 # Value: as data.frame
 #
 #   name type ncat  w  v0  v1  v2  v3   .....   v[ncat-1]
 #

 weight=read.table( infile, fill=1, blank.lines.skip=1, header=0, nrows=nrows
                    , na.strings=na.strings, skip=skip, flush=1, sep=sep
                    , stringsAsFactors=0
                    , col.names=paste("v",1:20,sep=""))

 nitems=nrow(weight)
 ncol=ncol(weight)

 header=0
 if( header == 0 ){
  rownames(weight)=weight[,1]
  colnames(weight)=c("name", "type", "ncat", "w"
                   , paste("v",0:(ncol(weight)-5),sep=""))
 }

 # remove excessive columns
 ncat=weight[,3]
 # weight[,3:ncol]=as.numeric(weight[,3:ncol])
 if( max(ncat)+4 < ncol ) weight=weight[,-((max(ncat)+5):ncol)]

 # item names and types
 weight[,1]=as.character(weight[,1])
 weight[,2]=toupper(weight[,2])
 dummy=matrix(as.numeric(unlist(weight[,3:ncol(weight)])),nrow(weight))
 weight[,3:ncol(weight)]=dummy


 if( print ){
  cat("\n\n The following item weight data frame was created from ",
      infile,".\n\n")
  Print(weight)
 }

 return( weight )

} # end of read.weight







