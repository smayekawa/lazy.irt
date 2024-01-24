irf <- function( param, theta=NULL, weight=NULL, zero=1, smallP=1e-6
               , thmin=-4, thmax=4, npoints=21, DinP=1
               , print=1, debug=0, plot=0 ){
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
  if( length(grep("^B[[:digit:]]*$", type[i])) > 0 & ncat[i]+loca+1 <= ncp) 
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
 ICRFP=matrix(0,npoints,sum(ncat))
 cnP=character(max(ncat))
 toP=cumsum(ncat)
 fromP=toP-ncat+1
 ll=mapply(seq,fromP,toP)  # simplify removed: 120928
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
 
 if( debug > 0 ) Print(locB,locG,locP,locN,ll,llB,llG,llP,llN)
 
 # for each item type, calculate icrf and store them in ICRFP
 if( nB > 0 ){
  # 3PLM
  pB=as.matrix( param[locB,4:6,drop=0] )
  rownames(pB)=iname[locB]
  colnames(pB)=c("a","b","c")
  iccB=icrfB( pB, theta, zero=1, smallP=smallP )  
  if( debug > 0 ) Print(pB,iccB, digits=2)
  ICRFP[,llB]=iccB;
  cnP[llB]=colnames(iccB) 
 }
 if( nP > 0 ){
  # Generalized Partial Credit with step parameters and DinP
  pP=as.matrix( param[locP,3:ncp,drop=0] )
  rownames(pP)=iname[locP]
  res=icrfP( pP, theta, smallP=smallP, DinP=DinP )
  iccP=res$P  
  if( debug > 0 ) Print(pP,iccP, digits=2)
  ICRFP[,llP]=iccP;
  cnP[llP]=colnames(iccP)  
 }
 if( nPN > 0 ){
  # Generalized Partial Credit in Nominal form wih 1.7
  pPN=as.matrix( param[locPN,3:ncp,drop=0] )
  rownames(pPN)=iname[locPN]
  res=icrfPN( pPN, theta, smallP=smallP )
  iccPN=res$P  
  if( debug > 0 ) Print(pPN,iccPN, digits=2)
  ICRFP[,llPN]=iccPN;
  cnP[llPN]=colnames(iccPN)  
 }
 if( nG > 0 ){
  # Graded Response Model (2PLM)
  pG=as.matrix( param[locG,3:ncp,drop=0] )
  rownames(pG)=iname[locG]
  res=icrfG( pG, theta, smallP=smallP )
  iccG=res$P  
  if( debug > 0 )Print(pG,iccG, digits=2)
  ICRFP[,llG]=iccG;
  cnP[llG]=colnames(iccG)  
 }
 if( nN > 0 ){
  # Nominal Response Model: *** not available yet ***
  pN=as.matrix( param[locN,3:ncp,drop=0] )
  rownames(pN)=iname[locN]
  res=icrfN( pN, theta, smallP=smallP )
  iccN=res$P  
  if( debug > 0 )Print(pN,,iccN, digits=2)
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
 
 if( zero != 1 ){
  # remove the 0-th category
  ICRFP=ICRFP[,-fromP]
  fromP=fromP-0:(nitems-1)
  toP=toP-1:nitems
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
  for( j in 1:nitems ){
   
   # titles
   main=paste("Plot of ICRF of "
              ,iname[j]," ( type = ",type[j],", ncat=",ncat[j]," )")
   sub=paste("param = ", paste(format(param[j,4:(ncat[j]+3)],digits=3)
                               , collapse=",  "))
   if( length(grep("^B[[:digit:]]*$", param$type[j])) > 0 )
    sub=paste("param = ", paste(format(param[j,4:(ncat[j]+4)],digits=3)
                                , collapse=",  "))
   if( param$type[j] == "N" )
    sub=paste("param = ", paste(format(param[j,4:(2*(ncat[j]-1)+4-1)],digits=3)
                                , collapse=",  "))
   
   # set up the plot
   plot(range(theta), c(0,1), type="n", xlab="theta", ylab="icrf" )
   colors <- rainbow(ncat[j])
   linetype <- c(1:ncat[j])
   plotchar <- seq(18,18+ncat[j],1)
   
   # add lines
   for (k in 1:ncat[j]) {
    lines(theta, ICRFP[,fromP[j]+k-1], type="b", lwd=1.5,
          lty=linetype[k], col=colors[k], pch=plotchar[k])
   }
   
   # add a title and subtitle
   title(main, sub)
   
   # add a legend
   legend(range(theta)[1], 0.8, (1:ncat[j])-1, cex=0.8, col=colors,
          pch=plotchar, lty=linetype, title="cat")
  }  
  if( plot >= 2 ){
   for( j in 1:(nitems) ){
    title=paste("Plot of IRF of "
                , iname[j]," ( type = ",type[j],", ncat=",ncat[j]
                , ", score range = [",minscore_i[j],",", maxscore_i[j],"] )")
    sub=paste("param = ", paste(format(param[j,4:(ncat[j]+3)],digits=3)
                                , collapse=",  "), "  (with weights)")
    if( length(grep("^B[[:digit:]]*$", param$type[j])) > 0 )
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
 }
 
 
 return( list(ICRF=ICRFP, IRF=irf, TRF=trf, fromP=fromP, toP=toP
            , vecv=vecv, minscore_i=minscore_i, maxscore_i=maxscore_i
            , minscore_t=minscore_t, maxscore_t=maxscore_t) )
 
 
} # end of irf






temp="
 name type ncat p1 p2 p3 p4
1 QB1 B3 2 1 0 .2 NA
2 QG1 G  3 1 -1 1 NA
3 QP1 P  4 1 -1 0 1
"; paramtest=lazy.matrix(temp,header=1)



npoints=21; minth=-4; maxth=4
theta=seq(minth,maxth,length.out=npoints)

res=irf( paramtest, theta, print=0, debug=0, plot=1)
