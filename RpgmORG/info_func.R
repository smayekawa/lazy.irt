info_func <- function( param, weight=NULL
                     , npoints=31, thmin=-4, thmax=4
                     , print=1, plot=0, debug=0 ){
 # calculation of information function
 # Shin-ichi Mayekawa
 # obscore modified: 20140329
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
 v=weight[,5:ncol(weight)]
 
 maxscore_i=apply(weight[,5:ncol(weight)],1,max,na.rm=1)
 maxscore_t=sum(w*maxscore_i)
 minscore_i=apply(weight[,5:ncol(weight)],1,min,na.rm=1)
 minscore_t=sum(w*minscore_i)
 
 if( debug > 0 ) Print(iname,nitems,nc,ncat,w, "/",maxscore_i, maxscore_t)
 
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
 
 # generate thata and prior theta dist
 theta=seq(thmin,thmax,length.out=npoints)
 thname=format(theta,digits=2)

 
 # calculate icrf and trf=cond. mean of X given theta
 temp=irf( param, theta, weight, print=0, debug=0, plot=0 ) 
 icrf=temp$ICRF
 trf=temp$TRF
 fromP=temp$fromP
 toP=temp$toP
 vecv=temp$vecv
 rm(temp)
 
 if( debug > 0 ) Print(icrf, fromP,toP)
  
 # first derivative of icrf
 temp=dirf( param, theta, weight=weight, print=0 )
 dicrf=temp$dICRF
 slope_trf=temp$dTRF
 dirf=temp$dIRF
 rm(temp)
 
 slopetrf=slope_trf
 
 # conditional variance of the weighted total score X;
 # This should be equal to the one calculated using Px_t.
 stdx_t=matrix(0,npoints,1)
 varu_t=matrix(0,npoints,nitems)
 for( k in 1:npoints ){
  Pk=icrf[k,,drop=0]
  dum=0
  for( j in 1:nitems ){
   vj=matrix(vecv[fromP[j]:toP[j]],ncat[j])*w[j]
   Pjk=Pk[fromP[j]:toP[j]]
   DD=diag(Pjk)-outer(Pjk,Pjk,"*")
   vDDv=t(vj)%*%DD%*%vj
   varu_t[k,j]=vDDv
   dum=dum+vDDv
  }
  stdx_t[k]=dum
 }
 stdx_t=sqrt(stdx_t);
 
 # information function associated with the obaserved weighted score
 info=(slopetrf^2)/(stdx_t^2)
 
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
 info_item_LOW=(dirf^2)/varu_t 
 info_LOW=rowSums(info_item_LOW)
 names(info_LOW)=thname

 if( print > 0 ){
  Print(theta,info,info_LO, info_LOW) 
  Print(info_item_LO)
  Print(info_item_LOW)
 }
 
 if( plot > 0 ){
  matplot( theta,cbind(info, info_LO, info_LOW), type="l"
         , main="Information Functions")
  legend( range(theta)[1], max(info_LO)-.1, c("info","info_LO","info_LOW")
        , cex=0.8, col=c("black", "red", "blue")
        , pch=1:3, lty=1, title="legend" )
 }
 
 res=list( theta=theta, info=info, info_LO=info_LO, info_LOW=info_LOW
         , info_item_LO=info_item_LO, info_item_LOW=info_item_LOW )
 
 
} # end of info_function



resInfo=info_func( param, npoints=15 )


 
 
 
 