irf <- function( param, theta=NULL, weight=NULL, zero=1, smallP=1e-6
               , thmin=-4, thmax=4, npoints=21, DinP=1
               , print=1, debug=0, plot=0, colors="black" ){
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
 
 ncat2=ncat
 if( zero != 1 ){
  # remove the 0-th category
  ICRFP=ICRFP[,-fromP]
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
  
  if( is.null(colors) ) colors <- rainbow(max(ncat2))
  else if( length(colors) == 1 ) colors=rep(colors,max(ncat2))
  
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
    sub=paste("param = ", paste(format(param[j,4:(2*(ncat[j]-1)+3)],digits=3)
                                , collapse=",  "))
   
   # set up the plot
   plot(range(theta), c(0,1), type="n", xlab="theta", ylab="icrf" )
   
   linetype <- c(1:ncat2[j])
   plotchar <- seq(18,18+ncat2[j],1)
   
   # add lines
   for (k in 1:ncat2[j]) {
    #################?????????? if( k == 1  & zero == 0 ) next
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
  
  if( length(a) > 1 ) da=diag(a)
  else da=matrix(a,1,1)
  Z=-1.7*outer(theta,b,"-")%*%da
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
  rownames(P)=format(theta,digits=3)
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
   Print(P, digits=3)
  }
  if( plot > 0 ){
   matplot(theta,P, type = "l", ylab="IRF"
         , xlim=c(min(theta),max(theta)), ylim=c(0,1)
         , main="Plot of ICRF of Binary Items")
  }
  return( P )
  
} # end of icrfB



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
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0,ncat[j]-1} Ez_{jk}(theta)
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
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0,ncat[j]-1} Ez_{jk}(theta)
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


icrfN <- function( param, theta, smallP=0, print=0, plot=0, debug=0, DinP=1 ){
 # calculation of Nominal Response Model ICRF 
 # Shin-ichi Mayekawa
 # 121123
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
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0,ncat[j]-1} Ez_{jk}(theta)
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
  cj=c[j,1:ncat1[j],drop=0]
  beta=c(0,a[j,1:ncat1[j]])
  alpha=cbind(0 , cj)

  if( debug ) Print(cj,beta,alpha)
  # z=theta*beta+j(npoints,1,1)*alpha
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
   sub=paste("param = ", paste(param[j,1:(2*ncat1[j])+1], collapse=",  "))
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
 
} # end of icrfN







dirf <- function( param, theta=NULL, weight=NULL, zero=1, smallP=0
                  , thmin=-4, thmax=4, npoints=21, DinP=1
                  , print=1, debug=0, plot=0, colors="black" ){
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
 # colors added: 20150315
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
 
 # for each item type, calculate icrf and store them in icrfP
 if( nB > 0 ){
  pB=as.matrix( param[locB,4:6,drop=0] )
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
  
  if( is.null(colors) ) colors <- rainbow(max(ncat))
  else if( length(colors) == 1 ) colors=rep(colors,max(ncat))
  
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






dicrfB <- function( param, theta, maxZ=700, zero=0, smallP=0
                  , print=0, plot=0 ){
 # derivative of Binary Logistic ICRF
 # Shin-ichi Mayekawa
 # 120214,15
 # bugfix: 120217
 # name as factor: 120224
 # checkparam: 120224,29
 # smallP: 121022
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
  rownames(dP)=format(theta,digits=3)
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
   Print(dP, digits=3)
  }
  if( plot > 0 ){
   matplot(theta,dP, type = "l", ylab="IRF"
         , xlim=c(min(theta),max(theta)), ylim=
          , main="Plot of the Derivatives of ICRF of Binary Items")
  }
  return( dP )
  
} # end of dicrfB



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
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0,ncat[j]-1} Ez_{jk}(theta)
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


dicrfP <- function( paramP, theta, DinP=1, smallP=0, print=0, plot=0, debug=0 ){
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
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0,ncat[j]-1} Ez_{jk}(theta)
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
  paramj0=paramj
  paramj0[,4:(4+length(paramnum)-1)]=paramnum
  res=irf( paramj0, theta, zero=0, print=0 )$ICRF
  res=matrix(res,,1)
  # vertically stack all categories
  return(  matrix( res, , 1 )  )
 }
 
 
 
 
 # const
 ncatj=paramj$ncat
 ncatj1=ncatj-1
 typej=paramj$type
 
 # # of parameters to be estimated
 if( typej == "B3" ) np=ncatj+1
 else if( typej == "N" ) np=2*ncatj1
 else np=ncatj
 
 if( is.null(theta) ){
  theta=seq(thmin,thmax,length.out=npoints)
 }
 else{
  if( is.matrix(theta) ) theta=as.vector(theta)
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
 if( typej == "G" ) b=cbind(b,9999)
 c=paramj$p3
 
 # storage
 Jack=matrix(0,ncatj1*lenth,np)
 colnames(Jack)=paste("p",1:np,sep="")
 
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
   temp=1.7*
    ( (theta-b[k])*PPj[,k]*(1-PPj[,k])-(theta-b[k+1])*PPj[,k+1]*(1-PPj[,k+1]) )
   Jack[((k-1)*lenth+1):(k*lenth),1]=temp
   for( z in 1:ncatj1 ){
    temp=0
    if( k == z ) temp=-1.7*a*PPj[,z]*(1-PPj[,z])
    else if( (k+1) == z ) temp=1.7*a*PPj[,z]*(1-PPj[,z])
    Jack[((k-1)*lenth+1):(k*lenth),1+z]=temp
   } # z
  } # k       
  
 } # end of type = "G"
 else if( typej == "PN" ){
 
  if( is.null(Pj) ){
   # icrf : no zero-th category
   Pj=icrfPN( paramj, theta, smallP=smallP )$P[,2:ncatj,drop=0]
  }
  # Jacobian
  thetab=matrix(theta,lenth,ncatj1)-matrix(1,lenth,1)%*%b
  Pqthetab=rowSums( (Pj%*%diag(1:ncatj1))*thetab )
  Jack=matrix(0,ncatj1*lenth,np)
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
  Jack=matrix(0,ncatj1*lenth,np)
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
  
#Print(Jack,digits=2)
#Jack=JacobianMat( c(a,c), icrfP00, eps = 1e-06 )
#Print(Jack,digits=2)    
  
 } # end of type = "N"
 
 
 
 if( zero == 1 ){
  # dP0dz = d(1-sum_k Pj)dz = 0 - sum_k dPjdz
  temp=0
  for( k in 1:ncatj1 ){
   temp=temp+Jack[((k-1)*lenth+1):(k*lenth),]
  } # k    
  Jack=rbind(-temp, Jack)
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
  Print(Pj,Jack, digits=2)
 }
 
 return( list(Jack=Jack, Pj=Pj, PPj=PPj) )
 
} # end of dirf_p













convP2PN <- function( param, print=0, DinP=1 ){
 # convert the original parametrization ( no 1.7 and step parameters)
 # to x-intercept parameters with 1.7
 # Shin-ichi Mayekawa
 # 20121015
 # DinP option:121109(London)
 # 121113
 # 
 #    ICRF of category k of item j at theta, namely, p_{jk}(theta)
 #    is defined as
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0,ncat[j]-1} Ez_{jk}(theta)
 #      where 
 #
 #    if type = "PN",
 #     Ez_{jk}(theta) = exp( 1.7 a_j k (theta - b_{jk}) )
 #                                              , k=0,1, ..., ncat[j]-1
 #    with b_{j0}=0.
 #
 #    or
 # 
 #    if type = "P",
 #     Ez_{jk}(theta) = exp( 1.7^DinP a*_j k (theta - sum_{h=0}^k b*_{jh}) )
 #                    = exp( 1.7^DinP a*_j k (theta - sum_{h=1}^k b*_{jh}) )
 #                                              , k=0,1, ..., ncat[j]-1
 #    with b*_{j0}=0.
 #
 #    Note that if DinP=1, 1.7 will be used.
 #
 #    This b*_{jk} is the original step parameter and it is 
 #    the value of theta where P_{jk-1} and P_{jk} intersect.
 #
 #
 locP=which( param$type == "P" )
 if( length(locP) == 0 )
  return( param )
 
 param0=param[locP,,drop=0]
 ncat=param0$ncat
 locparam=grep("^p[[:digit:]]$",colnames(param))
 if( length(locparam) >= 1 ) locparam=locparam[1]
 else locparam=grep("^a$",colnames(param))
 bs=param0[,(locparam+1):(locparam+max(ncat)-1),drop=0]
 b=bs
 b=t( apply( bs,1,cumsum)*(1/(1:ncol(bs))) )
 if( DinP == 0 ) param0[,locparam]=param0[,locparam]/1.7
 param0[,(locparam+1):(locparam+max(ncat)-1)]=b
 param0$type="PN"
 param[locP,]=param0
 
 if( print > 0 ){
  cat("\n item parameters with D in P =", DinP, "\n")
  Print(param0)
 }
 
 return( param )
 
} # end of convP2PN



convPN2P <- function( param, print=0, DinP=1 ){
 # convert the x-intercept parametrization ( with 1.7 and b-type parameters)
 # to the original step parameters w/o 1.7
 # Shin-ichi Mayekawa
 # 20121015
 # DinP option:121109(London)
 # 
 #    ICRF of category k of item j at theta, namely, p_{jk}(theta)
 #    is defined as
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0,ncat[j]-1} Ez_{jk}(theta)
 #      where 
 #
 #    if type = "PN",
 #     Ez_{jk}(theta) = exp( 1.7 a_j k (theta - b_{jk}) )
 #                                              , k=0,1, ..., ncat[j]-1
 #    with b_{j0}=0.
 #
 #    or
 # 
 #    if type = "P",
 #     Ez_{jk}(theta) = exp( 1.7^DinP a*_j k (theta - sum_{h=0}^k b*_{jh}) )
 #                    = exp( 1.7^DinP a*_j k (theta - sum_{h=1}^k b*_{jh}) )
 #                                              , k=0,1, ..., ncat[j]-1
 #    with b*_{j0}=0.
 #
 #    Note that if DinP=1, 1.7 will be used.
 #
 #    This b*_{jk} is the original step parameter and it is 
 #    the value of theta where P_{jk-1} and P_{jk} intersect.
 #
 #
 
 locPN=which( param$type == "PN" )
 if( length(locPN) == 0 )
  return( param )
 
 param0=param[locPN,,drop=0]
 ncat=param0$ncat
 locparam=grep("^p[[:digit:]]$",colnames(param))
 if( length(locparam) >= 1 ) locparam=locparam[1]
 else locparam=grep("^a$",colnames(param))
 b=param0[,(locparam+1):(locparam+max(ncat)-1)]
 bs=b
 for( j in 1:nrow(param0) ){
  for( k in 2:(ncat[j]-1) ){
   bs[j,k]=k*b[j,k]-(k-1)*b[j,k-1]
  }
 }
 
 if( DinP == 0 ) param0[,locparam]=1.7*param0[,locparam]
 param0[,(locparam+1):(locparam+max(ncat)-1)]=bs
 param0$type="P"
 param[locPN,]=param0
 
 if( print > 0 ){
  cat("\n item parameters with D in P =", DinP, "\n")
  Print(param0)
 }
 
 return( param )
 
} # end of convPN2P



convP2PN0 <- function( param, print=0, DinP=1 ){
 # convert the original parametrization ( no 1.7 and step parameters)
 # to a_j theta + c_jk type w/o 1.7
 # Shin-ichi Mayekawa
 # 20121122
 # 
 #    ICRF of category k of item j at theta, namely, p_{jk}(theta)
 #    is defined as
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0,ncat[j]-1} Ez_{jk}(theta)
 #      where 
 #
 #    if type = "PN",
 #     Ez_{jk}(theta) = exp( 1.7 a_j k (theta - b_{jk}) )
 #                                              , k=0,1, ..., ncat[j]-1
 #    with b_{j0}=0.
 #
 #    or
 # 
 #    if type = "P",
 #     Ez_{jk}(theta) = exp( 1.7^DinP a*_j k (theta - sum_{h=0}^k b*_{jh}) )
 #                    = exp( 1.7^DinP a*_j k (theta - sum_{h=1}^k b*_{jh}) )
 #                                              , k=0,1, ..., ncat[j]-1
 #    with b*_{j0}=0.
 #
 #    Note that if DinP=1, 1.7 will be used.
 #
 #    This b*_{jk} is the original step parameter and it is 
 #    the value of theta where P_{jk-1} and P_{jk} intersect.
 #
 #
 locP=which( param$type == "P" )
 if( length(locP) == 0 )
  return( param )
 
 param0=param[locP,,drop=0]
 ncat=param0$ncat
 locparam=grep("^p[[:digit:]]$",colnames(param))
 if( length(locparam) >= 1 ) locparam=locparam[1]
 else locparam=grep("^a$",colnames(param))
 bs=param0[,(locparam+1):(locparam+max(ncat)-1),drop=0]
 
 c=-t( apply( bs,1,cumsum) )*param0[,locparam]
 param0[,(locparam+1):(locparam+max(ncat)-1)]=c
 if( DinP == 1 ) 
  param0[,locparam:(locparam+max(ncat)-1)]=
  param0[,locparam:(locparam+max(ncat)-1)]*1.7
 param0$type="PN0"
 param0$type="PN"
 param[locP,]=param0
 
 if( print > 0 ){
  cat("\n item parameters with D in P =", DinP, "\n")
  Print(param0)
 }
 
 return( param )
 
} # end of convP2PN0



convP2N <- function( param, print=0, DinP=1 ){
 # convert the original parametrization ( no 1.7 and step parameters)
 # to a_jk theta + c_jk type nominal response model w/o 1.7
 # where a_jk = a_k * k
 # Shin-ichi Mayekawa
 # 20121122,1202
 # 
 #    ICRF of category k of item j at theta, namely, p_{jk}(theta)
 #    is defined as
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0,ncat[j]-1} Ez_{jk}(theta)
 #      where 
 #
 # 
 #    if type = "P",
 #     Ez_{jk}(theta) = exp( 1.7^DinP a*_j k (theta - sum_{h=0}^k b*_{jh}) )
 #                    = exp( 1.7^DinP a*_j k (theta - sum_{h=1}^k b*_{jh}) )
 #                                              , k=0,1, ..., ncat[j]-1
 #    with b*_{j0}=0.
 #
 #    Note that if DinP=1, 1.7 will be used.
 #
 #    This b*_{jk} is the original step parameter and it is 
 #    the value of theta where P_{jk-1} and P_{jk} intersect.
 #
 #
 
 locP=which( param$type == "P" )
 if( length(locP) == 0 )
  return( param )
 
 nitems=nrow(param)
 param0=param[locP,,drop=0]
 ncat=param0$ncat
 maxncat=max(ncat)
 locparam=grep("^p[[:digit:]]$",colnames(param))
 if( length(locparam) >= 1 ) locparam=locparam[1]
 else locparam=grep("^a$",colnames(param))
 bs=param0[,(locparam+1):(locparam+max(ncat)-1),drop=0]
 
 c=-t( apply( bs,1,cumsum) )*param0[,locparam]
 param0[,(locparam+1):(locparam+max(ncat)-1)]=c
 if( DinP == 1 ) 
  param0[,locparam:(locparam+max(ncat)-1)]=
  param0[,locparam:(locparam+max(ncat)-1)]*1.7
 a=param0[,locparam]
 param0$type="PN0"
 param[locP,]=param0
 
 np=max(2*(maxncat-1),max(param$ncat))
 param0=data.frame( param[,1:3]
                    ,matrix(NA,nitems,np,dimnames=list(NULL,paste("p",1:np,sep="")))
 )
 param0[,1:ncol(param)]=param
 for( j in 1:length(locP) ){
  param0[locP[j],4:(4+ncat[j]-2)]=a[j]*(1:(ncat[j]-1))
  param0[locP[j],(4+ncat[j]-1):(4+2*(ncat[j]-1)-1)]=
   param[locP[j],5:(5+ncat[j]-2)]
 }
 param0$type[locP]="N"
 param=param0
 
 if( print > 0 ){
  cat("\n item parameters with D in P =", DinP, "\n")
  Print(param0)
 }
 
 return( param )
 
} # end of convP2N






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
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0,ncat[j]-1} Ez_{jk}(theta)
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
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0,ncat[j]-1} Ez_{jk}(theta)
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




fitG2P <- function( paramP, wtype=1, DinP=1, dataframe=1
                  , npoints=21, thmin=-3, thmax=3
                  , maxabsv=1, print=1, plot=0, debug=0 ){
 # fitting Graded Response Model to Generalized Partial Credit Model
 # Shin-ichi Mayekawa
 # iml version: 080127 -- 080312
 # 120202
 # param data frame: 120214
 # PN renamed as PN0   b-type paramters: 120301
 # maxth -> thmax etc: 120307
 # param values in plot: 120309,12c
 # data frame output and DinP: 121113,14
 #
 #
 # Args:
 #
 #     paramP   parameter data frame for type = "P" or "PN" items with DinP
 #
 #     DinP = 1 to use 1.7 in Exp of partial credit model
 #
 #     npoints  # of theta points in [thmin, thmax]
 #     wtype = 1 to use normal weight, else uniform
 #
 #     dataframe = 1 to create parameter data frame, not matrix
 #
 #
 # Value
 #     paramG  = graded response model item parameter matrix
 #
 #
 # Needs:
 #   icrfG, icrfPN, convP2PN
 #
 
 # convert P to PN
 paramPN=convP2PN( paramP, DinP=1 )

 # argument name
 paramname=as.character(substitute(paramPN))
 isdf=0
 if( is.data.frame(paramPN) ){
  # param and weight given as data frame
  isdf=1
  paramPN=checkparam( paramPN, "PN", "fitG2P" )
  if( is.null(paramPN) ){
   cat("\n\nInput Item Parameter ERROR.\n\n")
   return()
  }
 }
 
 # const
 theta=seq(thmin,thmax,length.out=npoints)
 thetac=format(theta,digits=3)
 
 # weights for those elements with which abs(v) > maxabsv will be
 # set equal to min(W)

 nitems=nrow(paramPN)
 iname=rownames(paramPN)
 ncat=paramPN[,1]
 ncat1=ncat-1
 
 if( print >= 1 ){
  cat("\nFitting Graded Response Model to Generalized Partial Credit Model\n")
  cat("  # of items =", nitems, ", weight type =", wtype,"\n")
  cat("  # of theta points =", npoints, " in ["
       , thmin, " , ", thmax, "]\n")
  cat("  D in Partial Credit Model =", DinP,"\n")
  cat("  max absolute value of v =", maxabsv ,"\n")
 }
 
 # weight
 if( wtype == 0 ) W=matrix(0,npoints)
 else if( wtype == 1 ){
  W=exp(-0.5*theta*theta); W=W/sum(W)
 }
 Wmin=min(W); Wmax=max(W);

 # Gen. Partial Credit PP matrix
 res=icrfPN( paramPN, theta )
 PP=res$PP; P=res$P
 toP=res$toP; fromP=res$fromP
 toPP=res$toPP; fromPP=res$fromPP
 rm(res)
 PPname=colnames(PP)
 Pname=colnames(P)

 # storage
 paramG=matrix(0,nitems,max(ncat)+1
     , dimnames=list(iname,c("ncat","a",paste("b",1:(max(ncat)-1),sep=""))))
 paramG[,1]=ncat
 
 # for each item
 for( j in 1:nitems ){
  
  # generate design matrix U
  U=matrix(0,(ncat1[j])*npoints, ncat[j])
  U[,1]=rep(theta,ncat1[j])
  for( k in 1:(ncat1[j]) ){
   U[((k-1)*npoints+1):(k*npoints),k+1]=rep(1,npoints)
  }
  v=matrix( log(PP[,(fromPP[j]+1):(toPP[j]-1)]) -
                  log(1-PP[,(fromPP[j]+1):(toPP[j]-1)],,1) )
  # adjust weight according to the magnitude of logit
  Wj=rep(W,ncat1[j])
  Wj[which(v < -maxabsv)]=Wmin; Wj[which(v > maxabsv)]=Wmin
  Wj=Diag(Wj)
  if( debug ) Print(j, U,v, vecdiag(Wj), digits=3)
 
  # solve
  beta=(1/1.7) * solve(t(U)%*%Wj%*%U)%*%t(U)%*%Wj%*%v
  if( debug ) Print(beta)
  paramG[j,2]=beta[1]
  paramG[j,3:(3+ncat1[j]-1)]=-beta[2:(2+ncat1[j]-1)]/beta[1]
  
 } # end of j loop
 
 # truncate
 paramG[,2:ncol(paramG)]=trunc(paramG[,2:ncol(paramG)]*10000)/10000
 
 if( debug ) Print( paramPN, paramG )
 PG=icrfG( paramG, theta )$P
 diffP=matrix(colMeans(abs(P-PG)),,1,dimnames=list(colnames(PG),""))
 if( print ){
  cat("\nInput and Converted Item Parameters\n")
  cat("Input Partial Credit Item Parameters\n")
  print(paramP)
  cat("\nEstimated Graded Response Model Item Parameters\n")
  print(paramG)
  cat("\n Mean Absolute Differences between two models\n")
  print(diffP)
  if( print >= 2 ){
   cat(" ICRFs for the original Partial Credit Model\n")
   Print(P, digits=3)
   cat(" ICRFs for the estimated Graded Response Model\n")
   Print(PG, digits=3)
  }  
 }
 
 if( plot ){
  for( j in 1:(nitems) ){
   title=paste("Plot of ICRF of ",iname[j]," (ncat=",ncat[j],") [fit G to P]")
   sub=paste( "paramPN ="
         , paste(format(paramPN[j,2:(ncat[j]+1)],digits=3) , collapse=",  ")
         , "\nparamG ="
         , paste(format(paramG[j,2:(ncat[j]+1)],digits=3) , collapse=",  ")
             )   
   for( k in 1:(ncat[j]) ){
    plot(theta,P[,fromP[j]+k-1], xlab=""
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="")
    par(new=1)
   }
   for( k in 1:(ncat[j]-1) ){
    plot(theta,PG[,fromP[j]+k-1], xlab=""
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="b", ylab="")
    par(new=1)
   }
   plot(theta,PG[,toP[j]], main=title, sub=sub, xlab=""
        , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="b", ylab="ICRF")
   par(new=0)
  }
  if( plot >= 2 ){
   for( j in 1:(nitems) ){
    difP=P-PG
    maxd=max(difP); mind=min(difP)
    title=paste("Plot of ICRF Difference of ",iname[j]," (ncat=",ncat[j]
                ,") [fit G to P]")
    for( k in 1:(ncat1[j]-1) ){
     plot(theta,difP[,fromP[j]+k-1]
          , xlim=c(min(theta),max(theta)), ylim=c(mind,maxd)
          , type="l", ylab="")
     par(new=1)
    }
    plot(theta,difP[,toP[j]], main=title, sub="Partial - Graded(fitted)"
         , xlim=c(min(theta),max(theta)), ylim=c(mind,maxd), type="l"
         , ylab="difference")
    par(new=0)
   }
  }
 }
 
 # convert the result to data frame
 if( dataframe == 1 ){
  paramG=data.frame( name=paste("Q",1:nitems, sep="")
                   , type=rep("G",nitems), paramG, stringsAsFactors=0 )
 }
 
 
 return( paramG )
 
} # end of fitG2P




fitP2G <- function( paramG, wtype=1, DinP=1, dataframe=1
                  , npoints=21, thmin=-3, thmax=3
                  , maxabsv=3, print=1, plot=0, debug=0 ){
 # fitting  Generalized Partial Credit Model to Graded Response Model
 # Shin-ichi Mayekawa
 # iml version: 080127 -- 080312
 # 120202
 # param data frame: 120214
 # checkparam: 120229
 # PN renamed as PN0   b-type paramters: 120301,04
 # b-type parameters corrected: 120307
 # param values in plot: 120309,12c
 # data frame output and DinP: 121113,14
 #
 # Args:
 #
 #     paramG   parameter data frame for type = "G" items
 #     DinP = 1 to use 1.7 in Exp of partial credit model
 #
 #     npoints  # of theta points in [thmin, thmax]
 #     wtype = 1 to use normal weight, else uniform
 #
 #     dataframe = 1 to create parameter data frame, not matrix
 #
 # Value
 #   list of two types of item paramter matrices
 #     paramP   =  priginal parameter set with DinP
 #     paramPN  =  b type parameters
 #     paramPN0 =  c type paramters
 #
 #
 # Needs:
 #   icrfG, icrfPN
 #
 
 # argument name
 paramname=as.character(substitute(paramG))
 isdf=0
 if( is.data.frame(paramG) ){
  # param and weight given as data frame
  isdf=1
  paramG=checkparam( paramG, "G", "fitP2G" )
  if( is.null(paramG) ){
   cat("\n\nInput Item Parameter ERROR.\n\n")
   return()
  }
 }

 # const
 iname=rownames(paramG)
  
 # const
 theta=seq(thmin,thmax,length.out=npoints)
 thetac=format(theta,digits=3)
 
 # weights for those elements with which abs(v) > maxabsv will be
 # set equal to min(W)
 

 nitems=nrow(paramG)
 iname=rownames(paramG)
 ncat=paramG[,1]
 ncat1=ncat-1
 
 if( print >= 1 ){
  cat("\nFitting Generalized Partial Credit Model to Graded Response Model\n")
  cat("  # of items =", nitems, ", weight type =", wtype,"\n")
  cat("  # of theta points =", npoints, " in ["
       , thmin, " , ", thmax, "]\n")
  cat("  D in Partial Credit Model =", DinP,"\n")
  cat("  max absolute value of v =", maxabsv ,"\n")
 }
 if( debug ) Print(paramG, ncat1)
 
 # weight
 if( wtype == 0 ) W=matrix(1,npoints)
 else if( wtype == 1 ){
  W=exp(-0.5*theta*theta); W=W/sum(W)
 }
 Wmin=min(W); Wmax=max(W);

 # Graded P matrix
 res=icrfG( paramG, theta )
 P=res$P
 toP=res$toP; fromP=res$fromP
 rm(res)
 Pname=colnames(P)
 
 # storage
 paramPN0=matrix(NA,nitems,max(ncat)+1
             , dimnames=
                list(iname,c("ncat","a",paste("c",1:(max(ncat)-1),sep=""))))
 paramPN0[,1]=ncat
 
 # for each item
 for( j in 1:nitems ){
  
  # generate design matrix U
  U=matrix(0,ncat1[j]*npoints, ncat[j])
  U[,1]=(1:ncat1[j])%x%theta
  for( k in 1:ncat1[j] ){
   U[((k-1)*npoints+1):(k*npoints),k+1]=rep(1,npoints)
  }
  v=matrix( log( P[,(fromP[j]+1):toP[j]] / P[,fromP[j]] ) ,,1) 
  # adjust weight according to the magnitude of logit
  Wj=rep(W,ncat1[j])
  Wj[which(v < -maxabsv)]=Wmin; Wj[which(v > maxabsv)]=Wmin
  Wj=Diag(Wj)
  if( debug ) Print(j, U,v, vecdiag(Wj), digits=3)
 
  # solve
  beta=solve(t(U)%*%Wj%*%U)%*%t(U)%*%Wj%*%v
  if( debug ) Print(beta)
  paramPN0[j,2]=beta[1]
  paramPN0[j,3:(3+ncat1[j]-1)]=beta[2:(2+ncat1[j]-1)]
  
 } # end of j loop
 
 # truncate
 paramPN0[,2:ncol(paramPN0)]=trunc(paramPN0[,2:ncol(paramPN0)]*10000)/10000
 
 # c to b conversion (120307)
 paramPN=paramPN0
 paramPN[,2]=paramPN[,2]/1.7
 colnames(paramPN)=c("ncat","a",paste("b",1:(max(ncat)-1),sep=""))
 for( k in 3:ncol(paramPN) ){
  paramPN[,k]=-paramPN[,k]/( 1.7*(k-2)*paramPN[,2] )
 }
 # truncate
 paramPN[,2:ncol(paramPN)]=trunc(paramPN[,2:ncol(paramPN)]*10000)/10000
 
 # original set of paramters
 paramP=paramPN
 paramP[,2]=paramPN0[,2]
 if( DinP == 0 ) paramP[,2]=paramP[,2]/1.7
 b=paramPN[,3:ncol(paramPN)]
 bs=b
 for( j in 1:nrow(paramPN) ){
  for( k in 2:(ncat[j]-1) ){
   bs[j,k]=k*b[j,k]-(k-1)*b[j,k-1]
  }
 }
 paramP[,3:ncol(paramP)]=bs
 
 
 PP=icrfPN0( paramPN0, theta )$P
 diffP=matrix(colMeans(abs(P-PP)),,1,dimnames=list(colnames(PP),""))
 if( print ){
  cat("\nInput and  Converted Item Parameters\n")
  cat("\nInput and  Graded Response Item Parameters\n")
  Print(paramG, digits=3)
  cat("\n Converted c-type parameters\n")
  Print(paramPN0, digits=3)
  cat("\n  Converted b-type parameters\n")
  Print(paramPN, digits=3)
  cat("\n  Converted parameters in original format with D=1.7\n")
  Print(paramP, digits=3)
  cat("\n Mean Absolute Differences between two models\n")
  print(diffP, digits=3)
  if( print >= 2 ){
   cat(" ICRFs for the original Graded response Model\n")
   Print(P, digits=3)
   cat(" ICRFs for the estimated Partial Credit Model\n")
   Print(PP, digits=3)
  }  
 }
 
 if( plot ){
  for( j in 1:(nitems) ){
   title=paste("Plot of ICRF of ",iname[j]," (ncat=",ncat[j],") [fit P to G]")
   sub=paste( "paramG ="
           , paste(format(paramG[j,2:(ncat[j]+1)],digits=3) , collapse=",  ")
           , "\nparamPN ="
           , paste(format(paramPN[j,2:(ncat[j]+1)],digits=3) , collapse=",  ")
              )   
   for( k in 1:(ncat[j]) ){
    plot(theta,P[,fromP[j]+k-1], xlab=""
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="")
    par(new=1)
   }
   for( k in 1:(ncat[j]-1) ){
    plot(theta,PP[,fromP[j]+k-1], xlab=""
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="b", ylab="")
    par(new=1)
   }
   plot(theta,PP[,toP[j]], main=title, sub=sub, xlab=""
        , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="b", ylab="ICRF")
   par(new=0)
  }
  
  if( plot >= 2 ){
   for( j in 1:(nitems) ){
    difP=P-PP
    maxd=max(difP); mind=min(difP)
    title=paste("Plot of ICRF Difference of ",iname[j]," (ncat=",ncat[j]
                ,") [fit P to G]")
    for( k in 1:(ncat1[j]-1) ){
     plot(theta,difP[,fromP[j]+k-1]
          , xlim=c(min(theta),max(theta)), ylim=c(mind,maxd)
          , type="l", ylab="")
     par(new=1)
    }
    plot(theta,difP[,toP[j]], main=title, sub="Partial - Graded(fitted)"
         , xlim=c(min(theta),max(theta)), ylim=c(mind,maxd), type="l"
         , ylab="difference")
    par(new=0)
   }
  }
 }
 
 # convert the result to data frame
 if( dataframe == 1 ){
  colnames(paramPN)=c("ncat",paste("p",1:max(ncat),sep=""))
  colnames(paramP)=colnames(paramPN)
  colnames(paramPN0)=colnames(paramPN)
  paramPN=data.frame( name=paste("Q",1:nitems, sep="")
                    , type=rep("PN",nitems), paramPN, stringsAsFactors=0 )
  paramPN0=data.frame( name=paste("Q",1:nitems, sep="")
                     , type=rep("PN0",nitems), paramPN0, stringsAsFactors=0 )
  paramP=data.frame( name=paste("Q",1:nitems, sep="")
                   , type=rep("P",nitems), paramP, stringsAsFactors=0 )
 }
  
 return( list(paramPN=paramPN, paramPN0=paramPN0, paramP=paramP) )
 
} # end of fitP2G







gendataIRT <- function( Ntot, param, DinP=1
                        , theta=NULL, npoints=31, thmin=-3, thmax=3
                        , thd=NULL, thdist="NORMAL", thmean=0, thstd=1
                        , compress=0, sort=1 ){
 # generation of Item response data
 # Shin-ichi Mayekawa
 # 121017,18,24
 # DinP: 121110(London)
 # renamed as gendataIRT: 121110(London)
 # make the data integer when possible: 121126
 #
 #  
 # Args:
 #
 #   npoints  # of theta points:  or the length of theta
 #
 #   Ntot     # of total observations to be generated
 #            If Ntot == 1, one observation per each theta point
 #            else the # of observations per each theta point is
 #            proportional to round(Ntot*thd) where thd is
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
 #                npoints random numbers and thd nad Ntotis set equal to 1.
 #   compress = 1 to compress the data when Ntot = 1.
 #            See the Values section.
 #
 #
 # Values:
 #            list of U, fromP, toP, N, theta, thd
 #
 #            If compress = 0 and Ntot is large,
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
 #              In this case, fromP and toP are not compressed.
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
 
 # generate theta
 if( is.null(theta) ){
  theta=seq(thmin,thmax,length.out=npoints)
 }
 else{
  npoints=length(as.vector(theta))
 }
 
 # error check
 if( toupper(thdist) == "NORMAL"  ||  toupper(thdist) == "UNIFORM" )
  if( Ntot < 5*npoints ){
   cat("\n\nwarning(gendataIRT): Ntot may be too small. \n\n")
  }
 
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
  N[round(npoints/2):(round(npoints/2)-(Ntot-sum(N))+1)]=
   N[round(npoints/2):(round(npoints/2)-(Ntot-sum(N))+1)]+1
 }
 else N=rep(1,npoints)
 
 # generate irf
 res=irf( param, theta, DinP=DinP, zero=1, print=0, plot=0 )
 icrf=res$ICRF; fromP=res$fromP; toP=res$toP; ncat=res$maxscore_i+1
 rm(res)
 # Print(icrf, fromP, toP, ncat)
 
 # generate dummy expanded response according to icrf
 Y=matrix(NA,nrow(icrf),ncol(icrf))
 dimnames(Y)=dimnames(icrf)
 for( j in 1:nitems ){
  domj=c(0:(ncat[j]-1))
  Pj=icrf[,fromP[j]:toP[j]]
  # Print(ncat[j],Pj,digits=3)
  for( i in 1:npoints ){
   Y[i,fromP[j]:toP[j]]=t( rmultinom( 1, N[i], Pj[i,] ) )
  }
 }
 
 # compress if requested
 if( Ntot == 1 && compress == 1 ){
  U=matrix(NA,npoints,nitems)
  for( j in 1:nitems ){
   loc=apply( Y[,fromP[j]:toP[j]], 1, function(y) which( y == 1) )
   U[,j]=loc-1
  }
  Y=matrix(as.integer(U),nrow(U),ncol(U))
  rownames(Y)=rownames(icrf)
  colnames(Y)=param$name
 }
 
 return( list( U=Y, fromP=fromP, toP=toP
               , Ntot=Ntot, N=N, theta=theta, thd=thd ) )
 
} # end of gendataIRT




dummy_expand <- function( Uc, ncat=NULL, type=NULL, zero=1 ){
 # dummy expand the categorical valued variable u
 # Shin-ichi Mayekawa
 # 20121022
 # matrix input: 121024
 # missing value: 121025
 # make the data integer when possible: 121126
 
 nitems=ncol(Uc)
 npoints=nrow(Uc)
 iname=colnames(Uc)
 
 if( is.null(ncat) ){
  ncat=matrix(apply(Uc,2,max,na.rm=1),,1)+1
 }
 toP=cumsum(ncat)
 fromP=toP-ncat+1
 # Print(nitems,npoints,ncat,toP,fromP, iname)
 U=matrix(as.integer(0),npoints,toP[nitems])
 for( j in 1:nitems ){
  uj=Uc[,j]
  locobs=!is.na(uj)
  min=0; max=ncat[j]-1
  factuj=factor(uj, levels=min:max)
  Uj=model.matrix(~factuj-1)
  #Print(uj,factuj,Uj)
  #  attributes(Uj)$assign=NULL
  #  attributes(Uj)$contrasts=NULL
  #  colnames(Uj)=min:max
  U[locobs,fromP[j]:toP[j]]=as.integer(Uj)
 }
 if( is.null(iname) )
  iname=paste("Q",1:nitems,sep="")
 Pname=character(sum(ncat))
 catname=outer(paste(iname,"_",sep=""),0:max(ncat),paste,sep="")
 for( j in 1:nitems ){
  Pname[fromP[j]:toP[j]]=catname[j,1:ncat[j]]
 }
 colnames(U)=Pname
 rownames(U)=rownames(Uc)
 
 if( zero != 1 ){
  # remove the 0-th category
  U=U[,-fromP]
  fromP=fromP-0:(nitems-1)
  toP=toP-1:nitems
 }
 
 return( list(U=U, ncat=ncat, type=type, zero=zero, fromP=fromP, toP=toP) )
 
} # end of dummy_expand











checkparam <- function( param, type=NULL, modulename="checkparam" ){
 # check parameter data frame
 # Shin-ichi Mayekawa
 # 120224,27,29
 # name of giid: 120306
 # remove form/gfid: 120308
 # new types: 120911
 # type Bxx items: 120921
 #
 # Args:
 #
 #   param       parameter data frame name
 #   type        item type to be selected or NULL to select all
 #               B, G, P, PN or N  or ALL
 #   modulename  name of the module from which checkparam is called
 #
 #
 # Value:
 #
 #  parameter matrix of the form
 #     ncat p1 p2 ....
 #   or
 #  parameter data frame if type == "ALL"
 #
 #  If param is a matrix, param wll be returned as it is.
 #  If anything is wrong, NULL will be returned.
 #
 
 err=paste("\nerror(",modulename,")",sep="")
 
 if( is.null(param) ){
  cat(err, " Empty data frame or matrix:", substitute(param),"\n")
  return(NULL)
 }

 # param given as data frame or not
 if( is.data.frame(param) ){
  #
  #  Here, param is a data frame
  #
  # check if it is a genuine data frame
  if( is.null(param$name) & is.null(param$giid) ){
   cat(err, 
       " item param/weight data frame must have a variable named 'name'"
     , " or 'giid'.\n")
   return(NULL)
  }
  else if( is.null(param$name) ){
   param=data.frame(name=param$giid,param)
   param=param[,-which(colnames(param) == "giid")]
  }
  
  if( !is.null(param$gfid) ) param=param[,-which(colnames(param) == "gfid")]
  if( !is.null(param$form) ) param=param[,-which(colnames(param) == "form")]
  
  if( is.null(param$type) ){
   cat(err, 
       " item param/weight data frame must have a variable named 'type'.\n")
   return(NULL)
  }
  if( is.null(param$ncat) ){
   cat(err, 
       " item param/weight data frame must have a variable named 'ncat'.\n")
  return(NULL)
  }
  if( is.null(param$p1) ){
   cat(err, 
       " item param/weight data frame must have a variable named 'p1'.\n")
   return(NULL)
  }
  
  if( !is.null(type) ){
   if(  type == "ALL" ){
    # return the data frame as it is:   assuming checkparam is called from irf.
    return(param)
   }
   if( type == "B" ) loc=grep("^B[[:digit:]]*$", param$type)
   else              loc=which( param$type == type )
   if( length(loc) <= 0 ){
    cat(err, 
        " There are no type '",type,"' items in the data frame.\n")
    return(NULL)    
   }   
   param=param[loc,,drop=0]
  }
  
  # convert data frame to matrix
  locparamnum=grep("^p[[:digit:]]",colnames(param))
  lenp=ncol(param)
  parammat=as.matrix( cbind(param$ncat ,param[,locparamnum]) )
  colnames(parammat)=c("ncat",colnames(param)[locparamnum])

  if( !is.null(param$name) ) rownames(parammat)=param$name
  else rownames(parammat)=param$giid

 } # end of data frame
 
 else{
  # param is a matrix
  return(param)
 }
 
 return(parammat)
 
} # end of checkparam




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
 #     p2 is the intercept parameter c_{j1}
 #     p3 is the intercept parameter c_{j2}
 #
 #     p_ncat is the intercept parameter c_{j,[ncat-1]}
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
                   , stringsAsFactors=0
 )
 
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
 rsqrinfo=sqrt(1/info)
 if( debug > 0 ) Print(theta, stdx_t, info, rsqrinfo)
 
 # information function with the locally optimal category weights
 info_LO=rowSums( (dicrf^2)/icrf )
 rsqrinfo_LO=sqrt(1/info_LO)
 
 # information function with the locally optimal item weights
 info_LOW=rowSums( (dirf^2)/varu_t )
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
  matplot( theta, qt, type = "l", ylab="score"
           , xlim=c(min(theta),max(theta)), ylim=c(minscore_t,maxscore_t) 
           , main=title  )
  
  title=paste("Information Function defined as (slope of trf)^2 / var(x|theta)"
              ,", LO and LOW")
  minP=min(cbind(info,info_LO)); maxP=max(cbind(info,info_LO))
  matplot( theta, cbind(info,info_LO,info_LOW)
           , type = "l", ylab="information"
           , xlim=c(min(theta),max(theta)), ylim=c(minP,maxP) 
           , main=title  )
  title=paste("S.E. of Theta: inverse of info and info_LO")
  matplot( theta, cbind(rsqrinfo,rsqrinfo_LO), type = "l", ylab="score"
           , xlim=c(min(theta),max(theta)), ylim= 
            , main=title  )
 }
 
 # confidence interval
 # for each score, solve score = xU(theta) and score = xL(theta) for theta
 ci=matrix(NA,nscore_t,2)
 ci[,1]=approx(qt[,3],theta, xout=score, yleft=NA, yright=NA)$y
 ci[,2]=approx(qt[,1],theta, xout=score, yleft=NA, yright=NA)$y
 cid2=(ci[,2]-ci[,1])/2
 if( debug > 0 ) Print(score,ci,cid2)
 if( debug > 0 ) Print(qtd,slopetrf, digits=2)
 
 
 
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
           , type = "l", ylab="score"
           , ylim=c(thmin,thmax), xlim=c(minscore_t,maxscore_t) 
           , main=title  )
 }
 
 
 if( plot > 0 ){
  maint=paste("Posterior Mean and ", 100*(1-alpha)
              , "% CI of Theta given Observed Score")
  subt=paste("alpha =",alpha)
  matplot( score[locnmiss], cbind(meant_x,ci)[locnmiss,]
           , type = "l", ylab="score"
           , xlim=c(minscore_t,maxscore_t), ylim=c(thmin,thmax)
           , main=maint, sub=subt  )
  maint=paste(
   "Posterior Std and the Half Width of Post Quantile and CI")
  matplot( score[locnmiss], cbind(stdt_x, qthd2, cid2)[locnmiss,]
           , type = "l", ylab="score"
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
 theta_stat=as.data.frame( cbind(theta,Pt,trf,slopetrf,stdx_t,info
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
  Print(theta_stat, obs_stat, fmt=".2")
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



smn <- function( x, f, mu=NULL, sigma=NULL,    estmu=1, estsigma=1
               , maxiter=100, eps=1e-8, print=1, plot=0 ){
 # parametric scored multinomial: normal pdf model
 # Shin-ichi Mayekawa
 # 20121026
 # comments: 20121110(London)
 # plot: 121207
 #
 
 

 # const
 if( is.vector(x) ) m=length(x) else  m=nrow(x)
 n=sum(f)
 meanx=sum(x*f)/n
 stdx=sqrt( sum(x*x*f)/n-meanx^2 )
 P0=exp(-0.5*(x-meanx)^2/stdx^2)
 P0=P0/sum(P0)
 llh0=sum(f*log(P0))
 
 
 # initial values
 if( is.null(mu) )  mu=meanx
 if( is.null(sigma) ) sigma=stdx
 
 P=exp(-0.5*(x-mu)^2/sigma^2)
 P=P/sum(P)
 
 if( print > 0 ){
  cat("\n\nParametric Scored Multinomial Distribution: Normal PDF model\n\n")
  cat(" # of observations =", n,"\n")
  cat(" # of categories =", m,"\n")
  cat(" estimation of mu =", estmu, ",  estimation of sigma =", estsigma,"\n")
  cat(" max # of iterations =", maxiter," with eps =", eps,"\n")
  cat(" sample mean and std\n")
  cat("  mean =", meanx, ",  std =", stdx,"\n")  
  cat("  likelihood with the sample mean/std =", llh0,"\n")
 }
 
 
 llhp=llh0
 converged=0

 # main iteration
 for( llll in 1:maxiter ){
  
  # log-Jacobian matrix
  Jac=cbind( (x-mu)/sigma^2, (x-mu)^2/sigma^3 )
  
  # dispersion 
  Dr=n*( Diag(P)-P%*%t(P) )
  
  # gradient etc
  delta=matrix(t((f-n*P)),,1)
  g=t(Jac)%*%delta
  H=t(Jac)%*%Dr%*%Jac
  d=solve(H)%*%g
  
  # update parameters
  step=1; ok=0
  mu1=mu; sigma1=sigma;
  for( lllll in 1:20 ){
   if( estmu ) mu1=mu + step*d[1]
   if( estsigma ) sigma1=sigma + step*d[2]
   P1=exp(-0.5*(x-mu1)^2/sigma1^2)
   P1=P1/sum(P1)
   llh=sum(f*log(P1))
   if( llh >= llhp ){
    ok=1
    break
   }
   step=step/2
  }
  if( ok ){
   mu=mu1; sigma=sigma1
   P=P1   
  }
  
  # convergence
  llhimpr=(llh-llhp)/abs(llhp)
  maxag=max(abs(c(estmu,estsigma)*g))
  if( print >= 2 ) 
   Print(llll, llhp, llh, llhimpr,maxag)
  if( llhimpr <= eps  &&  maxag < eps ){
   converged=1
   break
  }
  
  # next iteration
  llhp=llh
  
  
 } # end of llll loop
 
 mean=sum(x*P); std=sqrt(sum(x*x*P)-mean^2)
 iter=llll
 
 # estimated density etc
 nP=n*P
 tab=cbind(x,f,nP,P)
 colnames(tab)=c("x","f","n*p","p")
 
 if( print >= 1 && converged )
  cat("\nIteration converged with ", iter, " iterations.  eps= ",eps,"\n")

 if( print > 0 ){
  cat("\n\nParametric Scored Multinomial Distribution: Normal PDF model\n\n")
  cat(" # of observations =", n,"\n")
  cat(" # of categories =", m,"\n")
  cat(" estimation of mu =", estmu, ",  estimation of sigma =", estsigma,"\n")
  cat("\n iterations required =", iter," with eps =", eps,"\n")
  cat(" log likelihood maximized =", llh, "\n")
  cat("\n estimated parameters of the Normal distribution\n")
  cat("  mu =", mu, ",  sigma =", sigma,"\n")
  cat("\n sample mean and std\n")
  cat("  mean =", meanx, ",  std =", stdx,"\n")  
  cat(" mean and std of the scored multinomial (preserved)\n")
  cat("  mean =", mean, ",  std =", std,"\n")
  cat("  likelihood with the sample mean/std =", llh0, "\n") 
  cat("\n final gradient\n")
  cat("  mu =", g[1], ",  sigma = ", g[2], "\n")
  if( print >= 3 ){
   cat("\ndata and the expected values (= n*prob)\n")
   print(tab)
  }
 }
 
 if( plot ){
  P=exp(-0.5*(x-mu1)^2/sigma1^2)
  sumP=sum(P)
  P=P/sum(P)
  x1=seq(min(x),max(x),length.out=51)
  P1=exp(-0.5*(x1-mu1)^2/sigma1^2)
  P1=P1/sumP
  nP1=n*P1
  maxy=max(nP,nP1,f)+5
  temp=barplot(height=f, names.arg=x, ylim=c(0,maxy)
      , main="Histogram of (x,f) with the Estimated Expected Values"
      , sub=paste("mu=",format(mu,digits=4)
                , ", sigma=", format(sigma,digits=4), "\n"
                , paste("  (mean=",format(meanx,digits=4)
                      , ", std=",format(stdx,digits=4), ")") )
                , offset=0  )
  par(new=1)
  plot( x, nP, ylim=c(0,maxy), xlim=c(min(x)-.5,max(x)+.5), type="p"
        , xlab="", ylab="", xaxt="n", yaxt="n" )
  par(new=1)
  plot( x1, nP1, ylim=c(0,maxy),  xlim=c(min(x)-.5,max(x)+.5), type="l"
        , xlab="", ylab="freq and n*prob", xaxt="n", yaxt="n" )
 }
 
 
 return( 
     list( P=P, x=x, mu=mu, sigma=sigma, estmu=estmu, estsigma=estsigma
         , llh=llh, mean=mean, std=std ) 
        ) 
  
 
} # end of smn




