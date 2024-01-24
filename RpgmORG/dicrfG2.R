dicrfG <- function( param, theta, print=0, plot=0, debug=0 ){
 # derivative of GRM ICRF
 # Shin-ichi Mayekawa
 # 120214
 # checkparam: 120224,29
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
  dPPj=dicrfB( paramj, theta, print=print1 )
  
  # icbrf
  dPP[,fromPP[j]]=0;
  dPP[,(fromPP[j]+1):(toPP[j]-1)]=dPPj
  dPP[,toPP[j]]=0;
  
  # dicrf as difference of icbrf
  #    pj=PP[,fromPP[j]:toPP[j]-1]-PP[,fromPP[j]+1:toPP[j]];
  dP[,fromP[j]:toP[j]]=dPP[,fromPP[j]:(toPP[j]-1)]-dPP[,(fromPP[j]+1):toPP[j]];
  
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



comments('
paramG="
    ncat   p1   p2 p3 p3
QG1    3 1.05 -1.40 0.97 NA
QG2    3 1.59 -0.52 0.92 NA
QG3    4 0.93 -2.97 0.55  1.5
"; paramG=lazy.matrix(paramG,header=1)


npoints=31; minth=-4; maxth=4
theta=seq(minth,maxth,length.out=npoints)

res=dicrfG( paramG, theta, print=2, plot=0, debug=0 )
')

npoints=31; minth=-4; maxth=4
theta=seq(minth,maxth,length.out=npoints)

  dicrfG( param, theta, print=2, plot=2, debug=0 )

