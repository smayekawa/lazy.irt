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
 # smallP: ???
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



