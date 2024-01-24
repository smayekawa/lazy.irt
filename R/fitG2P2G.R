#' Approximate Conversion of Partial Credit Items to Graded Response Items
#' Using Logit Transformation
#'
#' @param paramP Item Parameter Data Frame for Partial Credit Model
#' @param wtype = 1 to use normal weight, else uniform
#' @param DinP = 1 to include D=1.7 in logistic function
#' @param dataframe = 1 to create parameter data frame, not matrix
#' @param npoints # of discrete points for theta
#' @param thmin Minimum value of discrete thata value
#' @param thmax Maximum value of discrete thata value
#' @param maxabsv Maximum absolute value of the weight
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#' @param debug = 1 to print intemediate result
#'
#' @return
#' Graded Response Model Item Parameter Data Frame
#'
#' @examples
#' paramG1 <- fitG2P( paramS1[4,], npoints=21, plot=1, print=1 )
#'
#' @export
#'
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
 # bugfix: 201506011
 # legend, plot P not PN: 20150613
 # returns correct item name: 20180113
 # bugfix for wtype=0: 20180113
 # sort b parameters: 20180215
 # printerror: 20230329
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

 # store P w/o name and type
 paramPP=convPN2P( paramPN, DinP=1 )
 paramPP=paramPP[paramPP$type == "P",]
 paramPP=paramPP[,-(1:2)]

 # argument name
 paramname=as.character(substitute(paramPN))
 isdf=0
 if( is.data.frame(paramPN) ){
  # param and weight given as data frame
  isdf=1
  paramPN=checkparam( paramPN, "PN", "fitG2P", printerror=0 )
  if( is.null(paramPN) ){
   # cat("\n\nInput Item Parameter ERROR.\n\n")
   return()
  }
 }

 # const
 if( is.data.frame(paramPN) ) iname=paramPN$name
 else iname=rownames(paramPN)

 theta=seq(thmin,thmax,length.out=npoints)
 thetac=format(theta,digits=3)

 # weights for those elements with which abs(v) > maxabsv will be
 # set equal to min(W)

 nitems=nrow(paramPN)
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
 if( wtype == 0 ) W=matrix(1,npoints)
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
               , dimnames=list(iname,c("ncat",paste("p",1:max(ncat),sep=""))))
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
             log(1-PP[,(fromPP[j]+1):(toPP[j]-1)]),,1 )
  # adjust weight according to the magnitude of logit
  Wj=rep(W,ncat1[j])
  Wj[which(v < -maxabsv)]=Wmin; Wj[which(v > maxabsv)]=Wmin
  Wj=Diag(Wj)
  if( debug ) Print(j, U,v, vecdiag(Wj), digits=3)

  # solve
  beta=(1/1.7) * solve(t(U)%*%Wj%*%U)%*%t(U)%*%Wj%*%v
  if( debug ) Print(beta)
  paramG[j,2]=beta[1]
  bparam=sort( -beta[2:(2+ncat1[j]-1)]/beta[1] )
  paramG[j,3:(3+ncat1[j]-1)]=bparam

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
   sub=paste( "paramP ="
            , paste(format(paramPP[j,2:(ncat[j]+1)],digits=3)
                    , collapse=",  ")
            , "\nparamG ="
            , paste(format(paramG[j,2:(ncat[j]+1)],digits=3)
                    , collapse=",  ") )
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
   legend( range(theta)[1], 0.8, c("P","G"), cex=0.8, col=1:1
           ,  pch=c(NA,21), lty=c(1,2), title="model" )
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
  paramG=data.frame( name=iname
                     , type=rep("G",nitems), paramG, stringsAsFactors=0 )
 }


 return( paramG )

} # end of fitG2P





#' Approximate Conversion of Graded Response Items to Partial Credit Items
#' Using Logit Transformation
#'
#' @param paramG Item Parameter Data Frame for Graded Response Items
#' @param wtype = 1 to use normal weight, else uniform
#' @param DinP = 1 to include D=1.7 in logistic function
#' @param dataframe = 1 to create parameter data frame, not matrix
#' @param npoints # of discrete points for theta
#' @param thmin Minimum value of discrete thata value
#' @param thmax Maximum value of discrete thata value
#' @param maxabsv Maximum absolute value of the weight
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#' @param debug = 1 to print intemediate result
#'
#' @return A list of:
#' Partial Credit Model Item Parameter Data Frame in Nominal Format \cr
#' Partial Credit Model Item Parameter Data Frame in Nominal Format0 \cr
#' Partial Credit Model Item Parameter Data Frame in Standard Format
#'
#'
#' @examples
#' paramP1 <- fitP2G( paramS1[3,], npoints=21, plot=1, print=1 )$paramP
#' paramG1 <- fitG2P( paramP1, npoints=21, plot=1, print=1 )
#'
#'
#'
#' @export
#'

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
 # bugfix: 201506011,12
 # legend, plot P not PN: 20150613
 # returns correct item name: 20180113
 # printerror: 20230329
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
  paramG=checkparam( paramG, "G", "fitP2G", printerror=0 )
  if( is.null(paramG) ){
   # cat("\n\nInput Item Parameter ERROR.\n\n")
   return()
  }
 }

 # const
 if( is.data.frame(paramG) ) iname=paramG$name
 else iname=rownames(paramG)

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
 paramP[,2]=paramPN[,2]
 if( DinP == 0 ) paramP[,2]=paramP[,2]/1.7
 b=paramPN[,3:ncol(paramPN), drop=0]
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
              , paste(format(paramG[j,2:(ncat[j]+1)],digits=3)
                      , collapse=",  ")
              , "\nparamP ="
              , paste(format(paramP[j,2:(ncat[j]+1)],digits=3)
                      , collapse=",  ")
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
   legend( range(theta)[1], 0.8, c("G","P"), cex=0.8, col=1:1
           ,  pch=c(NA,21), lty=c(1,2), title="model" )
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
  paramPN=data.frame( name=iname
                      , type=rep("PN",nitems), paramPN, stringsAsFactors=0 )
  paramPN0=data.frame( name=iname
                     , type=rep("PN0",nitems), paramPN0, stringsAsFactors=0 )
  paramP=data.frame( name=iname
                   , type=rep("P",nitems), paramP, stringsAsFactors=0 )
 }

 return( list(paramPN=paramPN, paramPN0=paramPN0, paramP=paramP) )

} # end of fitP2G


