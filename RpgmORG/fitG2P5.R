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








paramG1=fitG2P( paramS2, npoints=21, plot=1, print=1, wtype=0 )




comments('
temp="
name type ncat p1 p2 p3 p4
Q1  G 3  1.0 -2  0 NA
Q2  G 3  1.0 -1  1 NA
Q3  G 3  1.0  0  2 NA
Q4  G 3  0.7 -3  1 NA
Q5  G 3  0.7 -2  2 NA
Q6  G 3  0.7 -1  3 NA
Q7  G 4  1.0 -2  0  2
Q8  G 4  0.7 -2  0  2
"; paramG=cards(temp,header=1)

#paramG=paramG[1,]


temp=fitP2G( paramG, npoints=21, plot=0, print=0 )
paramP1=temp$paramP

paramG1=fitG2P( paramP1, npoints=21, plot=1, print=1, debug=1 )


paramG11=fitG2P( paramPN1, npoints=21, plot=1, print=0, debug=0, wtype=1 )
paramG1=fitG2P( paramP1, npoints=21, plot=1, print=0, debug=0, wtype=1 )



paramP=data.frame(name="Q5", type="P", ncat=4, p1=0.9, p2=-2, p3=0, p4=2
                  , stringsAsFactors=0)
')
