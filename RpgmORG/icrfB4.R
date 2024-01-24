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




theta=seq(-4,4,length=21)
res2=icrfB( paramS1, theta, plot=1 )
