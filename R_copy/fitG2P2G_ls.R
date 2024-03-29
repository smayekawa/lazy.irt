#' Conversion of 3PLM Items to 2PLM Items
#'
#' @param param3 Item Parameter Data Frame for 3PLM Items
#' @param wtype = 1 to use dnorm(theta) as the weight
#' @param wmean The mean of normal distribution to be used as the weight
#' @param wsd The sd of normal distribution to be used as the weight
#' @param DinP = 1 to include D=1.7 in logistic function
#' @param npoints # of discrete points for theta
#' @param thmin Minimum value of discrete thata value
#' @param thmax Maximum value of discrete thata value
#' @param maxiter Maximum # of GN iterations
#' @param eps  Convergence criterion for the relative improvement of rmse
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#' @param debug = 1 to print intemediate result
#'
#' @return A list of: \cr
#' param3: Input 3PLM item parameter data frame (subset, type="P")\cr
#' param2: 2PLM Item Parameter Data Frame (type="G")\cr
#' rmse: Vector of sqrt(rss/length(theta)) for each item.
#'
#' @details
#' This function minimizes\cr
#' \code{ rss=sum( w*( vec(icrf_3(theta)) - vec(icrf_2(theta)) )^2 ) } \cr
#' with respec to the 2PLM item parameters, \cr
#' where \code{icrf_3(theta)} is the icrf of input 3PLM items and
#'  \code{icrf_2(theta)} is the icrf of fitted 2PLM items, \cr
#' and \code{w} is the weight vector normal or 1. \cr\cr
#' Unlike fitP2G_ls and fitG2P_ls, the icrf of the 0-th category is not used.
#' \cr
#' Weighted Gauss-Newton method is used for the minimization.
#' \cr
#'
#'
#' @examples
#' param2 <- fit223_ls( paramS2, plot=1, print=1 )
#' param21 <- fit223_ls( paramS2, plot=1, print=1, wtype=1 )
#'
#'
#' @export
#'

fit223_ls <- function( param3, wtype=0, wmean=0, wsd=1, DinP=1
                       , npoints=21, thmin=-3, thmax=3
                       , maxiter=100, eps=1e-6
                       , print=1, plot=0, debug=0 ){
 # fitting 2PLM to 3PLM
 # Shin-ichi Mayekawa
 # original version is fitG2P_ls
 # 223: 20180114
 # ncatj removed from the den of rmse: 20180116
 #
 # Args:
 #
 #     param2   parameter data frame for type = "G" items
 #     DinP = 1 to use 1.7 in Exp of partial credit model
 #
 #     npoints  # of theta points in [thmin, thmax]
 #     wtype = 1 to use normal weight, else uniform
 #
 #     dataframe = 1 to create parameter data frame, not matrix
 #
 # Value
 #   list of two types of item paramter matrices
 #     param3   =  priginal parameter set with DinP
 #
 #
 # Needs:
 #   irf, fitG2P
 #
 #
 #


 # argument name
 paramname=as.character(substitute(param3))

 # pick up 3PLM items
 isdf=0
 if( is.data.frame(param3) ){
  # param
  isdf=1
  param3=param3[param3[,"type"]=="B3",]
  if( is.null(param3) ){
   cat("\n\nInput Item Parameter ERROR.\n\n")
   return()
  }
 }
 param3=param3[,1:6]

 # const
 if( is.data.frame(param3) ) iname=param3$name
 else iname=rownames(param3)

 # const
 theta=seq(thmin,thmax,length.out=npoints)
 thetac=format(theta,digits=3)

 # weights for those elements with which abs(v) > maxabsv will be
 # set equal to min(W)

 nitems=nrow(param3)
 ncat=param3[,"ncat"]
 ncat1=ncat-1

 if( print >= 1 ){
  cat("\nFitting 2PLM to 3PLM\n")
  cat("  # of items =", nitems, ", weight type =", wtype, sep="" )
  if( wtype == 1 ) cat(" with N(",wmean,",",wsd,")", sep="")
  cat("\n")
  cat("  # of theta points =", npoints, " in ["
      , thmin, " , ", thmax, "]\n")
  cat("  D in Partial Credit Model =", DinP,"\n")
 }
 if( debug ) Print(param3, ncat)

 # weight
 if( wtype == 0 ) w=rep(1,npoints)
 else if( wtype == 1 ){
  w=dnorm((theta-wmean)/wsd)
 }
 w=w/sum(w)
 w2=sqrt(w)
 wmin=min(w); wmax=max(w);

 # GPCM icrf matrix:  This is Y
 res=irf( param3, theta, print=0, zero=0 )
 P=res$ICRF
 toP=res$toP; fromP=res$fromP
 rm(res)
 Pname=colnames(P)

 # initial
 param2=param3
 param2$p3=0
 param2$type="B"

 rmse=numeric(nitems)
 names(rmse)=iname
 grad=matrix(NA,max(ncat),nitems)
 rownames(grad)=paste("p",1:nrow(grad),sep="")
 colnames(grad)=iname
 grad=t(grad)

 # for each GPCM item
 for( j in 1: nitems ){

  if( print >= 1 ) cat("\nItem # ", j, "\n", sep="")

  ncatj=ncat[j]
  Pj=P[,fromP[j]:toP[j]]
  vpj=vec(Pj)
  paramjdf=param2[j,]
  paramj=paramjdf[1,3+(1:ncatj)]
  paramjp=paramj
  vpjhat=c( irf( paramjdf, theta, print=0, zero=0 )$ICRF )
  rmsej=sqrt( sum( w*(vpj-vpjhat)^2 )/npoints )
  rmsejp=rmsej


  # GN iteration
  for( llll in  1:maxiter ){

   # Jacobian:  d vec(Pj) / d paramj    ncatj*npoints x ncatj
   # Jac=JacobianMat( paramj, icrf, ..eps.. = 1e-06 )
   paramjdf[1,3+(1:ncatj)]=paramj
   Jac=w2*dirf_p( paramjdf, theta, print=0, zero=0 )$Jack
   # Print(Jac-Jac0)
   # gradient and info mat and GN direction
   g=-t( t(w*(vpj-vpjhat))%*%Jac )
   maxag=max(abs(g))
   H=t(Jac)%*%Jac
   d=Ginv(H, MASS=1)%*%g

   # Print(Jac,g,d)

   ok=0
   step=1
   for( lllll in 1:20 ){
    pj=paramj - step*d
    paramjdf[1,3+(1:ncatj)]=pj
    vpjhat=c( irf( paramjdf, theta, print=0, zero=0 )$ICRF )
    rmse1=sqrt( sum( w*(vpj-vpjhat)^2 )/npoints )
    # Print(llll, lllll, rmsej, rmse1)
    if( rmse1 <= rmsej ){
     ok=1
     paramj=pj
     rmsej=rmse1
     break
    }
    step=step/2
   }
   if( ok == 0 && print >= 1 ){
     Print("Halving failed! ",  rmsej, rmse1, max(abs(d)), fmt="12.7", "\n")
   }


   # check convergence
   rmseimpr=(rmsejp-rmsej)/rmsejp
   maxadp=max(abs(paramjp-paramj))
   if( print >= 2 )
    Print( llll, lllll, rmsej, rmsejp, rmseimpr, maxadp, maxag
           , fmt="i3 i2 .6")

   if( rmseimpr <= eps  &&  maxadp <= eps ) break

   # next iteration
   rmsejp=rmsej
   paramjp=paramj


  } # end of llll loop


  if( print >= 1 ){
   cat("\n Iteration terminated:\n")
   Print( llll, lllll, rmsej, rmsejp, rmseimpr, maxadp, maxag
          , fmt="i3 i2 .6")
  }

  # truncate
  paramj=trunc(paramj*10000)/10000

  param2[j,3+(1:ncatj)]=paramj
  rmse[j]=rmsej
  grad[j,1:ncat[j]]=g


 } # end of j loop for item




 rmsec=format(rmse,digits=5)
 P2=irf( param2, theta, print=0, zero=0 )$ICRF
 diffP=matrix(colMeans(abs(P-P2)),,1,dimnames=list(colnames(P2),""))
 if( print ){
  cat("\nInput and  Converted Item Parameters\n")
  cat("\nInput and  3PLM Item Parameters\n")
  Print(param3)
  cat("\n  Converted 2PLM D=1.7\n")
  Print(param2)
  cat("Gradients and RMSE minimized\n")
  Print(grad, rmse, fmt="8.5")
  if( print >= 2 ){
   cat("\n ICRFs of the original 3PLM\n")
   Print(P, fmt="5,3")
   cat(" ICRFs of the converted 2PLM\n")
   Print(P2, fmt="5,3")
  }
 }

 if( plot ){
  for( j in 1:(nitems) ){
   title=paste("ICRF of ",iname[j]," (ncat=",ncat[j]
               , ")  rmse=", rmsec[j]
               , " (wtype=", wtype, ")  [fit223_ls]", sep="")
   sub=paste( "param3 ="
              , paste(format(param3[j,3+(1:3)],digits=3)
                      , collapse=",  ")
              , "\nparam2 ="
              , paste(format(param2[j,3+(1:ncat[j])],digits=3)
                      , collapse=",  ")
   )
   plot(theta,P[,fromP[j]], xlab=""
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="")
   par(new=1)
   plot(theta,P2[,toP[j]], main=title, sub=sub, xlab=""
        , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="b", ylab="ICRF")
   legend( range(theta)[1]+0.5, 0.9, c("3","2"), cex=0.8, col=1:1
           ,  pch=c(NA,21), lty=c(1,2), title="model" )
   par(new=0)
  }

  if( plot >= 2 ){
   for( j in 1:(nitems) ){
    difP=P-P2
    maxd=max(difP); mind=min(difP)
    title=paste("ICRF Differences of ",iname[j]," (ncat=",ncat[j]
                , ")  rmse=", rmsec[j]
                , " (wtype=", wtype, ")  [fit223_ls]", sep="")
    plot(theta,difP[,toP[j]], main=title, sub="3PLM - 2PLM(fitted)"
         , xlim=c(min(theta),max(theta)), ylim=c(mind,maxd), type="l"
         , ylab="difference")
    par(new=0)
   }
  }
 }

 # convert the result to data frame

 return( named_list(param2, param3, rmse, grad, wtype, wmean, wsd) )

} # end of fit223_ls










#' Conversion of Partial Credit Items to Graded Response Items
#'
#' @param paramP Item Parameter Data Frame for GPCM Items
#' @param wtype = 1 to use dnorm(theta) as the weight
#' @param wmean The mean of normal distribution to be used as the weight
#' @param wsd The sd of normal distribution to be used as the weight
#' @param DinP = 1 to include D=1.7 in logistic function
#' @param npoints # of discrete points for theta
#' @param thmin Minimum value of discrete thata value
#' @param thmax Maximum value of discrete thata value
#' @param maxiter Maximum # of GN iterations
#' @param eps  Convergence criterion for the relative improvement of rmse
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#' @param debug = 1 to print intemediate result
#'
#' @return A list of: \cr
#' paramP: Input GPCM item parameter data frame (subset, type="P")\cr
#' paramG: GRM Item Parameter Data Frame (type="G")\cr
#' rmse: Vector or RMSEs
#'
#' @details
#' This function minimizes\cr
#' \code{ sum( w*( vec(icrf_P(theta)) - vec(icrf_GRM(theta)) )^2 ) } \cr
#' with respec to the GRM item parameters, \cr
#' where \code{icrf_P(theta)} is the icrf of input GPCM items and
#'  \code{icrf_GRM(theta)} is the icrf of fitted GRM items, \cr
#' and \code{w} is the weight vector normal or 1.
#' \cr
#' Weighted Gauss-Newton method is used for the minimization.
#' \cr
#'
#'
#' @examples
#' paramP1 <- fitP2G_ls( paramS2, plot=1, print=1 )$paramP
#' paramG1 <- fitG2P_ls( paramP1, plot=1, print=1 )
#'
#'
#'
#' @export
#'

fitG2P_ls <- function( paramP, wtype=0, wmean=0, wsd=1, DinP=1
                    , npoints=21, thmin=-3, thmax=3
                    , maxiter=100, eps=1e-6
                    , print=1, plot=0, debug=0 ){
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
 # LS version: 20180113
 # JacobianMat -> dirf_p: 2018114
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
 #
 #
 # Needs:
 #   irf, fitG2P
 #
 #
 #
 # Note that in this function, paramG and paramP are data frames,
 # not numeric matrices as used in fitG2P and fitG2P functions.
 #
 #


 icrf <- function( paramj ){
  # to be used in JacobianMat
  # everything else is from outside
  paramjdf[1,3+(1:ncatj)]=paramj
  vpjhat=w2*c( irf( paramjdf, theta, print=0 )$ICRF )
  return(vpjhat)
 } # end of icrf

 get_rmse <- function( paramj ){
  # everything else is from outside
  paramjdf[1,3+(1:ncatj)]=paramj
  vpjhat=c( irf( paramjdf, theta, print=0 )$ICRF )
  rmse=sqrt( sum( w*(vpj-vpjhat)^2 /npoints/ncatj ) )
  return( rmse )
 } # end of get_rmse



 # argument name
 paramname=as.character(substitute(paramP))

 # pick up GRM items
 isdf=0
 if( is.data.frame(paramP) ){
  # param
  isdf=1
  paramP=paramP[paramP[,"type"]=="P",]
  if( is.null(paramP) ){
   cat("\n\nInput Item Parameter ERROR.\n\n")
   return()
  }
 }

 # const
 if( is.data.frame(paramP) ) iname=paramP$name
 else iname=rownames(paramP)

 # const
 theta=seq(thmin,thmax,length.out=npoints)
 thetac=format(theta,digits=3)

 # weights for those elements with which abs(v) > maxabsv will be
 # set equal to min(W)

 nitems=nrow(paramP)
 ncat=paramP[,"ncat"]
 ncat1=ncat-1

 if( print >= 1 ){
  cat("\nFitting Graded Response Model to Generalized Partial Credit Model\n")
  cat("  # of items =", nitems, ", weight type =", wtype, sep="" )
  if( wtype == 1 ) cat(" with N(",wmean,",",wsd,")", sep="")
  cat("\n")
  cat("  # of theta points =", npoints, " in ["
      , thmin, " , ", thmax, "]\n")
  cat("  D in Partial Credit Model =", DinP,"\n")
 }
 if( debug ) Print(paramP, ncat)

 # weight
 if( wtype == 0 ) w=rep(1,npoints)
 else if( wtype == 1 ){
  w=dnorm((theta-wmean)/wsd)
  w=w/sum(w)
 }
 w2=sqrt(w)
 wmin=min(w); wmax=max(w);

 # GPCM icrf matrix:  This is Y
 res=icrfP( paramP, theta )
 P=res$P
 toP=res$toP; fromP=res$fromP
 rm(res)
 Pname=colnames(P)

 # initial
 paramG=fitG2P( paramP, npoints=21, plot=0, print=0, debug=0, wtype=wtype )

 rmse=numeric(nitems)
 names(rmse)=iname
 grad=matrix(NA,max(ncat),nitems)
 rownames(grad)=paste("p",1:nrow(grad),sep="")
 colnames(grad)=iname
 grad=t(grad)

 # for each GPCM item
 for( j in 1: nitems ){

  if( print >= 1 ) cat("\nItem # ", j, "\n", sep="")

  ncatj=ncat[j]
  Pj=P[,fromP[j]:toP[j]]
  vpj=vec(Pj)
  paramjdf=paramG[j,]
  paramj=paramjdf[1,3+(1:ncatj)]
  paramjp=paramj
  vpjhat=c( irf( paramjdf, theta, print=0 )$ICRF )
  rmsej=get_rmse( paramj )
  rmsejp=rmsej


  # GN iteration
  for( llll in  1:maxiter ){

   # Jacobian:  d vec(Pj) / d paramj    ncatj*npoints x ncatj
   # Jac=JacobianMat( paramj, icrf, ..eps.. = 1e-06 )
   paramjdf[1,3+(1:ncatj)]=paramj
   Jac=w2*dirf_p( paramjdf, theta, print=0, zero=1 )$Jack
   # Print(Jac-Jac0)
   # gradient and info mat and GN direction
   g=-t( t(w*(vpj-vpjhat))%*%Jac )
   maxag=max(abs(g))
   H=t(Jac)%*%Jac
   d=Ginv(H, MASS=1)%*%g

   # Print(Jac,g,d)

   ok=0
   step=1
   for( lllll in 1:20 ){
    pj=paramj - step*d
    paramjdf[1,3+(1:ncatj)]=pj
    vpjhat=c( irf( paramjdf, theta, print=0 )$ICRF )
    rmse1=sqrt( sum( w*(vpj-vpjhat)^2 )/npoints/ncatj )
    # Print(llll, lllll, rmsej, rmse1)
    if( rmse1 <= rmsej ){
     ok=1
     paramj=pj
     rmsej=rmse1
     break
    }
    step=step/2
   }
   if( ok == 0 && print >= 1 ){
     Print("Halving failed! ",  rmsej, rmse1, max(abs(d)), fmt="12.7", "\n")
   }


   # check convergence
   rmseimpr=(rmsejp-rmsej)/rmsejp
   maxadp=max(abs(paramjp-paramj))
   if( print >= 2 )
    Print( llll, lllll, rmsej, rmsejp, rmseimpr, maxadp, maxag
           , fmt="i3 i2 .6")

   if( rmseimpr <= eps  &&  maxadp <= eps ) break

   # next iteration
   rmsejp=rmsej
   paramjp=paramj


  } # end of llll loop


  if( print >= 1 ){
   cat("\n Iteration terminated:\n")
   Print( llll, lllll, rmsej, rmsejp, rmseimpr, maxadp, maxag
          , fmt="i3 i2 .6")
  }

  # truncate
  paramj=trunc(paramj*10000)/10000

  paramG[j,3+(1:ncatj)]=paramj
  rmse[j]=rmsej
  grad[j,1:ncat[j]]=g


 } # end of j loop for item




 rmsec=format(rmse,digits=5)
 PP=irf( paramP, theta, print=0 )$ICRF
 diffP=matrix(colMeans(abs(P-PP)),,1,dimnames=list(colnames(PP),""))
 if( print ){
  cat("\nInput and  Converted Item Parameters\n")
  cat("\nInput and  GPCM Item Parameters\n")
  Print(paramP)
  cat("\n  Converted GRM D=1.7\n")
  Print(paramG)
  cat("Gradients and RMSE minimized\n")
  Print(grad, rmse, fmt="8.5")
  if( print >= 2 ){
   cat("\n ICRFs of the original GPCM\n")
   Print(P, fmt="5,3")
   cat(" ICRFs of the converted GRM\n")
   Print(PP, fmt="5,3")
  }
 }

 if( plot ){
  for( j in 1:(nitems) ){
   title=paste("ICRF of ",iname[j]," (ncat=",ncat[j]
               , ")  rmse=", rmsec[j]
               , " (wtype=", wtype, ")  [fitG2P_ls]", sep="")
   sub=paste( "paramP ="
              , paste(format(paramP[j,3+(1:ncat[j])],digits=3)
                      , collapse=",  ")
              , "\nparamG ="
              , paste(format(paramG[j,3+(1:ncat[j])],digits=3)
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
   legend( range(theta)[1]+0.5, 0.9, c("P","G"), cex=0.8, col=1:1
           ,  pch=c(NA,21), lty=c(1,2), title="model" )
   par(new=0)
  }

  if( plot >= 2 ){
   for( j in 1:(nitems) ){
    difP=P-PP
    maxd=max(difP); mind=min(difP)
    title=paste("ICRF Differences of ",iname[j]," (ncat=",ncat[j]
                , ")  rmse=", rmsec[j]
                , " (wtype=", wtype, ")  [fitG2P_ls]", sep="")
    for( k in 1:(ncat1[j]-1) ){
     plot(theta,difP[,fromP[j]+k-1]
          , xlim=c(min(theta),max(theta)), ylim=c(mind,maxd)
          , type="l", ylab="")
     par(new=1)
    }
    plot(theta,difP[,toP[j]], main=title, sub="GPCM - GRM(fitted)"
         , xlim=c(min(theta),max(theta)), ylim=c(mind,maxd), type="l"
         , ylab="difference")
    par(new=0)
   }
  }
 }

 # convert the result to data frame

 return( named_list(paramG, paramP, rmse, grad, wtype, wmean, wsd) )

} # end of fitG2P_ls





#' Conversion of Graded Response Items to Partial Credit Items
#'
#' @param paramG Item Parameter Data Frame for Graded Response Items
#' @param wtype = 1 to use dnorm(theta) as the weight
#' @param wmean The mean of normal distribution to be used as the weight
#' @param wsd The sd of normal distribution to be used as the weight
#' @param DinP = 1 to include D=1.7 in logistic function
#' @param npoints # of discrete points for theta
#' @param thmin Minimum value of discrete thata value
#' @param thmax Maximum value of discrete thata value
#' @param maxiter Maximum # of GN iterations
#' @param eps  Convergence criterion for the relative improvement of rmse
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#' @param debug = 1 to print intemediate result
#'
#' @return A list of: \cr
#' paramG: Input GRM item parameter data frame (subset)\cr
#' paramP: GPCM Item Parameter Data Frame (type="P")\cr
#' rmse: Vector or RMSEs
#'
#' @details
#' This function minimizes\cr
#' \code{ sum( w*( vec(icrf_GRM(theta)) - vec(icrf_P(theta)) )^2 ) } \cr
#' with respec to the GPCM item parameters, \cr
#' where \code{icrf_GRM(theta)} is the icrf of input GRM items and
#'  \code{icrf_P(theta)} is the icrf of fitted GPCM items, \cr
#' and \code{w} is the weight vector normal or 1.
#' \cr
#' Weighted Gauss-Newton method is used for the minimization.
#' \cr
#'
#'
#' @examples
#' paramP1 <- fitP2G_ls( paramS2, plot=1, print=1 )$paramP
#' paramG1 <- fitG2P_ls( paramP1, plot=1, print=1 )
#'
#'
#'
#' @export
#'

fitP2G_ls <- function( paramG, wtype=0, wmean=0, wsd=1, DinP=1
                    , npoints=21, thmin=-3, thmax=3
                    , maxiter=100, eps=1e-6
                    , print=1, plot=0, debug=0 ){
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
 # LS version: 20180112,13
 # JacobianMat -> dirf_p: 2018114
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
 #
 #
 # Needs:
 #   irf, fitP2G
 #
 #
 #
 # Note that in this function, paramG and paramP are data frames,
 # not numeric matrices as used in fitP2G and fitG2P functions.
 #
 #


 icrf <- function( paramj ){
  # to be used in JacobianMat
  # everything else is from outside
  paramjdf[1,3+(1:ncatj)]=paramj
  vpjhat=w2*c( irf( paramjdf, theta, print=0 )$ICRF )
  return(vpjhat)
 } # end of icrf

 get_rmse <- function( paramj ){
  # everything else is from outside
  paramjdf[1,3+(1:ncatj)]=paramj
  vpjhat=c( irf( paramjdf, theta, print=0 )$ICRF )
  rmse=sqrt( sum( w*(vpj-vpjhat)^2 /npoints/ncatj ) )
  return( rmse )
 } # end of get_rmse



 # argument name
 paramname=as.character(substitute(paramG))

 # pick up GRM items
 isdf=0
 if( is.data.frame(paramG) ){
  # param
  isdf=1
  paramG=paramG[paramG[,"type"]=="G",]
  if( is.null(paramG) ){
   cat("\n\nInput Item Parameter ERROR.\n\n")
   return()
  }
 }

 # const
 if( is.data.frame(paramG) ) iname=paramG$name
 else iname=rownames(paramG)

 # const
 theta=seq(thmin,thmax,length.out=npoints)
 thetac=format(theta,digits=3)


 nitems=nrow(paramG)
 ncat=paramG[,"ncat"]
 ncat1=ncat-1

 if( print >= 1 ){
  cat("\nFitting Generalized Partial Credit Model to Graded Response Model\n")
  cat("  # of items =", nitems, ", weight type =", wtype, sep="" )
  if( wtype == 1 ) cat(" with N(",wmean,",",wsd,")", sep="")
  cat("\n")
  cat("  # of theta points =", npoints, " in ["
      , thmin, " , ", thmax, "]\n", sep="")
  cat("  D in Partial Credit Model =", DinP,"\n", sep="")
 }
 if( debug ) Print(paramG, ncat)

 # weight
 if( wtype == 0 ) w=rep(1,npoints)
 else if( wtype == 1 ){
  w=dnorm((theta-wmean)/wsd)
  w=w/sum(w)
 }
 w2=sqrt(w)
 wmin=min(w); wmax=max(w)

 # Graded icrf matrix:  This is Y
 res=icrfG( paramG, theta )
 P=res$P
 toP=res$toP; fromP=res$fromP
 rm(res)
 Pname=colnames(P)

 # initial
 temp=fitP2G( paramG, npoints=21, plot=0, print=0, debug=0, wtype=wtype )
 paramP=temp$paramP

 rmse=numeric(nitems)
 names(rmse)=iname
 grad=matrix(NA,max(ncat),nitems)
 rownames(grad)=paste("p",1:nrow(grad),sep="")
 colnames(grad)=iname
 grad=t(grad)


 # for each GRM item
 for( j in 1: nitems ){

  if( print >= 1 ) cat("\nItem # ", j, "\n", sep="")

  ncatj=ncat[j]
  Pj=P[,fromP[j]:toP[j]]
  vpj=vec(Pj)
  paramjdf=paramP[j,]
  paramj=paramjdf[1,3+(1:ncatj)]
  paramjp=paramj
  vpjhat=c( irf( paramjdf, theta, print=0 )$ICRF )
  rmsej=get_rmse( paramj )
  rmsejp=rmsej


  # GN iteration
  for( llll in  1:maxiter ){

   # Jacobian:  d vec(Pj) / d paramj    ncatj*npoints x ncatj
   # Jac=JacobianMat( paramj, icrf, ..eps.. = 1e-06 )
   paramjdf[1,3+(1:ncatj)]=paramj
   Jac=w2*dirf_p( paramjdf, theta, print=0, zero=1 )$Jack
   # gradient and info mat and GN direction
   g=-t( t(w*(vpj-vpjhat))%*%Jac )
   maxag=max(abs(g))
   H=t(Jac)%*%Jac
   d=Ginv(H, MASS=1)%*%g

   # Print(Jac,g,d)

   ok=0
   step=1
   for( lllll in 1:10 ){
    pj=paramj - step*d
    paramjdf[1,3+(1:ncatj)]=pj
    vpjhat=c( irf( paramjdf, theta, print=0 )$ICRF )
    rmse1=sqrt( sum( w*(vpj-vpjhat)^2 )/npoints/ncatj )
    # Print(llll, lllll, rmsej, rmse1)
    if( rmse1 <= rmsej ){
     ok=1
     paramj=pj
     rmsej=rmse1
     break
    }
    step=step/2
   }
   if( ok == 0 && print >= 1 ){
    Print("Halving failed! \n", rmsej, rmse1, fmt="12.7", "\n")
    # Print(X,W,rmse1)
   }


   # check convergence
   rmseimpr=(rmsejp-rmsej)/rmsejp
   maxadp=max(abs(paramjp-paramj))
   if( print >= 2 )
    Print( llll, lllll, rmsej, rmsejp, rmseimpr, maxadp, maxag
           , fmt="i3 i2 .6")

   if( rmseimpr <= eps  &&  maxadp <= eps ) break

   # next iteration
   rmsejp=rmsej
   paramjp=paramj


  } # end of llll loop


  if( print >= 1 ){
   cat("\n Iteration terminated:\n")
   Print( llll, lllll, rmsej, rmsejp, rmseimpr, maxadp, maxag
          , fmt="i3 i2 .6")
  }

  # truncate
  paramj=trunc(paramj*10000)/10000

  paramP[j,3+(1:ncatj)]=paramj
  rmse[j]=rmsej
  grad[j,1:ncat[j]]=g


 } # end of j loop for item




 rmsec=format(rmse,digits=5)
 PP=irf( paramP, theta, print=0 )$ICRF
 diffP=matrix(colMeans(abs(P-PP)),,1,dimnames=list(colnames(PP),""))
 if( print ){
  cat("\nInput and  Converted Item Parameters\n")
  cat("\nInput and  Graded Response Item Parameters\n")
  Print(paramG)
  cat("\n  Converted parameters in original format with D=1.7\n")
  Print(paramP)
  cat("Gradients and RMSE minimized\n")
  Print(grad, rmse, fmt="8.5")
  if( print >= 2 ){
   cat("\n ICRFs of the original Graded response Model\n")
   Print(P, fmt="5,3")
   cat(" ICRFs of the converted Partial Credit Model\n")
   Print(PP, fmt="5,3")
  }
 }

 if( plot ){
  for( j in 1:(nitems) ){
   title=paste("ICRF of ",iname[j]," (ncat=",ncat[j]
               , ")  rmse=", rmsec[j]
               , " (wtype=", wtype, ")  [fitP2G_ls]", sep="")
   sub=paste( "paramG ="
              , paste(format(paramG[j,3+(1:ncat[j])],digits=3)
                      , collapse=",  ")
              , "\nparamP ="
              , paste(format(paramP[j,3+(1:ncat[j])],digits=3)
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
   legend( range(theta)[1]+0.5, 0.9, c("G","P"), cex=0.8, col=1:1
           ,  pch=c(NA,21), lty=c(1,2), title="model" )
   par(new=0)
  }

  if( plot >= 2 ){
   for( j in 1:(nitems) ){
    difP=P-PP
    maxd=max(difP); mind=min(difP)
    title=paste("ICRF Differences of ",iname[j]," (ncat=",ncat[j]
                , ")  rmse=", rmsec[j]
                , " (wtype=", wtype, ")  [fitP2G_ls]", sep="")
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

 return( named_list(paramP, paramG, rmse, grad, wtype, wmean, wsd) )

} # end of fitP2G_ls


