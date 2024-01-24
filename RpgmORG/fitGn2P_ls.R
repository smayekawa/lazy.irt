#' Conversion of Partial Credit Items to Graded Response (2PNM) Items
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
#' paramG: GRM Item Parameter Data Frame (type="Gn")\cr
#' rmse: Vector or RMSEs \cr
#' grad: Gradient matrix \cr
#' wtype, wmean, wsd
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
#' paramGn1 <- fitGn2P_ls( paramP1, plot=1, print=1 )
#' paramG1 <- fitG2P_ls( paramP1, plot=1, print=1 )
#'
#'
#'
#' @export
#'

fitGn2P_ls <- function( paramP, wtype=0, wmean=0, wsd=1, DinP=1
                    , npoints=21, thmin=-3, thmax=3
                    , maxiter=100, eps=1e-6
                    , print=1, plot=0, debug=0 ){
 # fitting Graded Response Model to Generalized Partial Credit Model
 # Shin-ichi Mayekawa
 # fitG2L_ls: 20180113-20180114
 # fitG2L_ls modified: 20180211
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
  cat("\nFitting Graded Response Model (2PNM) to",
      " Generalized Partial Credit Model\n")
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
 paramG$type="Gn"

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
  cat("\n  Converted GRM \n")
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
               , " (wtype=", wtype, ")  [fitGn2P_ls]", sep="")
   sub=paste( "paramP ="
              , paste(format(paramP[j,3+(1:ncat[j])],digits=3)
                      , collapse=",  ")
              , "\nparamGn ="
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
                , " (wtype=", wtype, ")  [fitGn2P_ls]", sep="")
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

} # end of fitGn2P_ls


param=paramS2

res=fitG2P_ls( param, plot=1, print=1, maxiter=50, wtype=1 )

res=fitGn2P_ls( param, plot=1, print=1, maxiter=50, wtype=1 )





comments(
 '

 paramP=data.frame(name="P", type="P", ncat=3, p1=1, p2=-1, p3=1
, stringsAsFactors=0)



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

 paramG=paramG[c(7),]


 temp=fitG2P( paramG, npoints=21, plot=1, print=0, debug=0, wtype=1 )
 paramP1=temp$paramP; paramPN1=temp$paramPN; paramPN01=temp$paramPN0;
 paramG11=fitG2P( paramPN1, npoints=21, plot=1, print=0, debug=0, wtype=1 )
 paramG1=fitG2P( paramP1, npoints=21, plot=1, print=0, debug=0, wtype=1 )
 Print(paramG,paramG11,paramG1)







 res=fitG2P( paramG, npoints=21, plot=1, print=0, debug=1, wtype=1 )
 paramPN1=res$paramPN
 paramP1=res$paramP

 Print(paramP1,paramPN1,convP2PN(paramP1),convPN2P(paramPN1))

 res1=fitG2P( paramP1, npoints=21, plot=1, print=0 )
 Print(res1)
 res2=fitG2P( paramPN1, npoints=21, plot=1, print=0 )
 Print(res2)

 ')
