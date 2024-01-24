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
#' grad: Gradient matrix \cr
#' wtype, wmean, wsd
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
 # roxgen2: 20180211
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












res=fit223_ls( paramS2, plot=2, print=2, maxiter=50, wtype=0 )



comments(
 '

set.seed(1701)

n=100
p1=0.6+runif(n)
p2=rnorm(n)
p3=0.5*runif(n)

paramB3=data.frame( name=paste("Q",1:n,sep=""), type="B3", ncat=2,
                    p1=p1, p2=p2, p3=p3, stringsAsFactors=0 )

res2 <- fit223_ls( paramB3, plot=1, print=0, wtype=1 )

mean(res2$rmse)


 ')
