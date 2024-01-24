#' Conversion to Normal Graded Response Model
#' 
#' Conversion to Normal Graded Response Model using Weighted Least Squares \cr
#' Japanese help file: (\link[lazy.irt]{conv2Gn_JPH})
#'
#' @param param Item Parameter Data Frame with item types \cr
#'   "B3","Bn","Bn3","P", "G"
#' @param theta Vector of theta points
#' @param init  = 1 to use fitG2P, else use equally spaced b-parameters
#' @param paramG  initial parameter data frame.
#'  This has priority over \code{init}.
#' @param method = 1 to use item category response function (icrf)
#'  to calculate rmse (default) \cr
#' = 2 to use item category information function (icif) to calculate rmse \cr
#' = 3 to use item response function (irf) to calculate rmse
#'  (not yet available), \cr
#' = 4 to use item information function (iif) to calculate rmse,
#'  (not yet available), \cr
#' @param wtype = 0 not to use dnorm(theta) as the weight
#' @param wmean The mean of normal distribution to be used as the weight
#' @param wsd The sd of normal distribution to be used as the weight
#' @param DinP = 1 to include D=1.7 in logistic function
#' @param npoints # of discrete points for theta
#' @param thmin Minimum value of discrete theta value
#' @param thmax Maximum value of discrete theta value
#' @param printGN print level for lazy.mat::GN function
#' @param maxiter Maximum # of GN iterations
#' @param eps  Convergence criterion for the relative improvement of rmse
#' @param epsg  Convergence crit for the maximum absolute value of the gradient
#' @param epsx  Convergence crit for the maximum absolute change of
#' the parameter value
#' @param print >= 1 to print result
#' @param plot >= 1 to plot result
#'
#' @return A list of: \cr
#' paramNew: Fitted normal GRM item parameter data frame\cr
#'  2PLM or 2PNM items remain unchaged. \cr
#' paramP: Input GPCM Item Parameter Data Frame\cr
#' grad: Gradient matrix \cr
#' wtype, wmean, wsd, method, init \cr
#' rmse_p: rmse in terms of icrf (minimized with method=1)\cr
#' rmse_irf: rmse in terms of irf\cr
#' rmse_ii: rmse in terms of item information\cr
#' rmse_iic: rmse in terms of item category information
#'  (minimized with method=2)\cr
#' icrfNew, icifNew, iifNew, irfNew \cr
#' icrfOld, icifOld, iifOld, irfOld
#'
#' @details
#' This function finds the set of Normal GRM item parameters
#' which best fit the given icrfs or item info functions of the items
#' in the input parameter data frame. \cr\cr
#'
#' If method = 0, this function minimizes\cr
#' \code{ sum( w*( vec(icrf(theta)) - vec(icrf_GRM(theta|PARAM)) )^2 ) } \cr
#' with respect to the GRM item parameters, PARAM, \cr
#' where \code{icrf(theta)} is the icrf of input items, \cr
#'  \code{icrf_GRM(theta|PARAM)} is the icrf of fitted GRM items, \cr
#' and \code{w} is the weight vector ( \code{N(wmean,wsd^2)} or 1 ). \cr
#'
#' If method = 1, this function minimizes\cr
#' \code{ sum( w*( (info_i(theta) - info_i_GRM(theta|PARAM) )^2 ) } \cr
#' with respect to the GRM item parameters, PARAM, \cr
#' where \code{info_i(theta)} is the item information function
#'  of input items and \cr
#'  \code{info_i_GRM(theta|PARAM)} is the item information function of
#'  the fitted GRM items. \cr
#'
#' If method = 2, this function minimizes\cr
#' \code{ sum( w*( vec(info_ic(theta)) - vec(info_ic_GRM(theta|PARAM)) )^2 ) }
#'  \cr
#' with respect to the GRM item parameters, PARAM, \cr
#' where \code{info_ic(theta)} is the item category information function
#'  of input items and \cr
#'  \code{info_ic_GRM(theta|PARAM)} is the item category information function
#'  of the fitted GRM items. \cr
#'
#' \cr
#' When three parameter binary items are included, two parameter normal ogive
#' model will be fitted.
#' \cr
#' Weighted Gauss-Newton method (\code{lazy.mat::GN}) is used for
#' the minimization with the numerical
#' Jacobian matrix calculated by \code{lazy.mat::JacobianMat}.
#' \cr
#'
#'
#' @examples
#' # convert 3PLM and GPCM items to normal GRM items
#' param <- paramA1[c(2,5,8),]
#' theta <- seq(-4,4,length=51)#'
#' # maxiter for method=2 is too small.
#' res1 <- conv2Gn( param, theta, maxiter=20, plot=1, wtype=1, method=1 )
#' res2 <- conv2Gn( param, theta, maxiter=2, plot=1, wtype=1, method=2 )
#'
#' Print(res1$rmse_p, res1$rmse_iic, res1$rmse_ii)
#' Print(res2$rmse_p, res1$rmse_iic, res1$rmse_ii)
#'
#'
#' @export
#'

conv2Gn <- function( param, theta=NULL, init=1, paramG=NULL
                    , method=0, wtype=1, wmean=0, wsd=1, DinP=1
                    , npoints=21, thmin=-3, thmax=3, printGN=0
                    , maxiter=500, eps=1e-6, epsg=1e-6, epsx=1e-9
                    , print=1, plot=0 ){
 # fitting Normal Graded Response Model to Generalized Partial Credit Model
 # Shin-ichi Mayekawa
 # fitG2P: 20080127-20180113
 # fitGn2P_ls: 20180113-20180114
 # fitGn2P_ls modified to fitGn2P_ls (old): 20180211
 # fitGn2P_ls modified: 20180213dnc
 # fitGn2P_ls modified to fitGn2P_ls: 20180213dnc
 # P -> icrfP: 20180213dnc
 # Ginv -> solve -> matSwp: 20180215
 # order restrictions on b: 20180215
 # info -> method, wtype=1: 20180218
 # init: 20180221
 # denominator of rmse corrected: 20180221
 # max abs value of pj: 20180221
 # GN: 20221205dnc,06
 # rmse divided by npoints: 20221206
 # paramG as arg: 20221206
 # use iif: 20221206,07
 # roxygen2: 20221207
 # return variables renamed: 20221208
 # G->Gn: 20221209
 # bugfix for init=1: 20221210
 # epsg, epsx: 20230202
 # return type B2 items: 20230213
 # renamed as conv2Gn: 20230220,21
 # method=0 -> method=1: 20230220
 # paramP -> param: 20230221
 # output rmse_irf: 20230222
 # roxygen2: 20231027
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
 #     paramP   =  Original parameter set with DinP
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


 info_ic <- function( paramnum ){
  # item category information functions* of one item (vectorized)
  paramjdf[1,4:(4+length(paramnum)-1)]=paramnum
  infovv=vec( iif( paramjdf, theta, print=0 )$info_item_cat )
  return( infovv )
 } # end of info_ic

 info_i <- function( paramnum ){
  # item information function of one item
  paramjdf[1,4:(4+length(paramnum)-1)]=paramnum

  dicrf=dirf( paramjdf, theta, print=0, plot=0 )$dICRF
  icrf=irf( paramjdf, theta, print=0, plot=0 )$ICRF
  infovv=rowSums( dicrf^2/icrf )

  # This is slower: Better modify iif.
  # infovv=vec( iif( paramjdf, theta, print=0 )$info_item )

  return( infovv )
 } # end of info_i



 ymyhat <- function( paramnum, method ){
  # resudual for GN
  if( method == 1 ){
    paramjdf[1,3+(1:ncatj)]=paramnum
    resid=vpj-c( irf( paramjdf, theta, print=0, plot=0 )$ICRF )
  } else if( method == 3 ){
    resid=vpj-info_i(paramnum)
  } else{
    resid=vpj-info_ic(paramnum)
  }
   return( resid )
 } # end of ymyhat


 # param
 paramP=param
 paramP0=paramP

 # argument name
 paramname=as.character(substitute(paramP))
 types=c("B","B3","Bn3","P","G")


 # pick up GPCM items
 isdf=0
 if( is.data.frame(paramP) ){
  # param
  isdf=1
  paramP=paramP[paramP[,"type"] %in% types,]
  if( nrow(paramP) == 0 ){
    cat("\n\nerror1:(conv2Gn) Input Item Parameter ERROR.\n")
    cat(" Input parameter data frame must have at least one of"
        , " the following types of items:\n")
    cat(" B, B3, Bn3, P, G\n")
    return()
  }
 }

 # Here, paramP contains only GPCM items ("B","B3","Bn3","P").

 locP=which(paramP$type=="P")
 locB3=which(paramP$type %in% c("B3","Bn3"))

 # const
 if( is.data.frame(paramP) ) iname=paramP$name
 else iname=rownames(paramP) # No chance.
 itype=paramP$type

 # generate theta and prior theta dist
 if( is.null(theta) ){
  theta=seq(thmin,thmax,length.out=npoints)
 }
 else{
  if( is.matrix(theta) ) theta=as.vector(theta)
  npoints=length(theta)
  thmin=theta[1]; thmax=theta[npoints]
 }
 thetac=format(theta,digits=3)

 # weights for those elements with which abs(v) > maxabsv will be
 # set equal to min(W)

 nitems=nrow(paramP)
 ncat=paramP[,"ncat"]
 ncat1=ncat-1

 if( print >= 1 ){
  cat("\nConversion to Normal Graded Response Model\n")
  cat("  LS criterion in terms of ")
  if( method == 1 )  cat( "item category response functions\n" )
  else if( method == 3 )  cat( "item information function\n" )
  else cat( "item category informations*\n" )
  cat("  # of items = ", nitems, ", weight type = ", wtype, sep="" )
  if( wtype == 1 ) cat(" with N(",wmean,",",wsd,")", sep="")
  cat("\n")
  cat("  # of theta points =", npoints, "in ["
      , thmin, " , ", thmax, "]\n")
  cat("  D in Partial Credit Model =", DinP,"\n")
  if( is.null(paramG) ) cat("  initial =", init, "\n")
 }

 # weight
 if( wtype == 0 ) w=rep(1,npoints)
 else if( wtype == 1 ){
  w=dnorm((theta-wmean)/wsd)
  w=w/sum(w)
 }
 w2=sqrt(w)
 wmin=min(w); wmax=max(w);

 # GPCM icrf matrix
 # res=icrfP( paramP, theta )
 res=irf( paramP, theta, print=0 )
 icrfP=res$ICRF
 irfP=res$IRF
 toP=res$toP; fromP=res$fromP
 rm(res)
 Pname=colnames(icrfP)

 # item category info: ic_infoP is npoints x ncat
 dicrf=dirf( paramP, theta, zero=1, print=0, plot=0 )$dICRF
 ic_infoP=iif( paramP, theta, print=0 )$info_item_cat
 iteminfoP=mapply( function(x,y) rowSums(ic_infoP[,x:y])
                      , fromP, toP, SIMPLIFY=FALSE )
 iteminfoP=matrix(unlist(iteminfoP),,nitems)
 colnames(iteminfoP)=iname


 # initial
 if(  is.null(paramG) ){
   if( init == 1 ){
     paramG=paramP
     paramG$type="Gn"
     GG=fitG2P( paramP[locP,], npoints=21, plot=0, print=0, wtype=wtype )
     if( !is.null(GG) & length(locP) > 0 ){
       GG$type="Gn"
       paramG[locP,1:ncol( GG)]=GG
     }
   }
   else{
     # equally spaced b-paramters in p2, p3, ..., p_ncat
     paramG=paramP
     for( j in 1:nrow(paramG) ){
       ncatj=paramG[j,"ncat"]
       paramG[j,4+(1:(ncatj-1))]=1:(ncatj-1)-ncatj/2
     }
     paramG$type="Gn"
   }
   if( length(locB3) > 0 ){
     paramG[locB3,]$p3=0
     paramG[locB3,]$type="Bn"
   }
 } # end of initial

 # Here, paramG is the initial parameter data frame.

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

   # initial
   paramjdf=paramG[j,]
   paramj=paramjdf[1,3+(1:ncatj)]

   #
   if( method == 1 ){
    Pj=icrfP[,fromP[j]:toP[j]]
   } else if( method == 2 ){
    Pj=ic_infoP[,fromP[j]:toP[j]]    # This Pj is Item Category Info
   } else if( method == 3 ){
    Pj=rowSums(ic_infoP[,fromP[j]:toP[j]])  # This is Item Info
   }
   vpj=vec(Pj)
   if( method != 3 ){
    ww=rep(w,ncat[j])
   }
   else{
    ww=w
   }

   rss=sum(ww*ymyhat(paramj,method=method)^2)
   rmsejp=sqrt(rss/npoints)
   if( method != 3 ) rmsejp=sqrt(rss/npoints/ncat[j])

   res=GN( paramj, ymyhat, ww, method=method, print=printGN
           , maxiter=maxiter, eps=eps, epsg=epsg, epsx=epsx )
   paramj=res$par
   llll=res$iterations
   g=res$grad
   rss=res$objective
   rmsej=sqrt(rss/npoints)
   if( method != 3 ) rmsej=sqrt(rss/npoints/ncat[j])

   if( print >= 1 ){
    cat("\n Iteration terminated:\n")
    Print( llll, rmsej, rmsejp, g, fmt="i3 10.6")
   }

   # truncate
   paramj=trunc(paramj*10000)/10000

   paramG[j,3+(1:ncatj)]=paramj
   rmse[j]=rmsej
   grad[j,1:ncat[j]]=g

  } # end of j loop for item



 # fitted GRM icrf
 res=irf( paramG, theta, print=0 )
 icrfG=res$ICRF
 irfG=res$IRF
 ic_infoG=iif( paramG, theta, print=0 )$info_item_cat
 iteminfoG=mapply( function(x,y) rowSums(ic_infoG[,x:y])
                  , fromP, toP, SIMPLIFY=FALSE )
 iteminfoG=matrix(unlist(iteminfoG),,nitems)
 colnames(iteminfoG)=iname

 # rmse in terms of P
 ic_rss=colSums( w*(icrfP-icrfG)^2 )/npoints
 rmse_p=sqrt(
  unlist( mapply( function(x,y){sum(ic_rss[x:y])/(y-x+1)}
                     , fromP, toP, SIMPLIFY=FALSE ) ) )
 rmse_pF=formatC(rmse_p, digits=5, format="f")

 # rmse in terms of item category info
 ic_rss=colSums( w*(ic_infoP-ic_infoG)^2 )/npoints
 rmse_iic=sqrt(
  unlist( mapply( function(x,y){sum(ic_rss[x:y])/(y-x+1)}
                  , fromP, toP, SIMPLIFY=FALSE ) ) )
 rmse_iicF=formatC(rmse_iic, digits=5, format="f")

 # rmse in terms of irf
 rmse_irf=sqrt( colSums( w*(irfG-irfP)^2/npoints ) )
 rmse_irfF=formatC(rmse_irf, digits=5, format="f")
 
 # rmse in terms item info
 rmse_ii=sqrt( colSums( w*(iteminfoP-iteminfoG)^2/npoints ) )
 rmse_iiF=formatC(rmse_ii, digits=5, format="f")


 names(rmse_p)=iname
 names(rmse_irf)=iname
 names(rmse_iic)=iname
 names(rmse_ii)=iname



 if( print ){
  cat("\nInput and  Converted Item Parameters\n")
  cat("\n Input GPCM Item Parameters\n")
  Print(paramP)
  cat("\n Converted Normal GRM\n")
  Print(paramG)
  cat("Gradients and RMSE minimized:  method =", method, "\n")
  cat(" Depending on the value of method (1,2), resp,"
  , " rmse_p or rmse_iic was minimized.\n ")
  Print( rmse_p, rmse_irf)
  Print( rmse_iic, rmse_ii )
  Print(grad, rmse, fmt="8.5")
  if( print >= 2 ){
   cat("\n ICRFs of the original GPCM\n")
   Print(thetac, icrfP, fmt="5,3")
   cat(" ICRFs of the converted GRM\n")
   Print(thetac, icrfG, fmt="5,3")
  }
  if( print >= 3 ){
    cat("\n Item Information Function of the original GPCM\n")
    Print(thetac, iteminfoP, fmt="5,3")
    cat(" Item Information Function of the converted GRM\n")
    Print(thetac, iteminfoG, fmt="5,3")
  }
  if( print >= 4 ){
   cat("\n Item Category Informations of the original GPCM\n")
   Print(thetac, ic_infoP, fmt="5,3")
   cat(" Item Category Informations of the converted GRM\n")
   Print(thetac, ic_infoG, fmt="5,3")
  }
 }


 if( plot ){

  nparams=paramP$ncat
  nparams[locB3]=nparams[locB3]+1


  # icrf
  for( j in 1:(nitems) ){
   title=paste("ICRF of ",iname[j]," (type=", itype[j]," ncat=",ncat[j]
               , ")  rmse_p=", rmse_pF[j], " (method=", method
               , ", wtype=", wtype, ")", sep="")
   sub=paste( "paramP ="
              , paste(format(paramP[j,3+(1:nparams[j])],digits=3)
                      , collapse=",  ")
              , "\nparamGn ="
              , paste(format(paramG[j,3+(1:ncat[j])],digits=3)
                      , collapse=",  ")
   )
   for( k in 1:(ncat[j]) ){
    plot(theta,icrfP[,fromP[j]+k-1], xlab="", col=1, lty=1, lwd=2
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="")
    par(new=1)
   }
   for( k in 1:(ncat[j]-1) ){
    plot(theta,icrfG[,fromP[j]+k-1], xlab="", col=2, lty=1, lwd=2
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="b", ylab="")
    par(new=1)
   }
   plot(theta,icrfG[,toP[j]], main=title, sub=sub
        , xlab="", col=2, lty=1, lwd=2
        , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="b", ylab="ICRF")
   legend( range(theta)[1]+0.5, 0.9, c("P","G"), cex=0.8, col=1:2
           ,  pch=c(NA,21), lty=c(1,2), title="model" )
   par(new=0)
  }

  # irf
  maxy=max(c(irfG,irfP))
  for( j in 1:(nitems) ){
    ii=irfG[,j]; iih=irfP[,j]
    title=paste("Item Response Function of ",iname[j]
                ," (type=", itype[j]," ncat=",ncat[j]
                , ")  rmse_i=", rmse_iiF[j], " (method=", method
                , ", wtype=", wtype, ")", sep="")
    sub=paste( "paramG ="
               , paste(format(paramG[j,3+(1:nparams[j])],digits=3)
                       , collapse=",  ")
               , "\nparamP ="
               , paste(format(paramP[j,3+(1:ncat[j])],digits=3)
                       , collapse=",  ")
    )
    matplot( theta, cbind(ii,iih), main=title, sub=sub,  pch=c(NA,21), xlab=""
             , type="l", ylab="information", ylim=c(0,maxy), lty=c(1,2), lwd=2 )
    legend( range(theta)[1]+0.5, 0.9*maxy, c("G","P"), cex=1, col=1:2
            ,  pch=c(NA,NA), lty=c(1,1), title="model" )
  }
  
  # item info
  maxy=max(c(iteminfoP,iteminfoG))
  for( j in 1:(nitems) ){
   ii=rowSums( ic_infoP[,fromP[j]:toP[j]] )
   iih=rowSums( ic_infoG[,fromP[j]:toP[j]] )
   title=paste("Item Information Function of ",iname[j]
               ," (type=", itype[j]," ncat=",ncat[j]
               , ")  rmse_i=", rmse_iiF[j], " (method=", method
               , ", wtype=", wtype, ")", sep="")
   sub=paste( "paramP ="
              , paste(format(paramP[j,3+(1:nparams[j])],digits=3)
                      , collapse=",  ")
              , "\nparamGn ="
              , paste(format(paramG[j,3+(1:ncat[j])],digits=3)
                      , collapse=",  ")
   )
   matplot( theta, cbind(ii,iih), main=title, sub=sub,  pch=c(NA,21), xlab=""
          , type="l", ylab="information", ylim=c(0,maxy), lty=c(1,2), lwd=2 )
   legend( range(theta)[1]+0.5, 0.9*maxy, c("P","G"), cex=1, col=1:2
           ,  pch=c(NA,NA), lty=c(1,1), title="model" )
  }

  if( plot >= 3 ){
   # cat info
   maxy=max(c(ic_infoP),c(ic_infoG))
   for( j in 1:(nitems) ){
    title=paste("Item Category Information of "
                ,iname[j]," (type=", itype[j]," ncat=",ncat[j]
                , ")  rmse_ic=", rmse_iicF[j], " (method=", method
                , ", wtype=", wtype, ")", sep="")
    sub=paste( "paramP ="
               , paste(format(paramP[j,3+(1:nparams[j])],digits=3)
                       , collapse=",  ")
               , "\nparamGn ="
               , paste(format(paramG[j,3+(1:ncat[j])],digits=3)
                       , collapse=",  ")
    )
    for( k in 1:(ncat[j]) ){
     plot(theta,ic_infoP[,fromP[j]+k-1], xlab="", col=1, lty=1, lwd=2, pch=k
          , xlim=c(min(theta),max(theta)), ylim=c(0,maxy), type="b", ylab="")
     par(new=1)
    }
    for( k in 1:(ncat[j]-1) ){
     plot(theta,ic_infoG[,fromP[j]+k-1], xlab="", col=2, lty=1, lwd=2, pch=k
          , xlim=c(min(theta),max(theta)), ylim=c(0,maxy), type="b", ylab="")
     par(new=1)
    }
    plot(theta,ic_infoG[,toP[j]], main=title, sub=sub, xlab="", type="b"
         , col=2, lty=1, lwd=2, pch=ncat[j]
         , xlim=c(min(theta),max(theta)), ylim=c(0,maxy), ylab="ICRF")
    legend( range(theta)[1]+0.5, 0.9*maxy, c("P","G"), cex=0.8, col=1:2
            ,  pch=c(NA,21), lty=c(1,1), title="model" )
    par(new=0)
   } # end of j
  } # end of >= 3

 } # end of plot


 paramP0[paramP0$type %in% types,]=paramG
 paramG=paramP0


 return( named_list( paramNew=paramG, paramOld=param, method
                    , rmse_icrf=rmse_p, rmse_icif=rmse_iic, rmse_iif=rmse_ii
                    , rmse_p, rmse_iic, rmse_irf, rmse_ii
                    , grad, wtype, wmean, wsd, init
                    , theta, w, fromP, toP
                    , icrfNew=icrfG, icifNew=ic_infoG, iifNew=iteminfoG
                    , icrfOld=icrfP, icifOld=ic_infoP, iifOld=iteminfoP
                    , irfNew=irfP, irfOld=irfG ) )

} # end of conv2Gn


