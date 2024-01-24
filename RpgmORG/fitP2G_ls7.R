#' Conversion of Graded Response Items to Generalized Partial Credit Items
#'
#' @param paramG Item Parameter Data Frame for 3PLM or nomal orgive or GRM items
#' @param theta Vector of theta points
#' @param init  = 1 to use fitP2G, else use equally spaced b-parameters
#' @param paramP  initial parameter data frame.
#'  This has priority over \code{init}.
#' @param method = 0 to use icrf to calculate rmse (default) \cr
#' = 1 to use item info to calculate rmse, \cr
#' = 2 to use item category info* to calculate rmse
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
#' @param print >= 1 to print result
#' @param plot >= 1 to plot result
#'
#' @return A list of: \cr
#' paramP: Fitted GPCM item parameter data framer
#' paramG: Input GRM Item Parameter Data Frame (type="Gn")
#' grad: Gradient matrix \cr
#' wtype, wmean, wsd, method, init \cr
#' rmse_p: rmse in terms of icrf (method=0)\cr
#' rmse_ii: rmse in terms of item information (method=1)\cr
#' rmse_iic: rmse in terms of item category information* (method=2)\cr
#' icrfNew, icifNew, iifNew \cr
#' icrfOld, icifOld, iifOld
#'
#' @details
#' This function finds the set of GPCM item paramters
#' which best fit the given icrfs or item info functions of the items
#' in the input parameter data frame. \cr\cr
#'
#' If method = 0, this function minimizes\cr
#' \code{ sum( w*( vec(icrf(theta)) - vec(icrf_GPCM(theta|PARAM)) )^2 ) } \cr
#' with respect to the GPCM item parameters, PARAM, \cr
#' where \code{icrf(theta)} is the icrf of the input items, \cr
#'  \code{icrf_GPCM(theta|PARAM)} is the icrf of the fitted GPCM items, \cr
#' and \code{w} is the weight vector ( \code{N(wmean,wsd^2)} or 1 ). \cr
#'
#' If method = 1, this function minimizes\cr
#' \code{ sum( w*( info_i(theta) - info_i_GPCM(theta|PARAM) )^2 ) } \cr
#' with respect to the GPCM item parameters, PARAM, \cr
#' where \code{info_i(theta)} is the item information function
#'  of input items and \cr
#'  \code{info_i_GPCM(theta|PARAM)} is the item information function of
#'  the fitted GPCM items. \cr
#'
#' If method = 2, this function minimizes\cr
#' \code{ sum( w*( vec(info_ic(theta)) - vec(info_ic_GPCM(theta|PARAM)) )^2 ) }
#'  \cr
#' with respect to the GPCM item parameters, PARAM, \cr
#' where \code{info_ic(theta)} is the item category information function
#'  of input items and \cr
#'  \code{info_ic_GPCM(theta|PARAM)} is the item category information function
#'  of the fitted GPCM items. \cr
#'
#' \cr
#' When three parameter binary items are included, 2PLM will be fitted.
#' \cr
#' Weighted Gauss-Newton method (\code{lazy.mat::GN}) is used for
#' the minimization with the numerical
#' Jacobian matrix calculated by \code{lazy.mat::JacobianMat}.
#' \cr
#'
#'
#' @examples
#' paramP1 <- fitP2G_ls( paramS2, plot=1, print=1 )$paramNew
#' paramG1 <- fitGn2P_ls( paramP1, plot=1, print=1 )
#'
#' # convert 3PLM and GPCM items
#' param <- paramA1[c(2,5,8),]
#' theta <- seq(-4,4,length=51)
#'
#' # maxiter below is too small!!
#' res0 <- fitP2G_ls( param, theta, maxiter=20, plot=1, wtype=1, method=0 )
#' res1 <- fitP2G_ls( param, theta, maxiter=20, plot=1, wtype=1, method=1 )
#' res2 <- fitP2G_ls( param, theta, maxiter=3, plot=1, wtype=1, method=2 )
#'
#' Print(res0$rmse_p, res0$rmse_iic, res0$rmse_ii)
#' Print(res1$rmse_p, res1$rmse_iic, res1$rmse_ii)
#' Print(res2$rmse_p, res2$rmse_iic, res2$rmse_ii)
#'
#'
#' @export
#'

fitP2G_ls <- function( paramG, theta=NULL, init=1, paramP=NULL
                    , method=0, wtype=1, wmean=0, wsd=1, DinP=1
                    , npoints=21, thmin=-3, thmax=3, printGN=0
                    , maxiter=500, eps=1e-6, print=1, plot=0 ){
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
 # use iif: 20221206
 # converted from fitG2P_ls: 20221206
 # roxygen2: 20221207
 # plot caption bug fix: 20221208
 # return variables renamed: 20221208
 # bugfix for init=1: 20221210
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
  # matplot( theta, dicrf^2/icrf, type="l")
  infovv=rowSums( dicrf^2/icrf )
  return( infovv )
 } # end of info_i



 ymyhat <- function( paramnum, method ){
  # resudual for GN
  if( method == 0 ){
    paramjdf[1,3+(1:ncatj)]=paramnum
    resid=vpj-c( irf( paramjdf, theta, print=0, plot=0 )$ICRF )
  } else if( method == 1 ){
    resid=vpj-info_i(paramnum)
  } else{
    resid=vpj-info_ic(paramnum)
  }
   return( resid )
 } # end of ymyhat



 # argument name
 paramname=as.character(substitute(paramG))


 # pick up GRM items
 isdf=0
 if( is.data.frame(paramG) ){
  # param
  isdf=1
  paramG=paramG[paramG[,"type"] %in% c("Bn","B3","Bn3","G","Gn"),]
  if( nrow(paramG) == 0 ){
   cat("\n\nerror1:(fitP2G_ls) Input Item Parameter ERROR.\n")
    cat(" Input parameter data frame must have at least one of"
        , " the following types of items:\n")
    cat(" B, B3, Bn, Bn3, G, Gn\n")
   return()
  }
 }

 # Here, paramG contains only GRM items ("Bn","B3","Bn3","G","Gn").

 locG=which(paramG$type %in% c("G","Gn"))
 locB3=which(paramG$type %in% c("B3","Bn3"))

 # const
 if( is.data.frame(paramG) ) iname=paramP$name
 else iname=rownames(paramG) # No chance.
 itype=paramG$type

 # generate thata and prior theta dist
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

 nitems=nrow(paramG)
 ncat=paramG[,"ncat"]
 ncat1=ncat-1

 if( print >= 1 ){
  cat("\nFitting Generalized Partial Credit Model"
      ,"to Graded Response Model\n")
  cat("  LS criterion in terms of ")
  if( method == 0 )  cat( "item category response functions\n" )
  else if( method == 1 )  cat( "item information function\n" )
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
 res=irf( paramG, theta, print=0 )
 icrfG=res$ICRF
 toP=res$toP; fromP=res$fromP
 rm(res)
 Pname=colnames(icrfG)

 # item category info: ic_infoP is npoints x ncat
 dicrf=dirf( paramG, theta, zero=1, print=0, plot=0 )$dICRF
 icrf=irf( paramG, theta, zero=1, print=0, plot=0 )$ICRF
 ic_infoG=iif( paramG, theta, print=0 )$info_item_cat
 iteminfoG=mapply( function(x,y) rowSums(ic_infoG[,x:y])
                      , fromP, toP, SIMPLIFY=FALSE )
 iteminfoG=matrix(unlist(iteminfoG),,nitems)
 colnames(iteminfoG)=iname


 # initial
 if(  is.null(paramP) ){
   if( init == 1 ){
     paramP=paramG
     paramP$type="P"
     GG=fitP2G( paramG[locG,], npoints=21, plot=0, print=0, wtype=wtype )
     GG=GG$paramP
     if( !is.null(GG) &  length(locG) > 0 ){
       GG$type="P"
       paramP[locG,1:ncol( GG)]=GG
     }
   }
   else{
     # equally spaced b-paramters in p2, p3, ..., p_ncat
     paramP=paramG
     for( j in 1:nrow(paramP) ){
       ncatj=paramP[j,"ncat"]
       paramP[j,4+(1:(ncatj-1))]=1:(ncatj-1)-ncatj/2
     }
     paramP$type="P"
   }
   if( length(locB3) > 0 ){
     paramP[locB3,]$p3=0
     paramP[locB3,]$type="B"
   }
 } # end of initial

 # Here, paramG is the initial parameter data frame.

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

   # initial
   paramjdf=paramP[j,]
   paramj=paramjdf[1,3+(1:ncatj)]

   #
   if( method == 0 ){
    Pj=icrfG[,fromP[j]:toP[j]]
   } else if( method == 2 ){
    Pj=ic_infoG[,fromP[j]:toP[j]]    # This Pj is Item Category Info
   } else if( method == 1 ){
    Pj=rowSums(ic_infoG[,fromP[j]:toP[j]])  # This is Item Info
   }
   vpj=vec(Pj)
   if( method == 2 ){
    ww=rep(w,ncat[j])
   }
   else{
    ww=w
   }

   rss=sum(ww*ymyhat(paramj,method=method)^2)
   rmsejp=sqrt(rss/npoints)
   if( method != 1 ) rmsejp=sqrt(rss/npoints/ncat[j])
   res=GN( paramj, ymyhat, ww, method=method, print=printGN
           , maxiter=maxiter, eps=eps )
   paramj=res$par
   llll=res$iterations
   g=res$grad
   rss=res$objective
   rmsej=sqrt(rss/npoints)
   if( method != 1 ) rmsej=sqrt(rss/npoints/ncat[j])

   if( print >= 1 ){
    cat("\n Iteration terminated:\n")
    Print( llll, rmsej, rmsejp, g, fmt="i3 10.6")
   }

   # truncate
   paramj=trunc(paramj*10000)/10000

   paramP[j,3+(1:ncatj)]=paramj
   rmse[j]=rmsej
   grad[j,1:ncat[j]]=g

  } # end of j loop for item



 # fitted GPCM icrf
 icrfP=irf( paramP, theta, print=0 )$ICRF
 ic_infoP=iif( paramP, theta, print=0 )$info_item_cat
 iteminfoP=mapply( function(x,y) rowSums(ic_infoP[,x:y])
                  , fromP, toP, SIMPLIFY=FALSE )
 iteminfoP=matrix(unlist(iteminfoP),,nitems)
 colnames(iteminfoP)=iname

 # rmse in terms of P
 ic_rss=colSums( w*(icrfG-icrfP)^2 )/npoints
 rmse_p=sqrt(
  unlist( mapply( function(x,y){sum(ic_rss[x:y])/(y-x+1)}
                     , fromP, toP, SIMPLIFY=FALSE ) ) )
 rmse_pF=formatC(rmse_p, digits=5, format="f")

 # rmse in terms of item category info
 ic_rss=colSums( w*(ic_infoG-ic_infoP)^2 )/npoints
 rmse_iic=sqrt(
  unlist( mapply( function(x,y){sum(ic_rss[x:y])/(y-x+1)}
                  , fromP, toP, SIMPLIFY=FALSE ) ) )
 rmse_iicF=formatC(rmse_iic, digits=5, format="f")

 # rmse in terms item info
 rmse_ii=sqrt( colSums( w*(iteminfoG-iteminfoP)^2/npoints ) )
 rmse_iiF=formatC(rmse_ii, digits=5, format="f")


 names(rmse_p)=iname
 names(rmse_iic)=iname
 names(rmse_ii)=iname

 Print(rmse_p, rmse_iic,rmse_ii)


 if( print ){
   cat("\nInput and  Converted Item Parameters\n")
   cat("\n Input GRM Item Parameters\n")
   Print(paramG)
   cat("\n Converted GPCM\n")
   Print(paramP)
   cat("Gradients and RMSE minimized:  method =", method, "\n")
   cat(" Depending on the value of method (0,1,2), resp,"
       , " rmse_p, rmse_ii, or rmse_iic was minimized.\n ")
   Print( rmse_p, rmse_ii, rmse_iic )
   Print(grad, rmse, fmt="8.5")
   if( print >= 2 ){
     cat("\n ICRFs of the original GRM\n")
     Print(thetac, icrfG, fmt="5,3")
     cat(" ICRFs of the converted GPCM\n")
     Print(thetac, icrfP, fmt="5,3")
   }
   if( print >= 3 ){
     cat("\n Item Information Function of the original GRM\n")
     Print(thetac, iteminfoG, fmt="5,3")
     cat(" Item Information Function of the converted GPCM\n")
     Print(thetac, iteminfoP, fmt="5,3")
   }
   if( print >= 4 ){
     cat("\n Item Category Informations of the original GRM\n")
     Print(thetac, ic_infoG, fmt="5,3")
     cat(" Item Category Informations of the converted GPCM\n")
     Print(thetac, ic_infoP, fmt="5,3")
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
   sub=paste( "paramG ="
              , paste(format(paramG[j,3+(1:nparams[j])],digits=3)
                      , collapse=",  ")
              , "\nparamP ="
              , paste(format(paramP[j,3+(1:ncat[j])],digits=3)
                      , collapse=",  ")
   )
   for( k in 1:(ncat[j]) ){
    plot(theta,icrfG[,fromP[j]+k-1], xlab="", col=1, lty=1, lwd=2
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="")
    par(new=1)
   }
   for( k in 1:(ncat[j]-1) ){
    plot(theta,icrfP[,fromP[j]+k-1], xlab="", col=2, lty=1, lwd=2
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="b", ylab="")
    par(new=1)
   }
   plot(theta,icrfP[,toP[j]], main=title, sub=sub
        , xlab="", col=2, lty=1, lwd=2
        , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="b", ylab="ICRF")
   legend( range(theta)[1]+0.5, 0.9, c("G","P"), cex=0.8, col=1:2
           ,  pch=c(NA,21), lty=c(1,2), title="model" )
   par(new=0)
  }

  # item info
  maxy=max(c(iteminfoG,iteminfoP))
  for( j in 1:(nitems) ){
   ii=rowSums( ic_infoG[,fromP[j]:toP[j]] )
   iih=rowSums( ic_infoP[,fromP[j]:toP[j]] )
   title=paste("Item Information Function of ",iname[j]
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

  if( plot >= 3 ){
   # cat info
   maxy=max(c(ic_infoG),c(ic_infoP))
   for( j in 1:(nitems) ){
    title=paste("Item Category Information of "
                ,iname[j]," (type=", itype[j]," ncat=",ncat[j]
                , ")  rmse_ic=", rmse_iicF[j], " (method=", method
                , ", wtype=", wtype, ")", sep="")
    sub=paste( "paramG ="
               , paste(format(paramG[j,3+(1:nparams[j])],digits=3)
                       , collapse=",  ")
               , "\nparamP ="
               , paste(format(paramP[j,3+(1:ncat[j])],digits=3)
                       , collapse=",  ")
    )
    for( k in 1:(ncat[j]) ){
     plot(theta,ic_infoG[,fromP[j]+k-1], xlab="", col=1, lty=1, lwd=2, pch=k
          , xlim=c(min(theta),max(theta)), ylim=c(0,maxy), type="b", ylab="")
     par(new=1)
    }
    for( k in 1:(ncat[j]-1) ){
     plot(theta,ic_infoP[,fromP[j]+k-1], xlab="", col=2, lty=1, lwd=2, pch=k
          , xlim=c(min(theta),max(theta)), ylim=c(0,maxy), type="b", ylab="")
     par(new=1)
    }
    plot(theta,ic_infoP[,toP[j]], main=title, sub=sub, xlab="", type="b"
         , col=2, lty=1, lwd=2, pch=ncat[j]
         , xlim=c(min(theta),max(theta)), ylim=c(0,maxy), ylab="ICRF")
    legend( range(theta)[1]+0.5, 0.9*maxy, c("G","P"), cex=0.8, col=1:2
            ,  pch=c(NA,21), lty=c(1,1), title="model" )
    par(new=0)
   } # end of j
  } # end of >= 3

 } # end of plot



 return( named_list( paramNew=paramP, paramOld=paramG, rmse_p, method
                    , rmse_ii, rmse_iic
                    , grad, wtype, wmean, wsd, method, init
                    , theta, w, fromP, toP
                    , icrfNew=icrfP, icifNew=ic_infoP, iifNew=iteminfoP
                    , icrfOld=icrfG, icifOld=ic_infoG, iifOld=iteminfoG ) )

} # end of fitP2G_ls






comments(
'




res2=fitP2G_ls( paramA1[9,], plot=1, print=9, wtype=1, method=0, init=0 )


itemparam=paramA1[c(1,2,4),]

res0=fitP2G_ls( itemparam, theta, plot=4, print=1, wtype=1, method=0, init=0 )
rmse_p0=res0$rmse_p; rmse_ii0=res0$rmse_ii; rmse_iic0=res0$rmse_iic

res1=fitP2G_ls( itemparam, theta, plot=4, print=1, wtype=1, method=1, init=0 )
rmse_p1=res1$rmse_p; rmse_ii1=res1$rmse_ii; rmse_iic1=res1$rmse_iic

res2=fitP2G_ls( itemparam, theta, plot=4, print=1, wtype=1, method=2, init=0 )
rmse_p2=res2$rmse_p; rmse_ii2=res2$rmse_ii; rmse_iic2=res2$rmse_iic

rmse_p=rbind(rmse_p0,rmse_p1,rmse_p2)
rmse_ii=rbind(rmse_ii0,rmse_ii1,rmse_ii2)
rmse_iic=rbind(rmse_iic0,rmse_iic1,rmse_iic2)
rownames(rmse_p)=0:2; rownames(rmse_ii)=0:2; rownames(rmse_iic)=0:2

Print(rmse_p,rmse_ii,rmse_iic, fmt=".4")




parambad=cbind( data.frame( name="", type="P", stringsAsFactors=0 )
                ,statbad[,3:7] )
cn=colnames(parambad)
cn=gsub("p1","p", cn)
colnames(parambad)=cn
parambad$name=rownames(statbad)

res=fitP2G_ls( parambad, theta, plot=1, print=1, wtype=1, method=2 )

# res=fitG2P_ls( parambad[1,], theta, plot=1, print=3, wtype=1, method=1 )





param=paramA1[c(2,5,8),]
param=paramA1[2,]
param=paramA1[c(5,8),]

param=paramA1[c(2,5,8),]

theta=seq(-4,4,length=51)

# res=fitP2G_ls( param, theta, plot=1, print=2, maxiter=50, wtype=1, method=0 )
res=fitP2G_ls( param, theta, plot=1, print=1, wtype=1, method=1 )








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
