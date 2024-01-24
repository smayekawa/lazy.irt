#' Conversion of Partial Credit Items to Normal Graded Response Items
#'
#' @param paramP Item Parameter Data Frame for 2PLM or GPCM Items
#' @param theta Vector of theta points
#' @param init  = 1 to use fitP2G, else use equally spaced b-parameters
#' @param method = 0 to use icrf to calculate rmse (default) \cr
#' = 1 to use item info to calculate rmse, \cr
#' = 2 to use item category info* to calculate rmse
#' @param wtype = 0 not to use dnorm(theta) as the weight
#' @param wmean The mean of normal distribution to be used as the weight
#' @param wsd The sd of normal distribution to be used as the weight
#' @param DinP = 1 to include D=1.7 in logistic function
#' @param npoints # of discrete points for theta
#' @param thmin Minimum value of discrete thata value
#' @param thmax Maximum value of discrete thata value
#' @param maxiter Maximum # of GN iterations
#' @param eps  Convergence criterion for the relative improvement of rmse
#' @param print >= 1 to print result
#' @param plot >= 1 to plot result
#'
#' @return A list of: \cr
#' paramP: Input GPCM item parameter data frame (subset, type="P")\cr
#' paramG: GRM Item Parameter Data Frame (type="Gn")\cr
#' grad: Gradient matrix \cr
#' wtype, wmean, wsd, method, init \cr
#' rmse_p: rmse in terms of icrf (method=0)\cr
#' rmse_ii: rmse in terms of item infomation (method=1)\cr
#' rmse_iic: rmse in terms of item category infomation* (method=2)\cr
#' icrfP, infoP, iteminfoP \cr
#' icrfG, infoG, iteminfoG
#'
#' @details
#' This function finds the set of Normal GRM item paramters
#' which best fit the given icrfs or item info functions of the items
#' in the imput parameter data frame. \cr\cr
#'
#' If method = 0, this function minimizes\cr
#' \code{ sum( w*( vec(icrf_P(theta)) - vec(icrf_GRM(theta)) )^2 ) } \cr
#' with respect to the GRM item parameters, \cr
#' where \code{icrf_P(theta)} is the icrf of input GPCM items, \cr
#'  \code{icrf_GRM(theta)} is the icrf of fitted GRM items, \cr
#' and \code{w} is the weight vector ( N(wmean,wsd^2) or 1 ). \cr\cr
#' If method = 1, this function minimizes\cr
#' \code{ sum( w*( vec(info_i(theta)) - vec(info_i_GRM(theta)) )^2 ) } \cr
#' with respect to the GRM item parameters, \cr
#' where \code{info_i(theta)} is the item infomation function
#'  of input GPCM items and \cr
#'  \code{info_i_GRM(theta)} is the item infomation function of
#'  fitted GRM items. \cr
#' If method = 2, this function minimizes\cr
#' \code{ sum( w*( vec(info_ic(theta)) - vec(info_ic_GRM(theta)) )^2 ) } \cr
#' with respect to the GRM item parameters, \cr
#' where \code{info_ic(theta)} is the item category infomation function
#'  of input GPCM items and \cr
#'  \code{info_ic_GRM(theta)} is the item category infomation function* of
#'  fitted GRM items. \cr\cr
#'
#' \cr
#' When 3PLM items are included, 2PLM will be fitted.
#' \cr
#'
#' \cr
#' Weighted Gauss-Newton method is used for the minimization.
#' \cr
#'
#'
#' @examples
#' paramP1 <- fitP2G_ls( paramS2, plot=1, print=1 )$paramP
#' paramG1 <- fitGn2P_ls( paramP1, plot=1, print=1 )
#'
#' # convert 3PLM and GPCM items
#' param <- paramA1[c(2,5,8),]
#' theta <- seq(-4,4,length=51)
#'
#' # maxiter below is too small!!
#' res0 <- fitGn2P_ls( param, theta, maxiter=20, plot=1, wtype=1, method=0 )
#' res1 <- fitGn2P_ls( param, theta, maxiter=20, plot=1, wtype=1, method=1 )
#'
#' Print(res0$rmse_p, res0$rmse_iic, res0$rmse_ii)
#' Print(res1$rmse_p, res1$rmse_iic, res1$rmse_ii)
#'
#'
#' @export
#'

fitGn2P_ls <- function( paramP, theta=NULL, init=1
                    , method=0, wtype=1, wmean=0, wsd=1, DinP=1
                    , npoints=21, thmin=-3, thmax=3
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
  dicrf=dirf( paramjdf, theta, print=0, plot=0 )$dICRF
  icrf=irf( paramjdf, theta, print=0, plot=0 )$ICRF
  # matplot( theta, dicrf^2/icrf, type="l")
  infovv=vec( dicrf^2/icrf )
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




 # argument name
 paramname=as.character(substitute(paramP))

 # omit type "P" and ncat=2 items
 loc=which( paramP$type=="P" & paramP$ncat==2 )
 loc=which( paramP$type=="P" & paramP$ncat==200 )
 if( length(loc) > 0 ) paramP=paramP[-loc,]

 # pick up GPCM items
 isdf=0
 if( is.data.frame(paramP) ){
  # param
  isdf=1
  paramP=paramP[paramP[,"type"] %in% c("B","B3","Bn3","P"),]
  if( nrow(paramP) == 0 ){
   cat("\n\nerror1:(fitG2P_ls) Input Item Parameter ERROR.\n\n")
   return()
  }
 }
 locP=which(paramP$type=="P")
 locB3=which(paramP$type %in% c("B3","Bn3"))

 # const
 if( is.data.frame(paramP) ) iname=paramP$name
 else iname=rownames(paramP)

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

 nitems=nrow(paramP)
 ncat=paramP[,"ncat"]
 ncat1=ncat-1

 if( print >= 1 ){
  cat("\nFitting Normal Graded Response Model"
      ,"to Generalized Partial Credit Model\n")
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
  cat("  initial =", init, "\n")
 }

 # weight
 if( wtype == 0 ) w=rep(1,npoints)
 else if( wtype == 1 ){
  w=dnorm((theta-wmean)/wsd)
  w=w/sum(w)
 }
 w2=sqrt(w)
 wmin=min(w); wmax=max(w);

 # GPCM icrf matrix:  This is Y
 # res=icrfP( paramP, theta )
 res=irf( paramP, theta, print=0 )
 icrfP=res$ICRF
 toP=res$toP; fromP=res$fromP
 rm(res)
 Pname=colnames(icrfP)

 # item category info: infoP is npoints x ncat
 dicrf=dirf( paramP, theta, zero=1, print=0, plot=0 )$dICRF
 icrf=irf( paramP, theta, zero=1, print=0, plot=0 )$ICRF
 infoP=dicrf^2/icrf
 iteminfoP=mapply( function(x,y) rowSums(infoP[,x:y])
                      , fromP, toP, SIMPLIFY=FALSE )
 iteminfoP=matrix(unlist(iteminfoP),,nitems)
 colnames(iteminfoP)=iname

 # initial
 if( init == 1 ){
  paramG=paramP
  GG=fitG2P( paramP[locP,], npoints=21, plot=0, print=0, wtype=wtype )
  if( length(locP) > 0 ){
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
 rmse=numeric(nitems)
 names(rmse)=iname
 grad=matrix(NA,max(ncat),nitems)
 rownames(grad)=paste("p",1:nrow(grad),sep="")
 colnames(grad)=iname
 grad=t(grad)


 if( method == 0 ){
  ################### in terms of icrf ##################################


  # for each GPCM item
  for( j in 1: nitems ){

   if( print >= 1 ) cat("\nItem # ", j, "\n", sep="")

   ncatj=ncat[j]
   Pj=icrfP[,fromP[j]:toP[j]]
   vpj=vec(Pj)
   paramjdf=paramG[j,]
   paramj=paramjdf[1,3+(1:ncatj)]
   paramjp=paramj
   vpjhat=c( irf( paramjdf, theta, print=0, plot=0 )$ICRF )
   rmsej=sqrt( sum( w*(vpj-vpjhat)^2 )/ncatj )
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
    # d=Ginv(H, MASS=1)%*%g
    d=matSwp(H)%*%g

    # Print(Jac,g,d)

    ok=0
    step=1
    for( lllll in 1:20 ){
     pj=paramj - step*d
     paramjdf[1,3+(1:ncatj)]=pj
     vpjhat=c( irf( paramjdf, theta, print=0, plot=0 )$ICRF )
     rmse1=sqrt( sum( w*(vpj-vpjhat)^2 )/ncatj )
     # Print(llll, lllll, rmsej, rmse1)
     if( rmse1 <= rmsej ){
      pjb=pj[-1]
      if( any(sort(pjb) != pjb) ) break
      if( any(abs(pj) > 10 ) ) break
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
    if( rmsejp > 0 ) rmseimpr=(rmsejp-rmsej)/rmsejp
    else rmseimpr=0
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


 } # end of method == 0
 else{
  ################### in terms of info ##################################



  # for each GPCM item
  for( j in 1: nitems ){

   if( print >= 1 ) cat("\nItem # ", j, "\n", sep="")

   ncatj=ncat[j]
   Pj=infoP[,fromP[j]:toP[j]]
   if( method == 1 ){
    Pj=rowSums(infoP[,fromP[j]:toP[j]])
   }
   vpj=vec(Pj)
   paramjdf=paramG[j,]
   paramj=paramjdf[1,3+(1:ncatj)]
   paramjp=paramj
   if( method == 2 ){
    vpjhat=info_ic( paramj )
    rmsej=sqrt( sum( w*(vpj-vpjhat)^2 )/ncatj )
   }
   else{
    vpjhat=info_i( paramj )
    rmsej=sqrt( sum( w*(vpj-vpjhat)^2 ) )
   }
   rmsejp=rmsej


 if( method == 2 ){
  w=rep(w,ncat[j])
  we=rep(w2,ncat[j])
 }


   # GN iteration
   for( llll in  1:maxiter ){

    # Jacobian:  d vec(Pj) / d paramj    ncatj*npoints x ncatj
    # Jac=JacobianMat( paramj, icrf, ..eps.. = 1e-06 )
    paramjdf[1,3+(1:ncatj)]=paramj
    if( method == 2 ) Jac=JacobianMat( paramj, info_ic )
    else Jac=JacobianMat( paramj, info_i )

    # gradient and method mat and GN direction
    g=-t( t(w*(vpj-vpjhat))%*%Jac )
    maxag=max(abs(g))
    H=t(Jac)%*%(w*Jac)
    # d=Ginv(H, MASS=1)%*%g
    d=matSwp(H)%*%g

   Print(t(d)%*%H%*%d)

    ok=0
    step=1
    for( lllll in 1:20 ){
     pj=paramj - step*d
     paramjdf[1,3+(1:ncatj)]=pj
     if( method == 2 ){
      vpjhat=info_ic( paramj )
      rmse1=sqrt( sum( w*(vpj-vpjhat)^2 )/ncatj )
     }
     else{
      vpjhat=info_i( paramj )
      rmse1=sqrt( sum( w*(vpj-vpjhat)^2 ) )
     }
    Print(llll, lllll, rmsej, rmse1, paramj)
     if( rmse1 <= rmsej ){
      pjb=pj[-1]
      if( any(sort(pjb) != pjb) ) break
      if( any(abs(pj) > 10 ) ) break
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
    if( rmsejp > 0 ) rmseimpr=(rmsejp-rmsej)/rmsejp
    else rmseimpr=0
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


 } # end of method == 1 or 2



 icrfG=irf( paramG, theta, print=0 )$ICRF
 resd=info_func( paramG, theta, print=0, plot=0 )
 infoG=resd$A*resd$dicrf
 iteminfoG=mapply( function(x,y) rowSums(infoG[,x:y])
                  , fromP, toP, SIMPLIFY=FALSE )
 iteminfoG=matrix(unlist(iteminfoG),,nitems)
 colnames(iteminfoG)=iname

 rss=colSums( w*(icrfP-icrfG)^2 )
 rmse_p=sqrt(
  unlist( mapply( function(x,y){sum(rss[x:y])/(y-x+1)}
                     , fromP, toP, SIMPLIFY=FALSE ) ) )
 rmse_pF=formatC(rmse_p, digits=5, format="f")
 rss=colSums( w*(infoP-infoG)^2 )
 rmse_iic=sqrt(
  unlist( mapply( function(x,y){sum(rss[x:y])/(y-x+1)}
                  , fromP, toP, SIMPLIFY=FALSE ) ) )
 rmse_iicF=formatC(rmse_iic, digits=5, format="f")
 rmse_ii=sqrt( colSums( w*(iteminfoP-iteminfoG)^2 ) )
 rmse_iiF=formatC(rmse_ii, digits=5, format="f")

 names(rmse_p)=iname
 names(rmse_iic)=iname
 names(rmse_ii)=iname

 if( print ){
  cat("\nInput and  Converted Item Parameters\n")
  cat("\n Input GPCM Item Parameters\n")
  Print(paramP)
  cat("\n Converted Normal GRM\n")
  Print(paramG)
  cat("Gradients and RMSE minimized:  method =", method, "\n")
  cat(" If method==0, rmse_p is minimized, if method==1,"
     , "rmse_ii is minimized.\n")
  Print( rmse_p, rmse_ii, rmse_iic )
  Print(grad, rmse, fmt="8.5")
  if( print >= 2 ){
   cat("\n ICRFs of the original GPCM\n")
   Print(icrfP, iteminfoP, fmt="5,3")
   cat(" ICRFs of the converted GRM\n")
   Print(icrfG, iteminfoG, fmt="5,3")
  }
  if( print >= 2 ){
   cat("\n Item Category Informations of the original GPCM\n")
   Print(infoP, fmt="5,3")
   cat(" Item Category Informations of the converted GRM\n")
   Print(infoG, fmt="5,3")
  }
 }

 if( plot ){

  nparams=paramP$ncat
  nparams[locB3]=nparams[locB3]+1


  # icrf
  for( j in 1:(nitems) ){
   title=paste("ICRF of ",iname[j]," (ncat=",ncat[j]
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

  # item info
  maxy=max(c(iteminfoP,iteminfoG))
  for( j in 1:(nitems) ){
   ii=rowSums( infoP[,fromP[j]:toP[j]] )
   iih=rowSums( infoG[,fromP[j]:toP[j]] )
   title=paste("Item Information Function of ",iname[j]," (ncat=",ncat[j]
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
   maxy=max(c(infoP),c(infoG))
   for( j in 1:(nitems) ){
    title=paste("Item Category Information of ",iname[j]," (ncat=",ncat[j]
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
     plot(theta,infoP[,fromP[j]+k-1], xlab="", col=1, lty=1, lwd=2
          , xlim=c(min(theta),max(theta)), ylim=c(0,maxy), type="l", ylab="")
     par(new=1)
    }
    for( k in 1:(ncat[j]-1) ){
     plot(theta,infoG[,fromP[j]+k-1], xlab="", col=2, lty=1, lwd=2
          , xlim=c(min(theta),max(theta)), ylim=c(0,maxy), type="b", ylab="")
     par(new=1)
    }
    plot(theta,infoG[,toP[j]], main=title, sub=sub, xlab="", type="b"
         , col=2, lty=1, lwd=2
         , xlim=c(min(theta),max(theta)), ylim=c(0,maxy), ylab="ICRF")
    legend( range(theta)[1]+0.5, 0.9*maxy, c("P","G"), cex=0.8, col=1:2
            ,  pch=c(NA,21), lty=c(1,1), title="model" )
    par(new=0)
   }
  } # end of >= 3

 } # end of plot



 return( named_list( paramG, paramP, rmse_p, rmse_ii, rmse_iic
                    , grad, wtype, wmean, wsd, method, init
                    , icrfP, infoP, iteminfoP
                    , icrfG, infoG, iteminfoG ) )

} # end of fitGn2P_ls



res=fitGn2P_ls( paramA1, theta, plot=1, print=1, wtype=1, method=2, init=0 )




comments(
 '

parambad=cbind( data.frame( name="", type="P", stringsAsFactors=0 )
                ,statbad[,3:7] )
cn=colnames(parambad)
cn=gsub("p1","p", cn)
colnames(parambad)=cn
parambad$name=rownames(statbad)

res=fitGn2P_ls( parambad, theta, plot=1, print=1, wtype=1, method=2 )

# res=fitG2P_ls( parambad[1,], theta, plot=1, print=3, wtype=1, method=1 )





param=paramA1[c(2,5,8),]
param=paramA1[2,]
param=paramA1[c(5,8),]

param=paramA1[c(2,5,8),]

theta=seq(-4,4,length=51)

# res=fitGn2P_ls( param, theta, plot=1, print=2, maxiter=50, wtype=1, method=0 )
res=fitGn2P_ls( param, theta, plot=1, print=1, wtype=1, method=1 )








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
