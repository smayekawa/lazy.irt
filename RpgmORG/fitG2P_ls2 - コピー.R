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
#' paramG1 <- fitG2P_ls( paramP1, plot=1, print=1 )
#'
#'
#'
#' @export
#'

fitG2P_ls <- function( paramP, theta=NULL
                    , info=0, wtype=0, wmean=0, wsd=1, DinP=1
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
 # roxgen2: 20180211
 # info: 20180212
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
  # item category information function
  paramjdf[1,4:(4+length(paramnum)-1)]=paramnum
  dicrf=dirf( paramjdf, theta, print=0, plot=0 )$dICRF
  icrf=irf( paramjdf, theta, print=0, plot=0 )$ICRF
  # matplot( theta, dicrf^2/icrf, type="l")
  info=vec( dicrf^2/icrf )
  return( info )
 } # end of info_ic

 get_rmse <- function( paramj ){
  # everything else is from outside
  paramjdf[1,3+(1:ncatj)]=paramj
  vpjhat=c( irf( paramjdf, theta, print=0, plot=0 )$ICRF )
  rmse=sqrt( sum( w*(vpj-vpjhat)^2 /npoints/ncatj ) )
  return( rmse )
 } # end of get_rmse

 get_rmse_info <- function( paramj ){
  # everything else is from outside
  infohat=info_ic( paramj )
  rmse=sqrt( sum( w*(vpj-infohat)^2 /npoints/ncatj ) )
  return( rmse )
 } # end of get_rmse_info



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
  cat("\nFitting Graded Response Model to Generalized Partial Credit Model\n")
  cat("  LS criterion in terms of ")
  if( info ) cat( "item category informations\n" )
  else cat( "item category response functions\n" )
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

 # item category info: Info is npoints x ncat
 dicrf=dirf( param, theta, zero=1, print=0, plot=0 )$dICRF
 icrf=irf( param, theta, zero=1, print=0, plot=0 )$ICRF
 Info=dicrf^2/icrf
 iteminfo=mapply( function(x,y) rowSums(Info[,x:y])
                      , fromP, toP, SIMPLIFY=FALSE )
 iteminfo=matrix(unlist(iteminfo),,nitems)
 colnames(iteminfo)=iname

 # initial
 paramG=fitG2P( paramP, npoints=21, plot=0, print=0, debug=0, wtype=wtype )

 rmse=numeric(nitems)
 names(rmse)=iname
 grad=matrix(NA,max(ncat),nitems)
 rownames(grad)=paste("p",1:nrow(grad),sep="")
 colnames(grad)=iname
 grad=t(grad)

 if( info == 0 ){
  ################### in terms of icrf ##################################


  # for each GPCM item
  for( j in 1: nitems ){

   if( print >= 1 ) cat("\nItem # ", j, "\n", sep="")

   ncatj=ncat[j]
   Pj=P[,fromP[j]:toP[j]]
   vpj=vec(Pj)
   paramjdf=paramG[j,]
   paramj=paramjdf[1,3+(1:ncatj)]
   paramjp=paramj
   vpjhat=c( irf( paramjdf, theta, print=0, plot=0 )$ICRF )
   rmsej=sqrt( sum( w*(vpj-vpjhat)^2 )/npoints/ncatj )
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
    d=solve(H)%*%g

    # Print(Jac,g,d)

    ok=0
    step=1
    for( lllll in 1:20 ){
     pj=paramj - step*d
     paramjdf[1,3+(1:ncatj)]=pj
     vpjhat=c( irf( paramjdf, theta, print=0, plot=0 )$ICRF )
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


 } # end of info == 0
 else{
  ################### in terms of info ##################################



  # for each GPCM item
  for( j in 1: nitems ){

   if( print >= 1 ) cat("\nItem # ", j, "\n", sep="")

   ncatj=ncat[j]
   Pj=Info[,fromP[j]:toP[j]]
   vpj=vec(Pj)
   paramjdf=paramG[j,]
   paramj=paramjdf[1,3+(1:ncatj)]
   paramjp=paramj
   vpjhat=info_ic( paramj )
   rmsej=sqrt( sum( w*(vpj-vpjhat)^2 )/npoints/ncatj )
   rmsejp=rmsej

   # GN iteration
   for( llll in  1:maxiter ){

    # Jacobian:  d vec(Pj) / d paramj    ncatj*npoints x ncatj
    # Jac=JacobianMat( paramj, icrf, ..eps.. = 1e-06 )
    paramjdf[1,3+(1:ncatj)]=paramj
    Jac=w2*JacobianMat( paramj, info_ic )
#Print( vpj, vpjhat, paramj, Jac)
    # Print(Jac-Jac0)
    # gradient and info mat and GN direction
    g=-t( t(w*(vpj-vpjhat))%*%Jac )
    maxag=max(abs(g))
    H=t(Jac)%*%Jac
    # d=Ginv(H, MASS=1)%*%g
    d=solve(H)%*%g

    # Print(Jac,g,d)

    ok=0
    step=1
    for( lllll in 1:20 ){
     pj=paramj - step*d
     paramjdf[1,3+(1:ncatj)]=pj
     vpjhat=info_ic( paramj )
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


 } # end of info == 1



 PG=irf( paramG, theta, print=0 )$ICRF
 resd=info_func( paramG, theta, print=0, plot=0 )
 Infohat=resd$A*resd$dicrf
 iteminfohat=mapply( function(x,y) rowSums(Infohat[,x:y])
                  , fromP, toP, SIMPLIFY=FALSE )
 iteminfohat=matrix(unlist(iteminfohat),,nitems)
 colnames(iteminfohat)=iname

 rss=colSums( w*(P-PG)^2 )
 rmse_p=sqrt(
  unlist( mapply( function(x,y){sum(rss[x:y])/(y-x+1)/npoints}
                     , fromP, toP, SIMPLIFY=FALSE ) ) )
 rmse_pc=formatC(rmse_p, digits=5, format="f")

 rmse_i=sqrt( colSums( w*(Info-Infohat)^2 )/npoints/ncatj )
 rss=colSums( w*(Info-Infohat)^2 )
 rmse_i=sqrt(
  unlist( mapply( function(x,y){sum(rss[x:y])/(y-x+1)/npoints}
                  , fromP, toP, SIMPLIFY=FALSE ) ) )
 rmse_ic=formatC(rmse_i, digits=5, format="f")

 rmse_ii=sqrt( colSums( w*(iteminfo-iteminfohat) )/npoints )

 names(rmse_p)=iname
 names(rmse_i)=iname
 names(rmse_ii)=iname

 if( print ){
  cat("\nInput and  Converted Item Parameters\n")
  cat("\nInput and  GPCM Item Parameters\n")
  Print(paramP)
  cat("\n  Converted GRM D=1.7\n")
  Print(paramG)
  cat("Gradients and RMSE minimized:  info =", info, "\n")
  Print( rmse_p, rmse_i, rmse_ii )
  Print(grad, rmse, fmt="8.5")
  if( print >= 2 ){
   cat("\n ICRFs of the original GPCM\n")
   Print(P, iteminfo, fmt="5,3")
   cat(" ICRFs of the converted GRM\n")
   Print(PG, iteminfohat, fmt="5,3")
  }
  if( print >= 2 ){
   cat("\n Item Category Informations of the original GPCM\n")
   Print(Info, fmt="5,3")
   cat(" Item Category Informations of the converted GRM\n")
   Print(Infohat, fmt="5,3")
  }
 }

 if( plot ){

  # icrf
  for( j in 1:(nitems) ){
   title=paste("ICRF of ",iname[j]," (ncat=",ncat[j]
               , ")  rmse=", rmse_pc[j], " (info=", info
               , ", wtype=", wtype, ")", sep="")
   sub=paste( "paramP ="
              , paste(format(paramP[j,3+(1:ncat[j])],digits=3)
                      , collapse=",  ")
              , "\nparamG ="
              , paste(format(paramG[j,3+(1:ncat[j])],digits=3)
                      , collapse=",  ")
   )
   for( k in 1:(ncat[j]) ){
    plot(theta,P[,fromP[j]+k-1], xlab="", col=1, lty=1, lwd=2
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="")
    par(new=1)
   }
   for( k in 1:(ncat[j]-1) ){
    plot(theta,PG[,fromP[j]+k-1], xlab="", col=2, lty=1, lwd=2
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="b", ylab="")
    par(new=1)
   }
   plot(theta,PG[,toP[j]], main=title, sub=sub, xlab="", col=2, lty=1, lwd=2
        , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="b", ylab="ICRF")
   legend( range(theta)[1]+0.5, 0.9, c("P","G"), cex=0.8, col=1:2
           ,  pch=c(NA,21), lty=c(1,2), title="model" )
   par(new=0)
  }

  # item info
  maxy=max(c(iteminfo,iteminfohat))
  for( j in 1:(nitems) ){
   ii=rowSums( Info[,fromP[j]:toP[j]] )
   iih=rowSums( Infohat[,fromP[j]:toP[j]] )
   title=paste("Item Information Function of ",iname[j]," (ncat=",ncat[j]
               , ")  rmse=", rmse_ic[j], " (info=", info
               , ", wtype=", wtype, ")", sep="")
   sub=paste( "paramP ="
              , paste(format(paramP[j,3+(1:ncat[j])],digits=3)
                      , collapse=",  ")
              , "\nparamG ="
              , paste(format(paramG[j,3+(1:ncat[j])],digits=3)
                      , collapse=",  ")
   )
   matplot( theta, cbind(ii,iih), main=title, sub=sub,  pch=c(NA,21), xlab=""
          , type="l", ylab="information", ylim=c(0,maxy), lty=c(1,2), lwd=2 )
   legend( range(theta)[1]+0.5, 0.9*maxy, c("P","G"), cex=1, col=1:2
           ,  pch=c(NA,NA), lty=c(1,1), title="model" )
  }


  # cat info
  maxy=max(c(Info),c(Infohat))
  for( j in 1:(nitems) ){
   title=paste("Item Category Information of ",iname[j]," (ncat=",ncat[j]
               , ")  rmse=", rmse_ic[j], " (info=", info
               , ", wtype=", wtype, ")", sep="")
   sub=paste( "paramP ="
              , paste(format(paramP[j,3+(1:ncat[j])],digits=3)
                      , collapse=",  ")
              , "\nparamG ="
              , paste(format(paramG[j,3+(1:ncat[j])],digits=3)
                      , collapse=",  ")
   )
   for( k in 1:(ncat[j]) ){
    plot(theta,Info[,fromP[j]+k-1], xlab="", col=1, lty=1, lwd=2
         , xlim=c(min(theta),max(theta)), ylim=c(0,maxy), type="l", ylab="")
    par(new=1)
   }
   for( k in 1:(ncat[j]-1) ){
    plot(theta,Infohat[,fromP[j]+k-1], xlab="", col=2, lty=1, lwd=2
         , xlim=c(min(theta),max(theta)), ylim=c(0,maxy), type="b", ylab="")
    par(new=1)
   }
   plot(theta,Infohat[,toP[j]], main=title, sub=sub, xlab="", type="b"
        , col=2, lty=1, lwd=2
        , xlim=c(min(theta),max(theta)), ylim=c(0,maxy), ylab="ICRF")
   legend( range(theta)[1]+0.5, 0.9*maxy, c("P","G"), cex=0.8, col=1:2
           ,  pch=c(NA,21), lty=c(1,1), title="model" )
   par(new=0)
  }

 } # end of plot

 # convert the result to data frame

 return( named_list(paramG, paramP, rmse, grad, wtype, wmean, wsd, info) )

} # end of fitG2P_ls



param=paramA1[c(5,8),]

theta=seq(-4,4,length=51)

# res=fitG2P_ls( param, theta, plot=1, print=2, maxiter=50, wtype=1, info=0 )
res=fitG2P_ls( param, theta, plot=1, print=2, maxiter=50, wtype=1, info=1 )


comments(
 '

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
