
plot_result <- function( result ){
 # plot the result of fitG2P etc
 # Shin-ichi Mayekawa
 # 20221208,09
 #
 # with legend
 #

 #
 # layout should be used outside of this function.
 #

 paramNew=result$paramNew
 itypeNew=paramNew$type

 ncat=paramNew$ncat
 iname=paramNew$name
 method=result$method
 wtype=result$wtype
 theta=result$theta
 fromP=result$fromP
 toP=result$toP


 paramOld=result$paramOld
 itypeOld=paramOld$type

 nparam=ncat
 if( itypeOld %in% c("B3","Bn3") ) nparam=nparam+1

 rmse_p=result$rmse_p
 rmse_ii=result$rmse_ii
 rmse_iic=result$rmse_iic

 rmse_pF=formatC(rmse_p, digits=5, format="f")
 rmse_iiF=formatC(rmse_ii, digits=5, format="f")
 rmse_iicF=formatC(rmse_iic, digits=5, format="f")

 icrfNew=result$icrfNew
 icifNew=result$icifNew
 iifNew=result$iifNew[,,drop=0]
 icrfOld=result$icrfOld
 icifOld=result$icifOld
 iifOld=result$iifOld[,,drop=0]

 plot=4
 nitems=1




 Print(icrfNew,icrfOld, iifNew, iifOld, icifNew, icifOld, theta, fromP,toP)

  # icrf
  for( j in 1:(nitems) ){
   title=paste("ICRF of ",iname," (", itypeOld, " -> ", itypeNew
               ,")  ncat=", ncat
               , ")  rmse_p=", rmse_pF, " (method=", method
               , ", wtype=", wtype, ")", sep="")
   sub=paste( "paramOld(", itypeOld,") = "
              , paste(format(paramOld[j,3+(1:nparam)],digits=3)
                      , collapse=",  ")
              , "\nparamNew(", itypeNew,") = "
              , paste(format(paramNew[j,3+(1:nparam)],digits=3)
                      , collapse=",  "), sep=""    )
   for( k in 1:(ncat) ){
    plot(theta,icrfOld[,fromP[j]+k-1], xlab="", col=1, lty=1, lwd=2
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="l", ylab="")
    par(new=1)
   }
   for( k in 1:(ncat-1) ){
    plot(theta,icrfNew[,fromP+k-1], xlab="", col=2, lty=1, lwd=2
         , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="b", ylab="")
    par(new=1)
   }
   plot(theta,icrfNew[,toP], main=title, sub=sub
        , xlab="", col=2, lty=1, lwd=2
        , xlim=c(min(theta),max(theta)), ylim=c(0,1), type="b", ylab="ICRF")
   legend( range(theta)[1]+0.5, 0.9, c("New","Old"), cex=0.8, col=1:2
           ,  pch=c(NA,21), lty=c(1,2), title="model" )
   par(new=0)
  }

  # item info
  maxy=max(c(iifOld)*1.1)
  for( j in 1:(nitems) ){
   title=paste("Item Information Function of "
               ,iname," (", itypeOld, " -> ", itypeNew
               ,")  ncat=", ncat
               , ")  rmse_p=", rmse_pF, " (method=", method
               , ", wtype=", wtype, ")", sep="")
   sub=paste( "paramOld(", itypeOld,") = "
              , paste(format(paramOld[j,3+(1:nparam)],digits=3)
                      , collapse=",  ")
              , "\nparamNew(", itypeNew,") = "
              , paste(format(paramNew[j,3+(1:nparam)],digits=3)
                      , collapse=",  "), sep=""    )
   matplot( theta, cbind(iifOld,iifNew)
            , main=title, sub=sub,  pch=c(NA,21), xlab=""
            , type="l", ylab="information", ylim=c(0,maxy), lty=c(1,2), lwd=2 )
   legend( range(theta)[1]+0.5, 0.9*maxy, c("Old","New"), cex=1, col=1:2
           ,  pch=c(NA,NA), lty=c(1,1), title="model" )
  }

  if( plot >= 3 ){
   # cat info
   maxy=max(c(icifOld)*1.1)
  Print(maxy)
   for( j in 1:(nitems) ){
    title=paste("Item Category Information of "
                ,iname," (", itypeOld, " -> ", itypeNew
                ,")  ncat=", ncat
                , ")  rmse_p=", rmse_pF, " (method=", method
                , ", wtype=", wtype, ")", sep="")
    sub=paste( "paramOld(", itypeOld,") = "
               , paste(format(paramOld[j,3+(1:nparam)],digits=3)
                       , collapse=",  ")
               , "\nparamNew(", itypeNew,") = "
               , paste(format(paramNew[j,3+(1:nparam)],digits=3)
                       , collapse=",  "), sep=""    )
    for( k in 1:(ncat) ){
     plot(theta,icifOld[,fromP+k-1], xlab="", col=1, lty=1, lwd=2, pch=k
          , xlim=c(min(theta),max(theta)), ylim=c(0,maxy), type="b", ylab="")
     par(new=1)
    }
    for( k in 1:(ncat[j]-1) ){
     plot(theta,icifNew[,fromP+k-1], xlab="", col=2, lty=1, lwd=2, pch=k
          , xlim=c(min(theta),max(theta)), ylim=c(0,maxy), type="b", ylab="")
     par(new=1)
    }
    plot(theta,icifNew[,toP], main=title, sub=sub, xlab="", type="b"
         , col=2, lty=1, lwd=2, pch=ncat[j]
         , xlim=c(min(theta),max(theta)), ylim=c(0,maxy), ylab="ICRF")
    legend( range(theta)[1]+0.5, 0.9*maxy, c("Old","New"), cex=0.8, col=1:2
            ,  pch=c(NA,21), lty=c(1,1), title="model" )
    par(new=0)
   } # end of j
  } # end of >= 3


} # end of plot_result




par(oma=c(0,0,0,0))
mat <- matrix(1:9,3,3)
layout(mat)

plot_result( res0 )
plot_result( res1 )
plot_result( res2 )

par(oma=c(0,0,0,0))
layout(1)

