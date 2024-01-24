plot_uIRT <- function( irfetc, what="icrf" ){
 # plot irf etc
 # Shin-ichi Mayekawa
 #

 # trf

 matplot( theta, cbind(TRF, TRF_LO), type="l", col=1, lwd=2, lty=c(1,3)
          , main="TRF and TRF with the locally best item weights" )
 legend( "topleft", c("TRF","TRF with LO")
         , cex=0.8, col=1, lwd=3, lty=c(1,3)
         , pch=NA, title="legend" )


 # plot test info
 matplot( theta,cbind(info, info_LO, info_LOW), type="l"
          , xlim=c(min(theta),max(theta)), ylim=c(0,max(info_LO))
          , main="Test Information Functions", ylab="information")
 legend( range(theta)[1], max(info_LO)-.1, c("info","info_LO","info_LOW")
         , cex=0.8, col=c("black", "red", "green", "blue")
         , pch=NA, lty=c(1:3), title="legend" )

 # plot item info
 matplot( theta,cbind(info_item_LO), type="l"
          , col=c("black", "red", "green", "blue"), lty=c(1:3)
          , xlim=c(min(theta),max(theta)), ylim=c(0,max(info_item_LO))
          , main="Item Information Functions", ylab="information")
 if( legend > 0 ){
  legend( range(theta)[1], max(info_item_LO)-.1, iname
          , cex=0.8, col=c("black", "red", "green", "blue")
          , pch=NA, lty=c(1:3), title="legend" )
 }


 if( plot >= 2 ){

  for( j in 1:(nitems) ){
   title=paste("Item Information Function of "
               , iname[j]," ( type = ",type[j],", ncat=",ncat[j], " )")
   pp=param[j,4:(ncat[j]+3)]
   pp[is.na(pp)]=0
   sub=paste("param = ", paste(format(pp,digits=3)
                               , collapse=",  "), "  (with weights)")
   if( length(grep("^B[[:digit:]]*$", param$type[j])) > 0 ){
    pp=param[j,4:(ncat[j]+4)]
    pp[is.na(pp)]=0
    sub=paste("param = ", paste(format(pp,digits=3), collapse=",  "))
   }
   maxi=max(info_item_LO)
   if( maxinfo > 0 ) maxi=maxinfo
   plot(theta,info_item_LO[,j], main=title, sub=sub
        , xlim=c(min(theta),max(theta)), ylim=c(0,maxi)
        , type="l", ylab="information")
  }
 }

 matplot( theta,w_LO, type="l"
          , col=c("black", "red", "green", "blue"), lty=c(1:3)
          , xlim=c(min(theta),max(theta)), ylim=c(0,max(w_LO))
          , main="Locally Best Item Weights", ylab="weight")
 if( legend > 0 ){
  legend( range(theta)[2]-1, range(w_LO)[2]/4*3, iname
          , cex=0.8, col=c("black", "red", "green", "blue")
          , pch=NA, lty=c(1:3), title="legend" )
 }

 if( plot >= 3 ){

  for( j in 1:(nitems) ){
   title=paste("Locally Best Item Category Weights for"
               , iname[j]," ( type = ",type[j],", ncat=",ncat[j], " )")
   pp=param[j,4:(ncat[j]+3)]
   pp[is.na(pp)]=0
   sub=paste( "param = ", paste(format(pp,digits=3)
                                , collapse=",  ") )
   if( length(grep("^B[[:digit:]]*$", param$type[j])) > 0 ){
    pp=param[j,4:(ncat[j]+4)]
    pp[is.na(pp)]=0
    sub=paste("param = ", paste(format(pp,digits=3), collapse=",  "))
   }
   maxi=max(v_LO)
   #Print(j,theta,v_LO[,fromP[j]:toP[j]])
   matplot(theta,v_LO[,fromP[j]:toP[j]], main=title, sub=sub
           , xlim=c(min(theta),max(theta)), ylim=c(0,maxi)
           , type="l", ylab="weight")
  }

 }


} # end of plot_uIRT

