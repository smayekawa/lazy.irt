#' Approximate Conversion of Partial Credit Items to Graded Response Items
#' Using Logit Transformation
#'

fitL2N <- function( V, print=0, plot=0 ){
 # fitting 2PLM to NTT probability
 # Shin-ichi Mayekawa
 # 20170317
 #

 # V is nitem x nth

 nitem=nrow(V)
 nth=ncol(V)

 # logit transformation:  ltV is nth x nitem.
 ltV=logit(t(V))

 # column center the logit
 mean=colMeans(ltV)
 ltV=ltV-matrix(1,nth,1)%*%mean

 # pca of itV
 res=princ( ltV, ndim=1, print=0 )
 theta=res$F
 a=res$A

 # get b from mean and a
 b=-mean/a


 # parameter data frame
 param=data.frame( name=rownames(V), type="B", ncat=2, p1=a, p2=b, p3=0 )

 # calculate irf
 res=irf( param, theta, print=0 )
 irf=res$IRF
 Print(V,irf)
 rss=colSums( (t(V)-irf)^2 )
 rmse=sqrt(sum(rss/nth)/nitem)
 rmsej=sqrt(rss/nth)

 if( print ){
  cat("\nFitting 2PLM to NTT Probability\n")
  Print(nitem, nth, rmse)
  Print(param)
  Print(theta)
 }

 if( plot ){
  rmsef=formatC(rmse, wid=7, digits=5)
  title=paste("Fitting 2PLM to NTT Probabilities:  # of items =", nitem
  , "   rmse =", rmsef )
  matplot( theta, irf, type="b", main=title )

  for( j in 1:nitem ){
   rmsejf=formatC(rmsej[j], wid=7, digits=5)
   title=paste("item # ",j, "   rmse =", rmsejf)
   temp=cbind(V[j,],irf[,j])
   matplot( theta, temp, type=c("p","l"), pch=1, lty=1, col="black"
            , main=title)
  }

 }



 return( named_list( theta, param, rmse ) )


} # fitL2N



res=fitL2N( V, print=1, plot=1 )















