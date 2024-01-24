#' Approximate Conversion of Partial Credit Items to Graded Response Items
#' Using Logit Transformation
#'

fitL2N <- function( V, print=0, plot=0 ){
 # fitting 2PLM to LRT probability
 # Shin-ichi Mayekawa
 # 20170317
 # plot: 20170623
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
  cat("\nFitting 2PLM to LRT Probability\n")
  Print(nitem, nth, rmse)
  Print(param)
  Print(theta)
 }

 if( plot ){
  theta1=seq(theta[1],theta[nth],length=51)
  res=irf( param, theta1, print=0 )
  irf1=res$IRF
  rmsef=formatC(rmse, wid=7, digits=5)
  title=paste("Fitting 2PLM to LRT Probabilities:  # of items =", nitem
  , "   rmse =", rmsef )
  matplot( theta1, irf1, type="l", main=title )

  for( j in 1:nitem ){
   rmsejf=formatC(rmsej[j], wid=7, digits=5)
   title=paste("item # ",j, "   rmse =", rmsejf)
   temp=cbind(V[j,],irf[,j])
   plot( theta1, irf1[,j], type="l", pch=1, lty=1, col="black"
            , main=title, ylim=c(0,1))
   points( cbind(theta,V[j,]) )
  }

 }


 return( named_list( theta, param, rmse ) )


} # fitL2N


comments('

nclass=5
resm1 <- nIRT( Uc2, plot=1, print=1, nclass=nclass, estrho=1, monotone=1 )

# binary
V1=resm1$V[,seq(2,2*resm1$nitems,2)]
res=fitL2N( t(V1), print=1, plot=1 )
plot(theta0,resm1$theta,type="b")



nclass=7
resm2 <- nIRT( Uc2, plot=1, print=1, nclass=nclass, estrho=1, monotone=1 )

# binary
V1=resm2$V[,seq(2,2*resm2$nitems,2)]
res=fitL2N( t(V1), print=1, plot=1 )
plot(theta0,resm2$theta,type="b")





nclass=5
resm11 <- nIRT( Uc, plot=1, print=1, nclass=nclass, estrho=1, monotone=1 )

# binary
V1=resm11$V[,seq(2,2*resm11$nitems,2)]
res=fitL2N( t(V1), print=1, plot=1 )
plot(theta0,resm11$theta,type="b")


nclass=7
resm21 <- nIRT( Uc, plot=1, print=1, nclass=nclass, estrho=1, monotone=1 )

# binary
V1=resm21$V[,seq(2,2*resm21$nitems,2)]
res=fitL2N( t(V1), print=1, plot=1 )
plot(theta0,resm21$theta,type="b")


')













