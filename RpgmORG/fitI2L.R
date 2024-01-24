#' Approximate Conversion of Partial Credit Items to Graded Response Items
#' Using Logit Transformation
#'
#' @param V item x class probability matrix
#' @param print = 1 to print the result
#' @param plot = 1 to plot the result
#'
#' @examples
#' #
#' #### In the following examples, maxiter is set to 20 which is
#' #### not large enough to obtain convergence.
#' ####
#' #
#' #
#' #
#' set.seed(1701)
#'
#' param=paramB1[c(1:3,7:9,13:15),]
#' thmin=-2; thmax=2; npoint=5
#' N=1000
#'
#' # discrete theta
#' # theta0=seq(thmin,thmax,length=npoint)
#' theta0=c(-2, -1, 0, 2, 3)
#' theta=unlist(lapply( theta0, rep, round(N/npoint) ))
#' res2 <- gendataIRT( 1, paramB1, theta=theta, compress=1 )
#' Uc=as.data.frame(res2$U)
#' ncat=res2$ncat
#' type=res2$type
#'
#' # lrt parameters
#' nclass=5
#' resm1 <- uLRT( Uc, nclass=nclass, estrho=1, monotone=1, alpha=20
#'                                   , maxiter=20, plot=1, print=1 )
#' V1=resm1$V[,seq(2,2*resm1$nitems,2)]
#'
#' # conversion
#' res=fitI2L( t(V1), print=1, plot=1 )
#' plot(theta0,resm1$theta,type="b")
#'
#' @export
#'

fitI2L <- function( V, print=0, plot=0 ){
 # fitting 2PLM IRT to LRT probability
 # Shin-ichi Mayekawa
 # 20170317
 # plot: 20170623
 # fitL2N renamed as fit I2L: 20170803
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
  theta1=seq(theta[1],theta[nth],length=51)
  res=irf( param, theta1, print=0 )
  irf1=res$IRF
  rmsef=formatC(rmse, wid=7, digits=5)
  title=paste("Fitting 2PLM IRT to LRT Probabilities:  # of items =", nitem
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


} # end of fitI2L




comments('

nclass=5

resm1 <- uLRT( Uc2, maxiter=1000, nclass=nclass, estrho=1, monotone=1
             , print=1, plot=1, alpha=10*c(1,2,3,2,1))

resm1 <- uLRT( Uc2, maxiter=1000, nclass=nclass, estrho=1, monotone=1
             , print=1, plot=1, alpha=-10)

resm1 <- uLRT( Uc2, maxiter=1000, nclass=nclass, estrho=1, monotone=1
             , print=1, plot=1, alpha=5)

resm1 <- uLRT( Uc2, maxiter=1000, nclass=nclass, estrho=1, monotone=1
             , print=1, plot=1, alpha=30)

# binary
V1=resm1$V[,seq(2,2*resm1$nitems,2)]
res=fitI2L( t(V1), print=1, plot=1 )
plot(theta0,resm1$theta,type="b")








nclass=5
res1 <- uLRT( Uc2, plot=1, print=1, nclass=nclass, estrho=1, monotone=1 )

# binary
V1=res1$V[,seq(2,2*res1$nitems,2)]
res=fitI2L( t(V1), print=1, plot=1 )
plot(theta0,res$theta,type="b")





nclass=7
resm1 <- uLRT( Uc2, plot=1, print=1, nclass=nclass, estrho=1, monotone=1
           , alpha=10*c(10,20,40,60,40,20,10))

nclass=7
resm1 <- uLRT( Uc2, plot=1, print=1, nclass=nclass, estrho=1, monotone=1
         , alpha=-50)




')













