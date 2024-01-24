#' Conversion of LRT to IRT by Weighted LS
#'
#' ######################## NO GOOD!! ##################################
#'
#' @param V item x class probability matrix
#' @param print = 1 to print the estimated IRT item parameters \cr
#' = 2 to print the irf.
#' @param plot = 1 to plot the main result \cr
#' = 2 to plot irf of each item.
#' @param title title string
#' @param wtype = 1 to use dnomr(theta) as the weight \cr
#' = 0 to use no weight.
#' @param SQUAREM = 3 :   See the help of iSQUAREM in lazy.accel package.
#' @param nSQUAREM when to star iSQUAREM
#' @param minalpha = -999 :   See the help of iSQUAREM in lazy.accel package.
#' @param  maxalpha = -1 :   See the help of iSQUAREM in lazy.accel package.
#' @param  always = 0 :   See the help of iSQUAREM in lazy.accel package.
#' @param  reset1 = 1 :   See the help of iSQUAREM in lazy.accel package.
#' @param  reset2 = 2 :   See the help of iSQUAREM in lazy.accel package.
#'
#'
#' @details
#' The LS criterion: \cr
#' \code{ ssq( t(V) - irf(theta|item parameters) )} \cr
#' will be minimized  with respect to theta and item parameters.
#'
#'
#' @return A list of: \cr
#' theta The estimated theta value for each latent rank \cr
#' param IRT item parameter data.frame \cr
#' rmse The rmse stat. \cr
#' wtype The type of weight.
#'
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
#' param <- paramB1[c(1:3,7:9,13:15),]
#' thmin <- -2; thmax <- 2; npoint <- 5
#' N <- 1000
#'
#' # discrete theta
#' # theta0 <- seq(thmin,thmax,length=npoint)
#' theta0 <- c(-2, -1, 0, 2, 3)
#' theta <- unlist(lapply( theta0, rep, round(N/npoint) ))
#' res2 <- gendataIRT( 1, paramB1, theta=theta, compress=1 )
#' Uc <- as.data.frame(res2$U)
#' ncat <- res2$ncat
#' type <- res2$type
#'
#' # lrt parameters
#' nclass <- 5
#' resm1 <- uLRT( Uc, nclass=nclass, estrho=1, monotone=1, alpha=20
#'                                   , maxiter=20, plot=1, print=1 )
#' V1 <- resm1$V[,seq(2,2*resm1$nitems,2)]
#'
#' # conversion
#' res <- fitI2L_ls( t(V1), print=1, plot=1, wtype=1, maxiter=50 )
#' plot(theta0, res$theta,type="b", main="original theta vs recovered theta")
#'
#' @export
#'

fitI2L_ls <- function( V, print=0, plot=0, title=NULL, wtype=0
                       , maxiter=1000, eps=1e-6
                       , SQUAREM=3, nSQUAREM=1, minalpha=-999, maxalpha=-1
                       , always=1, reset1=0, reset2=1 ){
 # fitting 2PLM IRT to LRT probability
 # Shin-ichi Mayekawa
 # 20170317
 # plot: 20170623
 # fitL2N renamed as fit I2L: 20170803
 # title: 20171113
 # negative a: 20180111
 # LS version: 20180112
 # nlm -> optim: 20180112
 # combine param and theta: 20180112
 #

 # V is nitem x npoint
 #

 rssab <- function( paramj, tV, wtype=0 ){
  rss=( tV-icrfB( matrix(paramj,,3), theta ) )
  if( wtype == 1 ) rss=sum( dnorm(theta)*rss*rss )
  else  rss=sum( rss*rss )
  return( rss )
 } # end of rssab

 rssth <- function( theta, tV, wtype=0 ){
  rss=( tV-icrfB( param, theta ) )
  if( wtype == 1 ) rss=sum( dnorm(theta)*rss*rss )
  else  rss=sum( rss*rss )
  return( rss )
 } # end of rssth

 normth <- function( pp ){
  # normalize theta in iSQUAREM
  param[,1:2]=matrix(pp[1:(2*nitem)],,2)
  theta=pp[-(1:(2*nitem))]

  ms=mands(theta)
  mth=ms[1,1]; sth=ms[2,1]
  theta=(theta-mth)/sth
  param[,2]=(param[,2]-mth)/sth
  param[,1]=param[,1]*sth

  pp=c( param[,1:2], theta )
  return( pp )
 } # end of normth


 rmse <- function( pp ){
  # get rmse in iSQUAREM
  param[,1:2]=matrix(pp[1:(2*nitem)],,2)
  theta=pp[-(1:(2*nitem))]
  Print(theta,param)
  rss=rssab( param, tV, wtype=wtype )
  rmse=sqrt(rss/nitem/npoint)
  return( rmse )
 }

 nitem=nrow(V)
 npoint=ncol(V)


 itemname=rownames(V)
 if( is.null(itemname) ) itemname=1:nitem

 # initial
 resf=fitI2L( V, print=0, plot=0 )
 param=as.matrix( resf$param[,c("p1","p2","p3")] )
 theta=c(resf$theta)



 tV=t(V)

 rmsep=sqrt(rssab( param, tV, wtype=wtype )/nitem/npoint)

 if( print ){
  cat("\nFitting 2PLM to LRT Probability by LS\n")
  if( !is.null(title) ) Print(title)
  Print(nitem, npoint, rmsep, wtype)
  if( print >= 2 ){
   Print("initial")
   Print(param,theta)
  }
 }

 control=list(reltol=eps, maxit=maxiter)



 ###################################################################(1)
 # store param
 # initialize inline SQUAREM
 rmv(iSQUAREM)
 pp=c( param[,1:2], theta )
 iSQUAREM=generate_iSQUAREM( pp )
 param0=param
 theta0=theta
 ###################################################################(1)



 for( llll in 1:maxiter ){

  pp=c(param,theta)
  res=optim( pp, rmse, control=control, method="BFGS" )
  param[,1:2]=matrix(pp[1:(2*nitem)],,2)
  theta=pp[-(1:(2*nitem))]
Print("***", param,theta)
  # normalize
  ms=mands(theta)
  mth=ms[1,1]; sth=ms[2,1]
  theta=(theta-mth)/sth
  param[,2]=(param[,2]-mth)/sth
  param[,1]=param[,1]*sth


  rss=rssab( param, tV, wtype=wtype )
  # Print( llll, theta, rss )

  rmse=sqrt(rss/nitem/npoint)

  rmseimpr=(rmsep-rmse)/rmsep

  if( print >= 2 ){
   Print( llll, rmse, rmseimpr )
  }

  if( rmseimpr <= eps ) break



  ###################################################################(2)
  if( SQUAREM  &  llll >= nSQUAREM ){
   # inline SQUAREM
   pp=c( param[,1:2], theta )
   # Print(locSQEM,locms, nparam1, nparam, param)
   res=iSQUAREM( pp
                 , enforce_constraints=normth, badness_of_fit=rmse
                 , SQUAREM=SQUAREM, minalpha=minalpha, maxalpha=maxalpha
                 , bof_value=NULL, always=always, reset1=reset1, reset2=reset2
                 , print=0, debug=0 )
   #Print(pp,res$param)
   # reshape
   param[,1:2]=matrix(res$param[1:(2*nitem)],,2)
   theta=res$param[-(1:(2*nitem))]

  }
  ###################################################################(2)


  # next iteration
  rmsep=rmse


 } # end of llll loop


 if( print >= 1 ){
  cat("\niteration terminated with ", llll, " iterations.\n")
  Print( llll, rmse, rmseimpr )
 }


 resf$param[,c("p1","p2","p3")]=param
 param=resf$param


 # calculate irf
 res=irf( param, theta, print=0 )
 irf=res$IRF
 if( wtype == 1 ) rss=colSums( dnorm(theta)*(t(V)-irf)^2 )
 else rss=colSums( (t(V)-irf)^2 )
 rmse=sqrt(sum(rss/npoint)/nitem)
 rmsej=sqrt(rss/npoint)


 if( print ){
  cat("\nFitting 2PLM to LRT Probability by LS\n")
  if( !is.null(title) ) Print(title)
  Print(nitem, npoint, rmse, wtype)
  Print(param)
  Print(theta)
  if( print >= 2 ){
   Print(V)
   Print(irf)
  }
 }

 if( plot ){
  theta1=seq(theta[1],theta[npoint],length=51)
  res=irf( param, theta1, print=0 )
  irf1=res$IRF
  rmsef=formatC(rmse, wid=7, digits=5)
  title1=paste("Fitting 2PLM IRT to LRT Probabilities (LS):  # of items ="
               , nitem, " rmse =", rmsef )
  if( !is.null(title) ) title1=paste(title,"\n",title1)
  matplot( theta1, irf1, type="l", main=title1 )
  for( j in 1:nitem ){
   points( cbind(theta,V[j,]) )
  }

  if( plot >= 2 ){
   for( j in 1:nitem ){
    rmsejf=formatC(rmsej[j], wid=7, digits=5)
    title2=paste("item # ", j, "(",  itemname[j], ")   rmse =", rmsejf)
    temp=cbind(V[j,],irf[,j])
    title2=paste(title1,"\n", title2, sep="")
    plot( theta1, irf1[,j], type="l", pch=1, lty=1, col="black"
          , main=title2, ylim=c(0,1))
    points( cbind(theta,V[j,]) )
   }
  }
 }


 return( named_list( theta, param, rmse, wtype, llll, eps ) )


} # end of fitI2L_ls




res=fitI2L_ls( V, print=1, plot=1, maxiter=1000, wtype=1, SQUAREM=0 )






comments(
 '

 V0=read.csv("d:/RPGM/LRT/data2/V_24.csv", header=1, stringsAsFactors=0)
 prior_prob=c( 0.09, 0.15, 0.2, 0.225, 0.2, 0.15, 0.01 )
 prior_prob=prior_prob/sum(prior_prob)
 V=as.matrix(V0[,-c(1:3,9,10)])
 rownames(V)=V0[,2]
 nclass=ncol(V)
 classname=paste("class",1:nclass,sep="")

 # convert to IRT
 resL2I=fitI2L( V, print=1, plot=1 )
 paramL2I=resL2I$param










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
 res=fitI2L( t(V1), print=1, plot=1, title="testting title statement" )
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




 '
)













