#' Conversion of LRT to IRT by Weighted LS
#'
#' @param V item x class probability matrix
#' @param print = 1 to print the estimated IRT item parameters \cr
#' = 2 to print the irf.
#' @param plot = 1 to plot the main result \cr
#' = 2 to plot irf of each item.
#' @param title title string
#' @param wtype = 1 to use dnomr(theta) as the weight \cr
#' = 0 to use no weight.
#' @param maxiter Max # of iterations.
#' @param eps Convergence criterion for optim.
#' @param maxiter2 Max # of inner ietrations.
#' @param method Method for optim: NULL for the default method.
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
#' res <- fitI2L_ls( t(V1), print=1, plot=1 )
#' plot(theta0, res$theta,type="b", main="original theta vs recovered theta")
#'
#' @export
#'

fitI2L_ls <- function( V, print=0, plot=0, title=NULL, wtype=0
                       , maxiter=1000, eps=1e-6, maxiter2=9, method="BFGS"
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
 # optim for ab -> get_ab: 20180116
 # constant w: 20180116
 #

 # V is nitem x npoint
 #


 get_ab <- function( param2, theta, Pj, w
                     , maxiter=100, eps=1e-6, print=0 ){
  # fitting 2PLM to prob
  # Shin-ichi Mayekawa
  # 20180116
  #
  # Args:
  #
  #  param2 vector of the initial value consisting of p1 and p2
  #  theta  vector of theta points
  #  Pj     vector of icc evaluated at theta
  #  w      weight vector
  #

  npoints=length(theta)
  w2=sqrt(w)
  ncatj=2
  paramjdf=data.frame( name="temp", type="B", ncat=2
                       , p1=1, p2=0, p3=0, stringsAsFactors=0 )
  paramjdf[,4:5]=param2
  paramj=param2
  paramjp=paramj
  Pjhat=c( irf( paramjdf, theta, print=0, zero=0 )$ICRF )
  rmsej=sqrt( sum( w*(Pj-Pjhat)^2 )/npoints )
  rmsejp=rmsej

  if( rmsej <= eps ){
   return( param2 )
  }

  if( print > 0 ){
   cat("\nFitting 2PLM to the Given ICC\n")
   Print(param2,theta,Pj)
  }

  # GN iteration
  for( llll in  1:maxiter ){

   # Jacobian:  d vec(Pj) / d paramj    ncatj*npoints x ncatj
   # Jac=JacobianMat( paramj, icrf, ..eps.. = 1e-06 )
   paramjdf[1,3+(1:ncatj)]=paramj
   Jac=w2*dirf_p( paramjdf, theta, print=0, zero=0 )$Jack
   # Print(Jac-Jac0)
   # gradient and info mat and GN direction
   g=-t( t(w*(Pj-Pjhat))%*%Jac )
   maxag=max(abs(g))
   H=t(Jac)%*%Jac
   d=Ginv(H, MASS=1)%*%g

   # Print(Jac,g,d)

   ok=0
   step=1
   for( lllll in 1:20 ){
    pj=paramj - step*d
    paramjdf[1,3+(1:ncatj)]=pj
    Pjhat=c( irf( paramjdf, theta, print=0, zero=0 )$ICRF )
    rmse1=sqrt( sum( w*(Pj-Pjhat)^2 )/npoints )
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
   rmseimp=rmsejp-rmsej
   rmseimpr=rmseimp/rmsejp
   maxadp=max(abs(paramjp-paramj))
   if( print >= 2 )
    Print( llll, lllll, rmsej, rmsejp, rmseimp, rmseimpr, maxadp, maxag
           , fmt="i3 i2 .6")

   if( abs(rmseimpr-1) <= eps || rmseimpr <= eps  ||  maxadp <= eps ) break

   # next iteration
   rmsejp=rmsej
   paramjp=paramj


  } # end of llll loop


  if( print ){
   cat("\nEstimated 2PLM Item Parameter\n")
   Print(paramj)
  }


  return( paramj )

 } # end of get_ab









 rssab <- function( paramj, tV, w ){
  rss=tV-icrfB( matrix(paramj,,3), theta )
  rss=sum( w*rss*rss )
  return( rss )
 } # end of rssab

 rssth <- function( theta, tV, w ){
  rss=tV-icrfB( param, theta )
  rss=sum( w*rss*rss )
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
  rss=rssab( param, tV, w )
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

 # weight
 if( wtype == 1 ){
  w=dnorm(theta)
 }
 else{
  w=rep(1,npoint)
 }
 w=w/sum(w)

 tV=t(V)

 rmsep=sqrt(rssab( param, tV, w )/nitem/npoint)

 if( print ){
  cat("\nFitting 2PLM to LRT Probability by LS\n")
  if( !is.null(title) ) Print(title)
  Print(nitem, npoint, rmsep, wtype)
  Print( maxiter, eps, maxiter2 )
  if( print >= 2 ){
   Print("initial")
   Print(param,theta)
  }
 }

 control=list(reltol=eps, maxit=maxiter2)


 ###################################################################(1)
 # store param
 # initialize inline SQUAREM
 rmv(iSQUAREM)
 pp=c( param[,1:2], theta )
 iSQUAREM=generate_iSQUAREM( pp )
 param0=param
 theta0=theta
 ###################################################################(1)

 Print(param)

 for( llll in 1:maxiter ){

  # item parameters
  for( j in 1:nitem ){

   # paramj=param[j,1:2]
   # res=nlm( rssab, paramj, iterlim=5, gradtol=eps, tV=tV[,j], wtype=wtype )
   # paramj=res$estimate
   # param[j,]=paramj
   # if( is.null(method) ){
   #  res=optim( paramj, rssab, control=control, tV=tV[,j], wtype=wtype )
   # }
   # else{
   #  res=optim( paramj, rssab, control=control, tV=tV[,j], wtype=wtype
   #             , method=method )
   # }
   # paramj=res$par
   # param[j,1:2]=paramj
   # Print(param,paramj)

   paramj=param[j,1:2]
   paramj=c(1,0)
   vj=tV[,j]
   paramj=get_ab( paramj, theta, vj, w, maxiter=50, print=0 )
   param[j,1:2]=paramj

  } # end of j

  rss=rssab( param, tV, w )


  # theta
  if( 1 ){
   # theta points
   # res=nlm( rssth, theta, iterlim=5, gradtol=eps, tV=tV, w=w )
   # theta=res$estimate
   if( is.null(method) ){
    res=optim( theta, rssth, control=control, tV=tV, w=w )
   }
   else{
    res=optim( theta, rssth, control=control, tV=tV, w=w
               , method=method )
   }
   theta=res$par

  }
  #theta=seq(-3,3,length(ncol(V)))

  # normalize
  ms=mands(theta)
  mth=ms[1,1]; sth=ms[2,1]
  theta=(theta-mth)/sth
  param[,2]=(param[,2]-mth)/sth
  param[,1]=param[,1]*sth


  rss=rssab( param, tV, w )
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
 rss=colSums( w*(t(V)-irf)^2 )
 rmse=sqrt(sum(rss/npoint)/nitem)
 rmsej=sqrt(rss/npoint)


 if( print ){
  cat("\nFitting 2PLM to LRT Probability by LS\n")
  if( !is.null(title) ) Print(title)
  Print(nitem, npoint, rmse, wtype)
  Print( maxiter, eps, maxiter2 )
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



etime(0)
res=fitI2L_ls( V2, print=2, plot=1, wtype=0, SQUAREM=0 )
etime(1)





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













