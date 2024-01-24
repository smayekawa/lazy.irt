smn <- function( x, f, mu=NULL, sigma=NULL,    estmu=1, estsigma=1
               , maxiter=100, eps=1e-8, print=1, plot=0 ){
 # parametric scored multinomial: normal pdf model
 # Shin-ichi Mayekawa
 # 20121026
 # comments: 20121110(London)
 # plot: 121207
 # plot bugfix: 20150614
 #



 # const
 if( is.vector(x) ) m=length(x) else  m=nrow(x)
 n=sum(f)
 meanx=sum(x*f)/n
 stdx=sqrt( sum(x*x*f)/n-meanx^2 )
 P0=exp(-0.5*(x-meanx)^2/stdx^2)
 P0=P0/sum(P0)
 llh0=sum(f*log(P0))


 # initial values
 if( is.null(mu) )  mu=meanx
 if( is.null(sigma) ) sigma=stdx

 P=exp(-0.5*(x-mu)^2/sigma^2)
 P=P/sum(P)

 if( print > 0 ){
  cat("\n\nParametric Scored Multinomial Distribution: Normal PDF model\n\n")
  cat(" # of observations =", n,"\n")
  cat(" # of categories =", m,"\n")
  cat(" estimation of mu =", estmu, ",  estimation of sigma =", estsigma,"\n")
  cat(" max # of iterations =", maxiter," with eps =", eps,"\n")
  cat(" sample mean and std\n")
  cat("  mean =", meanx, ",  std =", stdx,"\n")
  cat("  likelihood with the sample mean/std =", llh0,"\n")
 }


 llhp=llh0
 converged=0

 # main iteration
 for( llll in 1:maxiter ){

  # log-Jacobian matrix
  Jac=cbind( (x-mu)/sigma^2, (x-mu)^2/sigma^3 )

  # dispersion
  Dr=n*( Diag(P)-P%*%t(P) )

  # gradient etc
  delta=matrix(t((f-n*P)),,1)
  g=t(Jac)%*%delta
  H=t(Jac)%*%Dr%*%Jac
  d=solve(H)%*%g

  # update parameters
  step=1; ok=0
  mu1=mu; sigma1=sigma;
  for( lllll in 1:20 ){
   if( estmu ) mu1=mu + step*d[1]
   if( estsigma ) sigma1=sigma + step*d[2]
   P1=exp(-0.5*(x-mu1)^2/sigma1^2)
   P1=P1/sum(P1)
   llh=sum(f*log(P1))
   if( llh >= llhp ){
    ok=1
    break
   }
   step=step/2
  }
  if( ok ){
   mu=mu1; sigma=sigma1
   P=P1
  }

  # convergence
  llhimpr=(llh-llhp)/abs(llhp)
  maxag=max(abs(c(estmu,estsigma)*g))
  if( print >= 2 )
   Print(llll, llhp, llh, llhimpr,maxag)
  if( llhimpr <= eps  &&  maxag < eps ){
   converged=1
   break
  }

  # next iteration
  llhp=llh


 } # end of llll loop

 mean=sum(x*P); std=sqrt(sum(x*x*P)-mean^2)
 iter=llll

 # estimated density etc
 nP=n*P
 tab=cbind(x,f,nP,P)
 colnames(tab)=c("x","f","n*p","p")

 if( print >= 1 && converged )
  cat("\nIteration converged with ", iter, " iterations.  eps= ",eps,"\n")

 if( print > 0 ){
  cat("\n\nParametric Scored Multinomial Distribution: Normal PDF model\n\n")
  cat(" # of observations =", n,"\n")
  cat(" # of categories =", m,"\n")
  cat(" estimation of mu =", estmu, ",  estimation of sigma =", estsigma,"\n")
  cat("\n iterations required =", iter," with eps =", eps,"\n")
  cat(" log likelihood maximized =", llh, "\n")
  cat("\n estimated parameters of the Normal distribution\n")
  cat("  mu =", mu, ",  sigma =", sigma,"\n")
  cat("\n sample mean and std\n")
  cat("  mean =", meanx, ",  std =", stdx,"\n")
  cat(" mean and std of the scored multinomial (preserved)\n")
  cat("  mean =", mean, ",  std =", std,"\n")
  cat("  likelihood with the sample mean/std =", llh0, "\n")
  cat("\n final gradient\n")
  cat("  mu =", g[1], ",  sigma = ", g[2], "\n")
  if( print >= 3 ){
   cat("\ndata and the expected values (= n*prob)\n")
   print(tab)
  }
 }

 if( plot ){
  P=exp(-0.5*(x-mu1)^2/sigma1^2)
  sumP=sum(P)
  P=P/sum(P)
  x1=seq(min(x),max(x),length.out=51)
  P1=exp(-0.5*(x1-mu1)^2/sigma1^2)
  P1=P1/sumP
  nP1=n*P1
  maxy=max(nP,nP1,f)+5
  temp=barplot(height=c(f), names.arg=c(x), ylim=c(0,maxy)
      , main="Histogram of (x,f) with the Estimated Expected Values"
      , sub=paste("mu=",format(mu,digits=4)
                , ", sigma=", format(sigma,digits=4), "\n"
                , paste("  (mean=",format(meanx,digits=4)
                      , ", std=",format(stdx,digits=4), ")") )
                , offset=0  )
  par(new=1)
  plot( x, nP, ylim=c(0,maxy), xlim=c(min(x)-.5,max(x)+.5), type="p"
        , xlab="", ylab="", xaxt="n", yaxt="n" )
  par(new=1)
  plot( x1, nP1, ylim=c(0,maxy),  xlim=c(min(x)-.5,max(x)+.5), type="l"
        , xlab="", ylab="freq and n*prob", xaxt="n", yaxt="n" )
 }


 return(
     list( P=P, x=x, mu=mu, sigma=sigma, estmu=estmu, estsigma=estsigma
         , llh=llh, mean=mean, std=std )
        )


} # end of smn







seed <- 1701
set.seed(seed)
n <- 1000
npoints <- 21
resg <- freqdist( x=rnorm(n), weights=1, npoints=npoints, min=-3, max=3 )
res <- smn( resg$midpoints, resg$freq, print=1, plot=1, mu=0, sigma=1 )
