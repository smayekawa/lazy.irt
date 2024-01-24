
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
   g=-t( t(w2*(Pj-Pjhat))%*%Jac )
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



npoints=21
theta=seq(-3,3,length=npoints)
w=dnorm(theta)
w=w/sum(w)
param=paramB1
irf=irf( param, theta, plot=0, print=0, zero=0 )

itemid=1
Pj=irf$ICRF[,itemid]
param2=c(1,0)

res1=get_ab( param2, theta, Pj, w, print=2 )





comments(
 '

set.seed(1701)

n=100
p1=0.6+runif(n)
p2=rnorm(n)
p3=0.5*runif(n)

paramB3=data.frame( name=paste("Q",1:n,sep=""), type="B3", ncat=2,
                    p1=p1, p2=p2, p3=p3, stringsAsFactors=0 )

res2 <- fit223_ls( paramB3, plot=1, print=0, wtype=1 )

mean(res2$rmse)


 ')
