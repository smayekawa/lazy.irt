param=paramB1[1:5,]
maxscore=sum(param$ncat-1)
param$p1=1
res <- obscore( param, npoints=51 )
theta_stat=res$theta_stat
theta=theta_stat$theta
TRF=theta_stat$TRF
stdx_TRF=theta_stat$stdx_t

plot(TRF,stdx_TRF)
plot(theta,stdx_TRF)

# Transform x-t scale to g(x)-g(t).
# Then,  the variance of g(x) is g'(t)^2 stdx_t(t)^2.
# We want above to be a constant, sigma^2.
# Therefore, g'(x) = sigma/stdx_t(x)
#

sigma=0.5



gdash <- function( t ){
 res=interpol( TRF,gd0, t )[,2]
 return( res )
} # gdash



gd0=sigma/stdx_TRF
t=seq(0,maxscore,0.1)
stdx_t=interpol( TRF,stdx_TRF, t )[,2]
g=Integrate( gdash, t, 0 )
plot(t,g)

gd=gdash( t )
gs=gd*stdx_t
plot( g, gs )
Print( t,stdx_t, g, gd, gs )
