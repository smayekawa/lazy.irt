
a=1

thmax=8
thmin=-thmax
theta=seq(thmin,thmax,length=25)
num=a*dnorm(a*theta)
den=1-pnorm(a*theta)
v=num/den
Print(theta,num,den,v, fmt="15.9")



param=paramS2[5:6,]
param$type="Gn"

resInfo=info_func( param, plot=4, print=3, theta=theta, numderiv=0, smallP=0 )


param=paramS1[2,]
param$p1=1
param$p3=0.01
param$type="Bn"

resInfo=info_func( param, plot=4, print=3, theta=theta, numderiv=0, smallP=0 )

