# sum of uniforms
#
#
seed=1701
set.seed(seed)

n=500
r=3

rho=0.5
R=(1-rho)*diag(r)+rho*matrix(1,r,r)
udv=svd(R)
RR=udv$u%*%diag(sqrt(udv$d))
X=matrix(rnorm(n*r),n)
X=X%*%t(RR)
X=X%*%diag(1/sqrt(diag(var(X))))
cor(X)

Y=pnorm(X)
cor(Y)

brks=seq(0,1,0.1)
Yc=apply( Y,2, cut, breaks=brks, labels=0:9 )
Yn=matrix(as.numeric(Yc),n)

Print(cor(X),cor(Yn), fmt="4.2")


Ys=rowSums(Yn)
midp=0:9
brks=midp2brk( midp )
apply( Yn, 2, hist, breaks=brks
       , main=paste("n=",n, ",  # of subtests=",r,",  corr=",rho,sep="") )
hist(Ys, main=paste("n=",n, ",  # of subtests=",r,",  corr=",rho,sep=""))

