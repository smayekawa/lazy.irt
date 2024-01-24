# sum of uniforms
# 20170615,16,0831
#
seed=1701
set.seed(seed)

# # of scores
n=500
# # of subtests
r=5

# X has equi-correlation matrix
rho=0.5
R=(1-rho)*diag(r)+rho*matrix(1,r,r)
udv=svd(R)
RR=udv$u%*%diag(sqrt(udv$d))
X=matrix(rnorm(n*r),n)
X=X%*%t(RR)
X=X%*%diag(1/sqrt(diag(var(X))))

# percent score
Y=pnorm(X)

Print(cor(X),cor(Y), fmt="4.2")

# graded percent score based on Y
brks=seq(0,1,0.1)
Yc=apply( Y,2, cut, breaks=brks, labels=0:9 )
Ycn=matrix(as.numeric(Yc),n)

Print(cor(X),cor(Ycn), fmt="4.2")

# histogram
Ys=rowSums(Ycn)
midp=0:9
brks=midp2brk( midp )
apply( Ycn, 2, hist, breaks=brks
       , main=paste("n=",n, ",  # of subtests=",r,",  corr=",rho,sep="") )

# histogram of sum of Ycn
hist(Ys, main=paste("n=",n, ",  sum of ",r," subtests.,  corr=",rho,sep=""))

