# sum of uniforms
#
#

n=500
r=5
X=matrix(sample(0:9,n*r,replace=1),n)
Xs=rowSums(X)
midp=0:9
brks=midp2brk( midp )
apply( X, 2, hist, breaks=brks )
hist(Xs)

