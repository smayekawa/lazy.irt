# comparing sumsmn to normal approx
# 20180703
#

# score range
min=0
max=1

x=min:max
ncat=length(x)

# rough mean/std
mean0=(max+min)/2
std0=(max-min)/2/2

# # of multinomial components to be added
npoints=2

# mean and std of each component
mm=gen_midp( min, max, npoints+2 )[-c(1,npoints+2)]
ss=rep(std0/sqrt(npoints), npoints)

Print( min, max, ncat, npoints, mean0, std0 )
Print( x, mm, ss )

# domain and prob
V=matrix(x,ncat,npoints)
P=mapply( function( m,s ) graded_prob(x,m,s, method=1), mm, ss )
Print(colSums(P))


# same weight for each component
w=rep(1,npoints)

Print(V, P, w)

# exact result
res=sumsmnw( P, V, w, compress=0, print=1, plot=1 )
ms=mands( res[,1], res[,2])
Print(ms)


