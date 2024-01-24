generate_prob <- function( param, theta ){
 # generate ntt probs
 # Shin-ichi Mayekawa
 # 20170316,17
 #

 V=irf( param, theta, print=0, plot=0 )$IRF

 return( t(V) )

} # generate_prob



seed=1701+5
set.seed(seed)

nth=5
# equally spaced theta
th=seq(-2,2, length=nth)
# modify theta
th[1]=-1.25
th[3]=0.8

# generate prob
V=generate_prob( paramB1, th )
# add error
error=0.2
V=V+error*(runif(nth)-0.5)
V[V>=1]=1-0.05
V[V<=0]=0.05

