#
set.seed(1701)

param=paramB1[c(1:3,7:9,13:15),]

thmin=-2;
thmax=2
npoint=5
N=1000

theta0=seq(thmin,thmax,length=npoint)
theta0=c(-2, -1, 0, 2, 3)
theta=unlist(lapply( theta0, rep, round(N/npoint) ))

res2 <- gendataIRT( 1, paramB1, theta=theta )
# Print(res2$N,res2$U,res2$theta,rowSums(res2$U[,seq(2,ncol(res2$U),2)]))

U=res2$U
ncat=res2$U
type=res2$type
Uc=as.data.frame( dummy_compress(U)$Uc )



# normal theta
Uc2 <- gendataIRT( 1, paramB1, npoints=1000, thdist="rnorm", compress=1 )$U
Uc2 <- as.data.frame(Uc2)
