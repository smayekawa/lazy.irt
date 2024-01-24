#
set.seed(1701)

param=paramB1[c(1:3,7:9,13:15),]
thmin=-2; thmax=2; npoint=5
N=1000

# discrete theta
# theta0=seq(thmin,thmax,length=npoint)
theta0=c(-2, -1, 0, 2, 3)
theta=unlist(lapply( theta0, rep, round(N/npoint) ))
res2 <- gendataIRT( 1, paramB1, theta=theta, compress=1 )
Uc=as.data.frame(res2$U)
ncat=res2$ncat
type=res2$type


# normal theta
Uc2 <- gendataIRT( 1, paramB1, npoints=N, thdist="rnorm", compress=1 )$U
Uc2 <- as.data.frame(Uc2)



comments('

set.seed(1701)

param=paramB1[c(1:3,7:9,13:15),]
thmin=-2; thmax=2; npoint=5
N=1000

# discrete theta
# theta0=seq(thmin,thmax,length=npoint)
theta0=c(-2, -1, 0, 2, 3)
theta=unlist(lapply( theta0, rep, round(N/npoint) ))
res2 <- gendataIRT( 1, paramB1, theta=theta, compress=1 )
Uc=as.data.frame(res2$U)
ncat=res2$ncat
type=res2$type

nclass=5
res1 <- uLRT( Uc, nclass=nclass, estrho=1, monotone=1, alpha=10
              , maxiter=20, plot=1, print=1 )

# normal theta
Uc2 <- gendataIRT( 1, paramB1, npoints=N, thdist="rnorm", compress=1 )$U
Uc2 <- as.data.frame(Uc2)

nclass=5
res1 <- uLRT( Uc2, nclass=nclass, estrho=1, monotone=1
              , maxiter=20, plot=1, print=1 )


')
