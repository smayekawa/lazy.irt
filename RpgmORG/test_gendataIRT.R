set.seed(1701)
temp <- gendataIRT( 1, paramS1, npoints=20, thdist="rnorm", compress=1 )
type=temp$type
ncat=temp$ncat
Uc <- as.data.frame(temp$U)
U=dummy_expand( indata )
Print(Uc)


set.seed(1701)
P=matrix(c(1,2,3, 3,2,1, 1,1,1), 3,3, byrow=1)
U=matrix(0,3,3)
for( i in 1:nrow(P) ){
 U[i,]=t( rmultinom( 1, 100, P[i,] ) )
}
Print(U)





set.seed(1701)
npoints=20
param=paramS1
nitems=nrow(param)
temp <- gendataIRT( 1, param, npoints=npoints, thdist="rnorm", compress=1 )
type=temp$type
ncat=temp$ncat
Uc <- as.data.frame(temp$U)
Print(Uc)


set.seed(1701)
npoints=20
param=paramS1[1:2,]
param[1,]$p2=-1
nitems=nrow(param)
Nmat=matrix(1000,npoints,nitems)
Nmat[10,]=0
Nmat[11,1]=0
temp <- gendataIRT( 1, param, Nmat=Nmat, npoints=npoints
                    , thdist="rnorm", compress=1 )
type=temp$type
ncat=temp$ncat
UcN <- as.data.frame(temp$U)
Print(UcN)






