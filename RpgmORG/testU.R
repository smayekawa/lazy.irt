temp=c("
 name type ncat p1 p2 p3 p4 p5 p6
 Q1 G 6  1, -2,-1,0,1,2
 Q2 G 6  1, -3,-1,0,1,2
 Q3 G 6  1, -2,-1,0,1,3
 Q4 G 6  1, -3,-1,0,1,3
 Q5 G 6  1, -2,-1.5,0,1,2
 Q6 G 6  1, -2,-1,0,1.5,2
 Q7 G 6  1, -2,-1.5,0,1.5,2
 Q8 G 6  1, -3,-1.5,0,1,2
 Q9 G 6  1, -3,-1,5,1,2
  ")
 param=cards(temp,header=1)


 set.seed(1701)
 npoints=1000
 nitems=nrow(param)
 ncat=param$ncat
 type=param$type
 Nmat=matrix(sample(1:19,npoints*nitems,replace=1),npoints,nitems)
 pmiss=0.2
 for( j in 1:nitems ){
  Nmat[sample(1:npoints,npoints*pmiss),j]=0
 }
 temp <- gendataIRT( 1, param, Nmat=Nmat, npoints=npoints
                     , thdist="rnorm", compress=0 )
 UN <- as.data.frame(temp$U)
 res3 <- uIRT( U=UN, ncat=ncat, type=type, maxiter=200 )

 temp=irf(res3$param, plot=1, print=0)


 res3 <- uIRT( U=UN, ncat=ncat, maxiter=200 )

 temp=irf(res3$param, plot=1, print=0)
