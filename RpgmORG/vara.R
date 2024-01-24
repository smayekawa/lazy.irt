set.seed(1701)

p=c(1,2,2)
p=p/sum(p)
ncat=length(p)

repmax=100

res=matrix(0,repmax,3)
resa=matrix(0,repmax,2*ncat)
for( i in 1:repmax ){

 p=runif(ncat)
 p=p/sum(p)

 a=sample(1:(ncat+1),ncat-1,replace=0)
 a=c(0,a)
 a=0:(ncat-1)
 ms=mands( a, p )
 mean=ms[1]
 var=ms[3]
 dd=var-mean^2
 res[i,]=c(var,mean^2,dd)
 resa[i,]=c(p,a)
 Print(i,p,a,dd)


}

locpos=which(res[,3] >= 0)
Print(res[locpos,],resa[locpos,])

locneg=which(res[,3] < 0)
Print(res[locneg,],resa[locneg,])



