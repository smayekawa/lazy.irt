seed=170100
set.seed(seed)

nitem=300
n=5000
iname=paste("Q",formatC(1:nitem,width=3,format="d",flag=0),sep="")
b=sample(c(-1,0,1),nitem,replace=1)
b=sort(b)

paramL=data.frame(name=iname, type="B2", ncat=2, p1=0.8, p2=b, p3=0 )


dataIRT=gendataIRT( Ntot=1, npoints=n, param=paramL, thdist="rnorm", zero=0
                 , compress=0, sort=1 )
Uc=data.frame(id=paste("s",formatC(1:n,width=4,format="d",flag=0),sep=""), dataIRT$U, stringsAsFactors = 0 )

res=uIRT( Uc=Uc, idvar="id" )
