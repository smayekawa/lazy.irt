





pp=rbind(paramB1,paramB1,paramB1)
pp$name=1:nrow(pp)
# pp$p1=.5

res=gendataIRT( 1, pp, npoints=10000, zero=0, thdist="rnorm" )
U=res$U
X=rowSums(U)
hist(X, nclass=55)


ppp=pp
ppp$p1=.5

res=gendataIRT( 1, ppp, zero=0, npoints=10000,thmin=-3, thmax=3, thdist="rnorm" )
U=res$U
X=rowSums(U)
hist(X, nclass=55)



#  out_obscore=obscore( ppp, weight=NULL, npoints=121, print=0, plot=3 )





ppp=data.frame(name=paste("Q",1:10,sep=""),type="B", ncat=2, p1=1,p2=0,p3=0
               , stringsAsFactors = 0)





theta=res$theta
icrf=res$icrf
icrf=icrf[,grep("_1",colnames(icrf))]

Print(icrf[,1],logistic(theta-ppp[1,"p2"]),logistic1(theta-ppp[1,"p2"]), fmt="4.2")
