Uc=matrix(sample(0:1,50,replace=1),25)
Uc=cbind(Uc,matrix(sample(0:2,50,replace=1),25))
table(Uc)

res=table(Uc[,1],Uc[,2],Uc[,3],Uc[,4])

# freqdist( Uc, npoints=c(2,2,3,4)) # no good

freqdist( Uc, midpoints=list(c(0,1),c(0,1),c(0,1,2),c(0,1,2,3)))

restx=xtabs(~Uc[,1]+Uc[,2]+Uc[,3]+Uc[,4])
restx=xtabs(~Uc[,1]+Uc[,2]+Uc[,3])
f=restx
namesf=dimnames(f) # This is the formatted value of c by xtabs.
dimf=dim(f)
midpindex=as.matrix( expand.grid(lapply(dimf,function(a) seq(1:a)-1)) )
dim(f)=prod(dimf)
f=unclass(as.matrix(f)) # removing the table attribute
Print(midpindex,f)



code=paste("Uc[,",1:4,"]",sep="")
code=paste("f=table(",paste(code,collapse=","),")", collapse="")
eval(parse(text=code))


code=paste("Uc$",1:4,sep="")
code=paste("f=table(",paste(code,collapse=","),")", collapse="")
eval(parse(text=code))


# This fails!!!!!!!!!!!!!!!!!
dd=data.frame( matrix(sample(0:1,10000,replace=1),,50) )
resf=freqdist( dd, midpoints=c(0,1), omitmiss=1 )



dd=data.frame( matrix(sample(0:1,100000000,replace=1),,500) )
ddu=unique(dd)
resf=freqdist( dd, midpoints=c(0,1), omitmiss=1 )
