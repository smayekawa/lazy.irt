library(lazy.irt)

pdata=read.csv( "RpgmORG/IRT_9type.csv", header=TRUE ,row.names=1)

#number of picher type
nitem=9

#response data
s=6
ncat=rep(s,nitem)

# new varaible names
cn0=paste("P",1:nitem,"_0",sep="")
cn1=paste("P",1:nitem,"_1",sep="")
cn2=paste("P",1:nitem,"_2",sep="")
cn3=paste("P",1:nitem,"_3",sep="")
cn4=paste("P",1:nitem,"_4",sep="")
cn5=paste("P",1:nitem,"_5",sep="")
cn=c( mapply(c,cn0,cn1,cn2,cn3,cn4,cn5) )
colnames(pdata)=cn

type=rep("G",nitem)
type=rep("P",nitem)

res=uIRT( U=pdata, ncat=ncat, type=type, print=2, minp1=0.1, maxabsparam=20 )
temp=irf( res$param,print=0,plot=1 )

mode=find_mode( res$param, print=1 )


theta=matrix(esttheta( U=pdata/1.3, param=res$param,
                       npoints=7 )$thetahat,,1)
rownames(theta)=rownames(pdata)
( locnan=which( is.nan(theta) ) )



