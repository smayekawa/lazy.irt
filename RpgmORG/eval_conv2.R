
esttheta <- function( U, P ){
 logP=log(P)
 ULP=U%*%t(logP)
 # ULP[ULP > 709]=709
 eULP=exp(ULP)*matrix(1,nrow(ULP))%*%t(thd)
 rseULP=rowSums(eULP)
 # rseULP[rseULP < 1e-307]=1e-307
 H=eULP/rseULP
 thetahat=H%*%theta
 return( as.vector(thetahat) )
} # end of esttheta


fitG <- function( info, print=1 ){

 # convert
 res2=fitG2P_ls( param1, theta, wtype=1, info=info, print=print )
 param2=res2$paramG

 # thetahat
 P2=res2$icrfG
 thetahat2=esttheta( U, P2 )
 rmse_t=sqrt( sum( (thetahat1-thetahat2)^2 )/ncat )
 rmse_tt=sqrt( sum( ((thetahat1-thetahat2)^2)[-c(1,ncat)] )/(ncat-1) )

 # accumulate
 stat0=c( info, unlist(param2[1,3+(1:ncat)])
          , res2$rmse_p, res2$rmse_ii, res2$rmse_iic
          , rmse_t, rmse_tt, thetahat2 )

 return( stat0 )

} # end of fitG


fitGn <- function( info, print=1 ){

 # convert
 res2=fitGn2P_ls( param1, theta, wtype=1, info=info, print=print )
 param2=res2$paramG

 # thetahat
 P2=res2$icrfG
 thetahat2=esttheta( U, P2 )
 rmse_t=sqrt( sum( (thetahat1-thetahat2)^2 )/ncat )
 rmse_tt=sqrt( sum( ((thetahat1-thetahat2)^2)[-c(1,ncat)] )/(ncat-1) )

 # accumulate
 stat0=c( info, unlist(param2[1,3+(1:ncat)])
          , res2$rmse_p, res2$rmse_ii, res2$rmse_iic
          , rmse_t, rmse_tt, thetahat2 )

 return( stat0 )

} # end of fitGn



fitP <- function( info, print=1 ){

 # convert
 res2=fitP2G_ls( param1, theta, wtype=1, info=info, print=print )
 param2=res2$paramP

 # thetahat
 P2=res2$icrfP
 thetahat2=esttheta( U, P2 )
 rmse_t=sqrt( sum( (thetahat1-thetahat2)^2 )/ncat )
 rmse_tt=sqrt( sum( ((thetahat1-thetahat2)^2)[-c(1,ncat)] )/(ncat-1) )

 # accumulate
 stat0=c( info, unlist(param2[1,3+(1:ncat)])
          , res2$rmse_p, res2$rmse_ii, res2$rmse_iic
          , rmse_t, rmse_tt, thetahat2 )

 return( stat0 )

} # end of fitGn



seed=1701
set.seed(seed)


mode="Gn2P"  #  G2P, Gn2P, P2G, P2Gn
ncat=4

param1=data.frame( name="Q1", type="x", ncat=ncat, stringsAsFactors=0 )
p=matrix(0,1,ncat)
colnames(p)=paste("p",1:ncat,sep="")
param1=cbind(param1,p)

# which model
{
 if( mode == "G2P" ){
  fitfunc=fitG
  param1$type="P"
 } else if( mode == "Gn2P" ){
  param1$type="P"
  fitfunc=fitGn
 }
 else if( mode == "P2G" ){
  param1$type="G"
  fitfunc=fitP
 }
 else if( mode == "P2Gn" ){
  param1$type="Gn"
  fitfunc=fitP
 }
}

npoint=51
theta=seq(-4,4,length=npoint)
thd=dnorm(theta)
thd=thd/sum(thd)

U=diag(ncat)

nrep=50

sname0=c( "ncat", paste("p1",1:ncat,sep=""), paste("th1",0:(ncat-1),sep="") )
sname=c( "info", paste("p2",1:ncat,sep="") )
sname=c( sname, "rmse_p","rmse_ii","rmse_iic","rmse_t","rmse_tt" )
sname=c( sname, paste("th2",0:(ncat-1) ,sep="") )
stat=matrix(0,nrep*3,length(sname0)+length(sname))

for( llll in 1:nrep ){

 # generate param
 a=1
 b=rnorm(ncat-1)
 b=sort(b)
 param1[,3+(1:ncat)]=c(a,b)


 Print(llll, param1)


 # original stat
 res=irf( param1, theta, print=0, plot=0 )
 P1=res$ICRF
 thetahat1=esttheta( U, P1 )

 ss=c( ncat, unlist(param1[1,3+(1:ncat)]), thetahat1 )


 # conversion: P to G

 # info=0
 stat0=fitfunc( info=0, print=0 )
 stat[3*(llll-1)+1,1:(2*ncat+1)]=ss
 stat[3*(llll-1)+1,(length(ss)+1):(length(ss)+length(stat0))]=stat0

 # info=1
 stat0=fitfunc( info=1, print=0 )
 stat[3*(llll-1)+2,1:(2*ncat+1)]=ss
 stat[3*(llll-1)+2,(length(ss)+1):(length(ss)+length(stat0))]=stat0

 # info=2
 stat0=fitfunc( info=2, print=0 )
 stat[3*(llll-1)+3,1:(2*ncat+1)]=ss
 stat[3*(llll-1)+3,(length(ss)+1):(length(ss)+length(stat0))]=stat0


} # end of llll

colnames(stat)=c(sname0,sname)


statdf=cbind(data.frame(mode=mode, seed=seed, stringsAsFactors=0),stat)
Print(statdf)

title=paste(mode, ":  rmse_")
boxplot( rmse_t~info, data=statdf, main=paste(title,"theta",sep="") )
boxplot( rmse_tt~info, data=stat, main=paste(title,"theta trimmed",sep="") )
boxplot( rmse_p~info, data=stat, main=paste(title,"prob",sep="") )
boxplot( rmse_ii~info, data=stat, main=paste(title,"item info",sep="") )
boxplot( rmse_iic~info, data=stat, main=paste(title,"item cat info",sep="") )




locbad=which(stat[,"rmse_ii"] > 1)

statbad=statdf[locbad,]
save(statbad, file="RpgmORG/statbad.RData")
