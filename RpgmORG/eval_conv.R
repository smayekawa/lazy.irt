
esttheta <- function( U, P ){
 logP=log(P)
 ULP=U%*%t(logP)
 ULP[ULP > 709]=709
 eULP=exp(ULP)*matrix(1,nrow(ULP))%*%t(thd)
 rseULP=rowSums(eULP)
 rseULP[rseULP < 1e-307]=1e-307
 H=eULP/rseULP
 thetahat=H%*%theta
 return( as.vector(thetahat) )
} # end of esttheta


seed=1701
set.seed(seed)


ncat=4
param1=data.frame( name="Q1", type="P", ncat=ncat, stringsAsFactors=0 )
p=matrix(0,1,ncat)
colnames(p)=paste("p",1:ncat,sep="")
param1=cbind(param1,p)

npoint=51
theta=seq(-4,4,length=npoint)
thd=dnorm(theta)
thd=thd/sum(thd)

U=diag(ncat)

nrep=1


stat=NULL
sname=c("info","rmse_p","rmse_ii","rmse_iic","rmse_t")


for( llll in 1:nrep ){

 # generate param
 a=1
 b=rnorm(ncat-1)
 b=sort(b)
 param1[,3+(1:ncat)]=c(a,b)

 # original stat
 res=irf( param1, theta, print=0, plot=1 )
 P1=res$ICRF
 thetahat1=esttheta( U, P1 )

 Print(llll, thetahat1)


 # conversion: P to G

 # info=0
 res2=fitG2P_ls( param, theta, wtype=1, info=0 )
 param2=res2$paramG
 #res=irf( param2, theta, print=0, plot=0 )
 P2=res2$icrfG
 thetahat2=esttheta( U, P2 )

 # stat
 rmse_t=sqrt(sum( (thetahat1-thetahat2)^2 )/npoint )

 stat0=matrix(
  c( res2$info, res2$rmse_p, res2$rmse_ii, res2$rmse_iic, rmse_t )
  , 1, )
 stat=rbind(stat,stat0)


 # info=1
 res2=fitG2P_ls( param, theta, wtype=1, info=1 )
 param2=res2$paramG
 #res=irf( param2, theta, print=0, plot=0 )
 P2=res2$icrfG
 thetahat2=esttheta( U, P2 )

 # stat
 rmse_t=sqrt(sum( (thetahat1-thetahat2)^2 )/npoint )

 stat0=matrix(
  c( res2$info, res2$rmse_p, res2$rmse_ii, res2$rmse_iic, rmse_t )
  , 1, )
 stat=rbind(stat,stat0)


 # info=2
 res2=fitG2P_ls( param, theta, wtype=1, info=2 )
 param2=res2$paramG
 #res=irf( param2, theta, print=0, plot=0 )
 P2=res2$icrfG
 thetahat2=esttheta( U, P2 )

 # stat
 rmse_t=sqrt(sum( (thetahat1-thetahat2)^2 )/npoint )

 stat0=matrix(
  c( res2$info, res2$rmse_p, res2$rmse_ii, res2$rmse_iic, rmse_t )
  , 1, )
 stat=rbind(stat,stat0)




} # end of llll

colnames(stat)=sname









