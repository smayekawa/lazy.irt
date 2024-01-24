

doit <- function( mode, outroot="RpgmORG/stat/", ncat=3
                  , nrep=50, npoint=51, seed=1701 ){

 # evaluate plytomous conversion
 # 20180215,17,18,21
 #

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


 fitG <- function( method, print=1, init=1 ){

  # convert
  res2=fitG2P_ls( param1, theta, init=init, method=method, print=print )
  param2=res2$paramG

  # thetahat
  P2=res2$icrfG
  thetahat2=esttheta( U, P2 )
  rmse_t=sqrt( sum( (thetahat1-thetahat2)^2 )/ncat )
  rmse_tt=sqrt( sum( ((thetahat1-thetahat2)^2)[-c(1,ncat)] )/(ncat-1) )

  resInfo=info_func( param2, theta, plot=0, print=0 )
  a_iinfo2=resInfo$a_info_item_LO

  # accumulate
  stat0=c( method, unlist(param2[1,3+(1:ncat)])
           , res2$rmse_p, res2$rmse_ii, res2$rmse_iic
           , rmse_t, rmse_tt, thetahat2, a_iinfo2 )
  names(stat0)=sname2

  return( stat0 )

 } # end of fitG


 fitGn <- function( method, print=1, init=1 ){

  # convert
  res2=fitGn2P_ls( param1, theta, init=init, method=method, print=print )
  param2=res2$paramG

  # thetahat
  P2=res2$icrfG
  thetahat2=esttheta( U, P2 )
  rmse_t=sqrt( sum( (thetahat1-thetahat2)^2 )/ncat )
  rmse_tt=sqrt( sum( ((thetahat1-thetahat2)^2)[-c(1,ncat)] )/(ncat-1) )

  resInfo=info_func( param2, theta, plot=0, print=0 )
  a_iinfo2=resInfo$a_info_item_LO

  # accumulate
  stat0=c( method, unlist(param2[1,3+(1:ncat)])
           , res2$rmse_p, res2$rmse_ii, res2$rmse_iic
           , rmse_t, rmse_tt, thetahat2, a_iinfo2 )
  names(stat0)=sname2

  return( stat0 )

 } # end of fitGn



 fitP <- function( method, print=1, init=1 ){

  # convert
  res2=fitP2G_ls( param1, theta, init=init, method=method, print=print )
  param2=res2$paramP

  # thetahat
  P2=res2$icrfP
  thetahat2=esttheta( U, P2 )
  rmse_t=sqrt( sum( (thetahat1-thetahat2)^2 )/ncat )
  rmse_tt=sqrt( sum( ((thetahat1-thetahat2)^2)[-c(1,ncat)] )/(ncat-1) )

  resInfo=info_func( param2, theta, plot=0, print=0 )
  a_iinfo2=resInfo$a_info_item_LO

  # accumulate
  stat0=c( method, unlist(param2[1,3+(1:ncat)])
           , res2$rmse_p, res2$rmse_ii, res2$rmse_iic
           , rmse_t, rmse_tt, thetahat2, a_iinfo2 )
  names(stat0)=sname2

  return( stat0 )

 } # end of fitGn



 set.seed(seed)



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

 theta=seq(-4,4,length=npoint)
 thd=dnorm(theta)
 thd=thd/sum(thd)

 U=diag(ncat)


 sname0=c( "ncat", paste("p1",1:ncat,sep=""), paste("th1",0:(ncat-1),sep="")
           , "a_iinfo1" )
 sname2=c( "method", paste("p2",1:ncat,sep="") )
 sname2=c( sname2, "rmse_p","rmse_ii","rmse_iic","rmse_t","rmse_tt" )
 sname2=c( sname2, paste("th2",0:(ncat-1) ,sep=""), "a_iinfo2" )
 sname=c(sname0,sname2)
 stat=matrix(0,nrep*3,length(sname))

 for( llll in 1:nrep ){

  # generate param
  a=1
  b=rnorm(ncat-1)
  b=sort(b)
  param1[,3+(1:ncat)]=c(a,b)


  Print(llll, mode, param1)


  # original stat
  res=irf( param1, theta, print=0, plot=0 )
  P1=res$ICRF
  thetahat1=esttheta( U, P1 )
  resInfo=info_func( param1, theta, plot=0, print=0 )
  a_iinfo1=resInfo$a_info_item_LO
  ss=c( ncat, unlist(param1[1,3+(1:ncat)]), thetahat1, a_iinfo1 )


  # conversion: P to G

  # method=0
  stat0=fitfunc( method=0, print=0 )
  if( max(abs(stat0["rmse_p"])) > 0.1 ){
   stat0=fitfunc( method=0, print=0, init=0 )
  }
  stat[3*(llll-1)+1,1:(2*ncat+1+1)]=ss
  stat[3*(llll-1)+1,(length(ss)+1):(length(ss)+length(stat0))]=stat0

  # method=1
  stat0=fitfunc( method=1, print=0 )
  if( substr(mode,1,1) == "P"  &&  max(abs(stat0["rmse_ii"])) > 0.1 ){
   stat0=fitfunc( method=1, print=0, init=0 )
  }

  stat[3*(llll-1)+2,1:(2*ncat+1+1)]=ss
  stat[3*(llll-1)+2,(length(ss)+1):(length(ss)+length(stat0))]=stat0

  # method=2
  stat0=fitfunc( method=2, print=0 )
  if( max(abs(stat0["rmse_iic"])) > 0.1 ){
   stat0=fitfunc( method=2, print=0, init=0 )
  }
  stat[3*(llll-1)+3,1:(2*ncat+1+1)]=ss
  stat[3*(llll-1)+3,(length(ss)+1):(length(ss)+length(stat0))]=stat0


 } # end of llll

 colnames(stat)=sname

 difai=stat[,"a_iinfo1"]-stat[,"a_iinfo2"]
 loct1=which("th10" == colnames(stat))
 loct2=which("th20" == colnames(stat))
 dift=stat[,loct1-1+(1:ncat)]-stat[,loct2-1+(1:ncat)]
 stat=cbind(stat,difai,dift)
 colnames(stat)=c(sname, "difai", paste("dth",0:(ncat-1),sep=""))


 statdf=cbind(data.frame(mode=mode, seed=seed, stringsAsFactors=0),stat)

 locs=c(loct1:(loct1+ncat-1),(loct2-5):(loct2+2*ncat-1))
 Print(colMeans(stat[,locs]))

 title=paste(mode, ":  rmse_")
 maxt=max(statdf$rmse_t,statdf$rmse_tt)
 maxi=max(statdf$rmse_ii,statdf$rmse_iic)
 boxplot( rmse_t~method, data=statdf, main=paste(title,"theta",sep="")
          , ylim=c(0,maxt) )
 boxplot( rmse_tt~method, data=stat, main=paste(title,"theta trimmed",sep="")
          , ylim=c(0,maxt) )
 boxplot( rmse_p~method, data=stat, main=paste(title,"prob",sep="") )
 boxplot( rmse_ii~method, data=stat, main=paste(title,"item info",sep="")
          , ylim=c(0,maxi) )
 boxplot( rmse_iic~method, data=stat, main=paste(title,"item cat info",sep="")
          , ylim=c(0,maxi) )


 # outfile=paste(outroot, "stat_", mode, ".csv", sep="" )
 # write.table( statdf, outfile, row.names=FALSE, sep=", " )

 return( statdf )

} # end of doit




outroot="RpgmORG/stat/"
mode="G2P"  #  G2P, Gn2P, P2G, P2Gn
ncat=3
nrep=100

seed=1702

res1=doit( mode="G2P", outroot=outroot, ncat=ncat,  nrep=nrep, seed=seed )
res2=doit( mode="Gn2P", outroot=outroot, ncat=ncat,  nrep=nrep, seed=seed )
res3=doit( mode="P2G", outroot=outroot, ncat=ncat,  nrep=nrep, seed=seed )
res4=doit( mode="P2Gn", outroot=outroot, ncat=ncat,  nrep=nrep, seed=seed )

statdf=rbind(res1,res2,res3,res4)

outfile=paste(outroot, "stat_all_",  seed, "_", ncat, ".csv", sep="" )
write.table( statdf, outfile, row.names=FALSE, sep=", " )

savefile=paste(outroot, "stat_all_", seed, "_", ncat, ".Rdata", sep="" )
save(statdf,file=savefile)

try( detach(statdf) )
attach(statdf)

boxplot( rmse_p~mode+method, data=statdf, main="rmse_p" )
boxplot( rmse_ii~mode+method, data=statdf, main="rmse_ii" )
boxplot( rmse_t~mode+method, data=statdf, main="rmse_t" )







outroot="RpgmORG/stat/"
mode="G2P"  #  G2P, Gn2P, P2G, P2Gn
ncat=5
nrep=300

seed=1701

res1=doit( mode="G2P", outroot=outroot, ncat=ncat,  nrep=nrep, seed=seed )
res2=doit( mode="Gn2P", outroot=outroot, ncat=ncat,  nrep=nrep, seed=seed )
res3=doit( mode="P2G", outroot=outroot, ncat=ncat,  nrep=nrep, seed=seed )
res4=doit( mode="P2Gn", outroot=outroot, ncat=ncat,  nrep=nrep, seed=seed )

statdf=rbind(res1,res2,res3,res4)

outfile=paste(outroot, "stat_all_",  seed, "_", ncat, ".csv", sep="" )
write.table( statdf, outfile, row.names=FALSE, sep=", " )

savefile=paste(outroot, "stat_all_", seed, "_", ncat, ".Rdata", sep="" )
save(statdf,file=savefile)

try( detach(statdf) )
attach(statdf)

boxplot( rmse_p~mode+method, data=statdf, main="rmse_p" )
boxplot( rmse_ii~mode+method, data=statdf, main="rmse_ii" )
boxplot( rmse_t~mode+method, data=statdf, main="rmse_t" )





comments(
 '



 try( detach(statdf) )
 load(file="RpgmORG/stat/stat_all_1701_3.RData")
 statdfall3=statdf
 load(file="RpgmORG/stat/stat_all_1702_3.RData")
 statdfall3=rbind(statdfall3,statdf)
 load(file="RpgmORG/stat/stat_all_1703_3.RData")
 statdfall3=rbind(statdfall3,statdf)
 # str(statdfall3)
 attach(statdfall3)

 boxplot( rmse_p~mode+method, data=statdf, main="rmse_p" )
 boxplot( rmse_ii~mode+method, data=statdf, main="rmse_ii" )
 boxplot( rmse_t~mode+method, data=statdf, main="rmse_t" )

 anova( lm( rmse_p~mode*method) )
 anova( lm( rmse_ii~mode*method) )
 anova( lm( rmse_t~mode*method) )




 try( detach(statdf) )
 load(file="RpgmORG/stat/stat_all_1701_4.RData")
 statdfall4=statdf
 load(file="RpgmORG/stat/stat_all_1702_4.RData")
 statdfall4=rbind(statdfall4,statdf)
 load(file="RpgmORG/stat/stat_all_1703_4.RData")
 statdfall4=rbind(statdfall4,statdf)
 # str(statdfall4)
 attach(statdfall4)

 boxplot( rmse_p~mode+method, data=statdf, main="rmse_p" )
 boxplot( rmse_ii~mode+method, data=statdf, main="rmse_ii" )
 boxplot( rmse_t~mode+method, data=statdf, main="rmse_t" )

 anova( lm( rmse_p~mode*method) )
 anova( lm( rmse_ii~mode*method) )
 anova( lm( rmse_t~mode*method) )




 try( detach(statdf) )
 load(file="RpgmORG/stat/stat_all_1701_5.RData")
 statdfall5=statdf
 attach(statdfall5)







statdf345=stack_df(statdfall3,statdfall4,statdfall5)
statdf345=statdf345[,c("mode", "seed", "ncat", "p11", "p12", "p13", "p14", "p15", "th10", "th11", "th12", "th13", "th14", "a_iinfo1", "method", "p21", "p22", "p23", "p24", "p25", "rmse_p", "rmse_ii", "rmse_iic", "rmse_t", "rmse_tt", "th20", "th21", "th22", "th23", "th24", "a_iinfo2", "difai", "dth0", "dth1", "dth2", "dth3", "dth4")]











 boxplot( rmse_p~mode+method, data=statdf, main="rmse_p" )
 boxplot( rmse_ii~mode+method, data=statdf, main="rmse_ii" )
 boxplot( rmse_t~mode+method, data=statdf, main="rmse_t" )

 anova( lm( rmse_p~mode*method) )
 anova( lm( rmse_ii~mode*method) )
 anova( lm( rmse_t~mode*method) )



 boxplot( rmse_p~mode, data=statdf, main="rmse_p" )
 boxplot( rmse_ii~mode, data=statdf, main="rmse_ii" )
 boxplot( rmse_t~mode, data=statdf, main="rmse_t" )

 boxplot( rmse_p~method, data=statdf, main="rmse_p" )
 boxplot( rmse_ii~method, data=statdf, main="rmse_ii" )
 boxplot( rmse_t~method, data=statdf, main="rmse_t" )





temp=mands(  statdf[,vartype(statdf)=="numeric"], by=statdf[,"method"]  )
temp=mands(  statdf[,vartype(statdf)=="numeric"], by=statdf[,c("mode","method")]  )
get_valL( temp, ,,1)










 # get bad fit params
 statdfs=sortit( statdf, cols=-which(colnames(statdf)=="rmse_ii") )

  param0=data.frame( name="Q1", type="P", stringsAsFactors=0 )
  param1=cbind(param0,statdfs[1:10,c("ncat", "p11", "p12", "p13", "p14")])
  param1$type=substr( unlist(
    lapply( strsplit(statdfs[1:10,"mode"],"2"), "[[", 2 ) ), 1,1 )
  colnames(param1)=gsub("p1","p", colnames(param1))

 npoint=51
 theta=seq(-4,4,length=npoint)
 thd=dnorm(theta)
 thd=thd/sum(thd)

 res=fitGn2P_ls( param1, theta, print=1, plot=1, method=1 )




 '
)

















