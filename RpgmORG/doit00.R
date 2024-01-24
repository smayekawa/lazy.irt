# 2PLM -> 2PNM -> 2PLM
itemparam=paramA1[1,]; ncat0=itemparam$ncat
res0=fitGn2P_ls( itemparam, plot=0, print=0, wtype=1, method=0, init=0 )
res00=fitP2G_ls( res0$paramG, plot=4, print=1, wtype=1, method=0, init=0 )
param_p0=res00$paramP[,3+(1:ncat0)]

res1=fitGn2P_ls( itemparam, plot=0, print=0, wtype=1, method=1, init=0 )
res11=fitP2G_ls( res1$paramG, plot=4, print=1, wtype=1, method=1, init=0 )
param_p1=res11$paramP[,3+(1:ncat0)]

res2=fitGn2P_ls( itemparam, plot=0, print=0, wtype=1, method=2, init=0 )
res22=fitP2G_ls( res2$paramG, plot=4, print=1, wtype=1, method=2, init=0 )
param_p2=res22$paramP[,3+(1:ncat0)]

paramall=rbind(itemparam[1,3+(1:ncat0)],param_p0,param_p1,param_p2)
rownames(paramall)=c("2PLM",0:2)
Print(paramall)



