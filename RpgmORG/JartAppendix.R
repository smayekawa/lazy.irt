# 必要なパッケージの読み込み
library(lazy.irt)
# library(lazy.tools)
# library(lazy.mat)
# library(lazy.accel)

# 3PLM を 2PLM に変換 （conv2G を用いても良い。）
title="3PLM -> 2PLM"
itemparam=paramA1[2,]; ncat0=itemparam$ncat
nparam0=ncat0; if( itemparam$type %in% c("B3","Bn3") ) nparam0=nparam0+1

# method=1 conversion
res1=conv2P( itemparam, plot=4, print=1, npoints=51, method=1, init=1 )
rmse_icrf1=res1$rmse_icrf; rmse_irf1=res1$rmse_irf; rmse_icif1=res1$rmse_icif; rmse_iif1=res1$rmse_iif
param1=res1$paramNew[,3+(1:nparam0)]

# method=2 conversion
res2=conv2P( itemparam, plot=4, print=1, npoints=51, method=2, init=1 )
rmse_icrf2=res2$rmse_icrf; rmse_irf2=res2$rmse_irf; rmse_icif2=res2$rmse_icif; rmse_iif2=res2$rmse_iif
param2=res2$paramNew[,3+(1:nparam0)]

# accumulate results
rmse_icrf=rbind(rmse_icrf1,rmse_icrf2); rmse_irf=rbind(rmse_irf1,rmse_irf2)
rmse_icif=rbind(rmse_icif1,rmse_icif2); rmse_iif=rbind(rmse_iif1,rmse_iif2)
rownames(rmse_icrf)=1:2; rownames(rmse_irf)=1:2; rownames(rmse_icif)=1:2; rownames(rmse_iif)=1:2
paramall=rbind(itemparam[1,3+(1:nparam0)],param1,param2)
rownames(paramall)=c("3PLM",1:2)

cat("\n",title,"\n")
Print(rmse_icrf,rmse_irf,rmse_icif,rmse_iif, fmt=".4")
Print(paramall)


# GPCM を GRM-L に変換
title="P -> G4"
itemparam=paramA1[8,]; ncat0=itemparam$ncat
nparam0=ncat0; if( itemparam$type %in% c("B3","Bn3") ) nparam0=nparam0+1

# method=1 conversion
res1=conv2G( itemparam, plot=4, print=1, npoints=51, method=1, init=1 )
rmse_icrf1=res1$rmse_icrf; rmse_irf1=res1$rmse_irf; rmse_icif1=res1$rmse_icif; rmse_iif1=res1$rmse_iif
param1=res1$paramNew[,3+(1:ncat0)]

# method=2 conversion
res2=conv2G( itemparam, plot=4, print=1, npoints=51, method=2, init=1 )
rmse_icrf2=res2$rmse_icrf; rmse_irf2=res2$rmse_irf; rmse_icif2=res2$rmse_icif; rmse_iif2=res2$rmse_iif
param2=res2$paramNew[,3+(1:ncat0)]

# accumulate results
rmse_icrf=rbind(rmse_icrf1,rmse_icrf2); rmse_irf=rbind(rmse_irf1,rmse_irf2)
rmse_icif=rbind(rmse_icif1,rmse_icif2); rmse_iif=rbind(rmse_iif1,rmse_iif2)
rownames(rmse_icrf)=1:2; rownames(rmse_irf)=1:2; rownames(rmse_icif)=1:2; rownames(rmse_iif)=1:2
paramall=rbind(itemparam[1,3+(1:nparam0)],param1,param2)
rownames(paramall)=c("GPCM",1:2)

cat("\n",title,"\n")
Print(rmse_icrf,rmse_irf,rmse_icif,rmse_iif, fmt=".4")
Print(paramall)




#  GRM-L を GPCM に変換
title="G -> P4"
itemparam=paramA1[9,]; ncat0=itemparam$ncat
itemparam[,4:7]=c(.9155, -2.0507, 0, 2.0507)
nparam0=ncat0; if( itemparam$type %in% c("B3","Bn3") ) nparam0=nparam0+1

# method=1 conversion
res1=conv2P( itemparam, plot=4, print=1, npoints=51, method=1, init=1 )
rmse_icrf1=res1$rmse_icrf; rmse_irf1=res1$rmse_irf; rmse_icif1=res1$rmse_icif; rmse_iif1=res1$rmse_iif
param1=res1$paramNew[,3+(1:ncat0)]

# method=2 conversion
res2=conv2P( itemparam, plot=4, print=1, npoints=51, method=2, init=1 )
rmse_icrf2=res2$rmse_icrf; rmse_irf2=res2$rmse_irf; rmse_icif2=res2$rmse_icif; rmse_iif2=res2$rmse_iif
param2=res2$paramNew[,3+(1:ncat0)]

# accumulate results
rmse_icrf=rbind(rmse_icrf1,rmse_icrf2); rmse_irf=rbind(rmse_irf1,rmse_irf2)
rmse_icif=rbind(rmse_icif1,rmse_icif2); rmse_iif=rbind(rmse_iif1,rmse_iif2)
rownames(rmse_icrf)=1:2; rownames(rmse_irf)=1:2; rownames(rmse_icif)=1:2; rownames(rmse_iif)=1:2
paramall=rbind(itemparam[1,3+(1:nparam0)],param1,param2)
rownames(paramall)=c("GPCM",1:2)

cat("\n",title,"\n")
Print(rmse_icrf,rmse_irf,rmse_icif,rmse_iif, fmt=".4")
Print(paramall)
