

title="2PLM -> 2PNM"
title="B2-B2n"
itemparam=paramA1[1,]; ncat0=itemparam$ncat
res0=fitGn2P_ls( itemparam, plot=4, print=1, npoints=51, method=0, init=1 )
rmse_p0=res0$rmse_p; rmse_ii0=res0$rmse_ii; rmse_iic0=res0$rmse_iic
param0=res0$paramNew[,3+(1:ncat0)]

res1=fitGn2P_ls( itemparam, plot=4, print=1, npoints=51, method=1, init=1 )
rmse_p1=res1$rmse_p; rmse_ii1=res1$rmse_ii; rmse_iic1=res1$rmse_iic
param1=res1$paramNew[,3+(1:ncat0)]

res2=fitGn2P_ls( itemparam, plot=4, print=1, npoints=51, method=2, init=1 )
rmse_p2=res2$rmse_p; rmse_ii2=res2$rmse_ii; rmse_iic2=res2$rmse_iic
param2=res2$paramNew[,3+(1:ncat0)]

rmse_p=rbind(rmse_p0,rmse_p1,rmse_p2)
rmse_ii=rbind(rmse_ii0,rmse_ii1,rmse_ii2)
rmse_iic=rbind(rmse_iic0,rmse_iic1,rmse_iic2)
rownames(rmse_p)=0:2; rownames(rmse_ii)=0:2; rownames(rmse_iic)=0:2

paramall=rbind(itemparam[1,3+(1:ncat0)],param0,param1,param2)
rownames(paramall)=c("2PLM",0:2)

cat("\n",title,"\n")
Print(rmse_p,rmse_ii,rmse_iic, fmt=".4")
Print(paramall)

graph2pdf( outpdf=paste("RpgmORG/conv", title, ".pdf", sep="") )
plot_res()
graph2pdf( close=1 )




title="2PNM -> 2PLM"
title="B2n-B2"
itemparam=paramA1[3,]; ncat0=itemparam$ncat
res0=fitP2G_ls( itemparam, plot=4, print=1, npoints=51, method=0, init=1 )
rmse_p0=res0$rmse_p; rmse_ii0=res0$rmse_ii; rmse_iic0=res0$rmse_iic
param0=res0$paramNew[,3+(1:ncat0)]

res1=fitP2G_ls( itemparam, plot=4, print=1, npoints=51, method=1, init=1 )
rmse_p1=res1$rmse_p; rmse_ii1=res1$rmse_ii; rmse_iic1=res1$rmse_iic
param1=res1$paramNew[,3+(1:ncat0)]

res2=fitP2G_ls( itemparam, plot=4, print=1, npoints=51, method=2, init=1 )
rmse_p2=res2$rmse_p; rmse_ii2=res2$rmse_ii; rmse_iic2=res2$rmse_iic
param2=res2$paramNew[,3+(1:ncat0)]

rmse_p=rbind(rmse_p0,rmse_p1,rmse_p2)
rmse_ii=rbind(rmse_ii0,rmse_ii1,rmse_ii2)
rmse_iic=rbind(rmse_iic0,rmse_iic1,rmse_iic2)
rownames(rmse_p)=0:2; rownames(rmse_ii)=0:2; rownames(rmse_iic)=0:2

paramall=rbind(itemparam[1,3+(1:ncat0)], param0 ,param1, param2)
rownames(paramall)=c("2PNM",0:2)

cat("\n",title,"\n")
Print(rmse_p,rmse_ii,rmse_iic, fmt=".4")
Print(paramall)


graph2pdf( outpdf=paste("RpgmORG/conv", title, ".pdf", sep="") )
plot_res()
graph2pdf( close=1 )








#################################################################

title="3PLM -> 2PNM"
title="B3-B2n"
itemparam=paramA1[2,]; ncat0=itemparam$ncat
nparam0=ncat0; if( itemparam$type %in% c("B3","Bn3") ) nparam0=nparam0+1
res0=fitGn2P_ls( itemparam, plot=4, print=1, npoints=51, method=0, init=1 )
rmse_p0=res0$rmse_p; rmse_ii0=res0$rmse_ii; rmse_iic0=res0$rmse_iic
param0=res0$paramNew[,3+(1:nparam0)]

res1=fitGn2P_ls( itemparam, plot=4, print=1, npoints=51, method=1, init=1 )
rmse_p1=res1$rmse_p; rmse_ii1=res1$rmse_ii; rmse_iic1=res1$rmse_iic
param1=res1$paramNew[,3+(1:nparam0)]

res2=fitGn2P_ls( itemparam, plot=4, print=1, npoints=51, method=2, init=1 )
rmse_p2=res2$rmse_p; rmse_ii2=res2$rmse_ii; rmse_iic2=res2$rmse_iic
param2=res2$paramNew[,3+(1:nparam0)]

rmse_p=rbind(rmse_p0,rmse_p1,rmse_p2)
rmse_ii=rbind(rmse_ii0,rmse_ii1,rmse_ii2)
rmse_iic=rbind(rmse_iic0,rmse_iic1,rmse_iic2)
rownames(rmse_p)=0:2; rownames(rmse_ii)=0:2; rownames(rmse_iic)=0:2

paramall=rbind(itemparam[1,3+(1:nparam0)],param0,param1,param2)
rownames(paramall)=c("3PLM",0:2)

cat("\n",title,"\n")
Print(rmse_p,rmse_ii,rmse_iic, fmt=".4")
Print(paramall)


graph2pdf( outpdf=paste("RpgmORG/conv", title, ".pdf", sep="") )
plot_res()
graph2pdf( close=1 )






title="3PLM -> 2PLM"
title="B3-B2"
itemparam=paramA1[2,]; ncat0=itemparam$ncat
nparam0=ncat0; if( itemparam$type %in% c("B3","Bn3") ) nparam0=nparam0+1
res0=fitP2G_ls( itemparam, plot=4, print=1, npoints=51, method=0, init=1 )
rmse_p0=res0$rmse_p; rmse_ii0=res0$rmse_ii; rmse_iic0=res0$rmse_iic
param0=res0$paramNew[,3+(1:nparam0)]

res1=fitP2G_ls( itemparam, plot=4, print=1, npoints=51, method=1, init=1 )
rmse_p1=res1$rmse_p; rmse_ii1=res1$rmse_ii; rmse_iic1=res1$rmse_iic
param1=res1$paramNew[,3+(1:nparam0)]

res2=fitP2G_ls( itemparam, plot=4, print=1, npoints=51, method=2, init=1 )
rmse_p2=res2$rmse_p; rmse_ii2=res2$rmse_ii; rmse_iic2=res2$rmse_iic
param2=res2$paramNew[,3+(1:nparam0)]

rmse_p=rbind(rmse_p0,rmse_p1,rmse_p2)
rmse_ii=rbind(rmse_ii0,rmse_ii1,rmse_ii2)
rmse_iic=rbind(rmse_iic0,rmse_iic1,rmse_iic2)
rownames(rmse_p)=0:2; rownames(rmse_ii)=0:2; rownames(rmse_iic)=0:2

paramall=rbind(itemparam[1,3+(1:nparam0)], param0 ,param1, param2)
rownames(paramall)=c("3PLM",0:2)


cat("\n",title,"\n")
Print(rmse_p,rmse_ii,rmse_iic, fmt=".4")
Print(paramall)


graph2pdf( outpdf=paste("RpgmORG/conv", title, ".pdf", sep="") )
plot_res()
graph2pdf( close=1 )








#################################################################

title="GPCM -> GRM-N   ncat=3"
title="P-Gn3"
itemparam=paramA1[5,]; ncat0=itemparam$ncat
res0=fitGn2P_ls( itemparam, plot=4, print=1, npoints=51, method=0, init=1 )
rmse_p0=res0$rmse_p; rmse_ii0=res0$rmse_ii; rmse_iic0=res0$rmse_iic
param0=res0$paramNew[,3+(1:ncat0)]

res1=fitGn2P_ls( itemparam, plot=4, print=1, npoints=51, method=1, init=1 )
rmse_p1=res1$rmse_p; rmse_ii1=res1$rmse_ii; rmse_iic1=res1$rmse_iic
param1=res1$paramNew[,3+(1:ncat0)]


res2=fitGn2P_ls( itemparam, plot=4, print=1, npoints=51, method=2, init=1 )
rmse_p2=res2$rmse_p; rmse_ii2=res2$rmse_ii; rmse_iic2=res2$rmse_iic
param2=res2$paramNew[,3+(1:ncat0)]

rmse_p=rbind(rmse_p0,rmse_p1,rmse_p2)
rmse_ii=rbind(rmse_ii0,rmse_ii1,rmse_ii2)
rmse_iic=rbind(rmse_iic0,rmse_iic1,rmse_iic2)
rownames(rmse_p)=0:2; rownames(rmse_ii)=0:2; rownames(rmse_iic)=0:2

paramall=rbind(itemparam[1,3+(1:ncat0)],param0,param1,param2)
rownames(paramall)=c("GPCM",0:2)


cat("\n",title,"\n")
Print(rmse_p,rmse_ii,rmse_iic, fmt=".4")
Print(paramall)


graph2pdf( outpdf=paste("RpgmORG/conv", title, ".pdf", sep="") )
plot_res()
graph2pdf( close=1 )






title="GPCM -> GRM-L   ncat=3"
title="P-G3"
itemparam=paramA1[5,]; ncat0=itemparam$ncat
res0=fitG2P_ls( itemparam, plot=4, print=1, npoints=51, method=0, init=1 )
rmse_p0=res0$rmse_p; rmse_ii0=res0$rmse_ii; rmse_iic0=res0$rmse_iic
param0=res0$paramNew[,3+(1:ncat0)]

res1=fitG2P_ls( itemparam, plot=4, print=1, npoints=51, method=1, init=1 )
rmse_p1=res1$rmse_p; rmse_ii1=res1$rmse_ii; rmse_iic1=res1$rmse_iic
param1=res1$paramNew[,3+(1:ncat0)]


res2=fitG2P_ls( itemparam, plot=4, print=1, npoints=51, method=2, init=1 )
rmse_p2=res2$rmse_p; rmse_ii2=res2$rmse_ii; rmse_iic2=res2$rmse_iic
param2=res2$paramNew[,3+(1:ncat0)]

rmse_p=rbind(rmse_p0,rmse_p1,rmse_p2)
rmse_ii=rbind(rmse_ii0,rmse_ii1,rmse_ii2)
rmse_iic=rbind(rmse_iic0,rmse_iic1,rmse_iic2)
rownames(rmse_p)=0:2; rownames(rmse_ii)=0:2; rownames(rmse_iic)=0:2

paramall=rbind(itemparam[1,3+(1:ncat0)],param0,param1,param2)
rownames(paramall)=c("GPCM",0:2)


cat("\n",title,"\n")
Print(rmse_p,rmse_ii,rmse_iic, fmt=".4")
Print(paramall)


graph2pdf( outpdf=paste("RpgmORG/conv", title, ".pdf", sep="") )
plot_res()
graph2pdf( close=1 )






#################################################################

title="GPCM -> GRM-N   ncat=4"
title="P-Gn4"
itemparam=paramA1[8,]; ncat0=itemparam$ncat
res0=fitGn2P_ls( itemparam, plot=4, print=1, npoints=51, method=0, init=1 )
rmse_p0=res0$rmse_p; rmse_ii0=res0$rmse_ii; rmse_iic0=res0$rmse_iic
param0=res0$paramNew[,3+(1:ncat0)]

res1=fitGn2P_ls( itemparam, plot=4, print=1, npoints=51, method=1, init=1 )
rmse_p1=res1$rmse_p; rmse_ii1=res1$rmse_ii; rmse_iic1=res1$rmse_iic
param1=res1$paramNew[,3+(1:ncat0)]


res2=fitGn2P_ls( itemparam, plot=4, print=1, npoints=51, method=2, init=1 )
rmse_p2=res2$rmse_p; rmse_ii2=res2$rmse_ii; rmse_iic2=res2$rmse_iic
param2=res2$paramNew[,3+(1:ncat0)]

rmse_p=rbind(rmse_p0,rmse_p1,rmse_p2)
rmse_ii=rbind(rmse_ii0,rmse_ii1,rmse_ii2)
rmse_iic=rbind(rmse_iic0,rmse_iic1,rmse_iic2)
rownames(rmse_p)=0:2; rownames(rmse_ii)=0:2; rownames(rmse_iic)=0:2

paramall=rbind(itemparam[1,3+(1:ncat0)],param0,param1,param2)
rownames(paramall)=c("GPCM",0:2)


cat("\n",title,"\n")
Print(rmse_p,rmse_ii,rmse_iic, fmt=".4")
Print(paramall)


graph2pdf( outpdf=paste("RpgmORG/conv", title, ".pdf", sep="") )
plot_res()
graph2pdf( close=1 )






title="GPCM -> GRM-L   ncat=4"
title="P-G4"
itemparam=paramA1[8,]; ncat0=itemparam$ncat
res0=fitG2P_ls( itemparam, plot=4, print=1, npoints=51, method=0, init=1 )
rmse_p0=res0$rmse_p; rmse_ii0=res0$rmse_ii; rmse_iic0=res0$rmse_iic
param0=res0$paramNew[,3+(1:ncat0)]

res1=fitG2P_ls( itemparam, plot=4, print=1, npoints=51, method=1, init=1 )
rmse_p1=res1$rmse_p; rmse_ii1=res1$rmse_ii; rmse_iic1=res1$rmse_iic
param1=res1$paramNew[,3+(1:ncat0)]


res2=fitG2P_ls( itemparam, plot=4, print=1, npoints=51, method=2, init=1 )
rmse_p2=res2$rmse_p; rmse_ii2=res2$rmse_ii; rmse_iic2=res2$rmse_iic
param2=res2$paramNew[,3+(1:ncat0)]

rmse_p=rbind(rmse_p0,rmse_p1,rmse_p2)
rmse_ii=rbind(rmse_ii0,rmse_ii1,rmse_ii2)
rmse_iic=rbind(rmse_iic0,rmse_iic1,rmse_iic2)
rownames(rmse_p)=0:2; rownames(rmse_ii)=0:2; rownames(rmse_iic)=0:2

paramall=rbind(itemparam[1,3+(1:ncat0)],param0,param1,param2)
rownames(paramall)=c("GPCM",0:2)


cat("\n",title,"\n")
Print(rmse_p,rmse_ii,rmse_iic, fmt=".4")
Print(paramall)



graph2pdf( outpdf=paste("RpgmORG/conv", title, ".pdf", sep="") )
plot_res()
graph2pdf( close=1 )







#################################################################

title="GRM-N -> GPCM   ncat=3"
title="Gn-P3"
itemparam=paramA1[7,]; ncat0=itemparam$ncat
res0=fitP2G_ls( itemparam, plot=4, print=1, npoints=51, method=0, init=1 )
rmse_p0=res0$rmse_p; rmse_ii0=res0$rmse_ii; rmse_iic0=res0$rmse_iic
param0=res0$paramNew[,3+(1:ncat0)]

res1=fitP2G_ls( itemparam, plot=4, print=1, npoints=51, method=1, init=1 )
rmse_p1=res1$rmse_p; rmse_ii1=res1$rmse_ii; rmse_iic1=res1$rmse_iic
param1=res1$paramNew[,3+(1:ncat0)]


res2=fitP2G_ls( itemparam, plot=4, print=1, npoints=51, method=2, init=1 )
rmse_p2=res2$rmse_p; rmse_ii2=res2$rmse_ii; rmse_iic2=res2$rmse_iic
param2=res2$paramNew[,3+(1:ncat0)]

rmse_p=rbind(rmse_p0,rmse_p1,rmse_p2)
rmse_ii=rbind(rmse_ii0,rmse_ii1,rmse_ii2)
rmse_iic=rbind(rmse_iic0,rmse_iic1,rmse_iic2)
rownames(rmse_p)=0:2; rownames(rmse_ii)=0:2; rownames(rmse_iic)=0:2

paramall=rbind(itemparam[1,3+(1:ncat0)],param0,param1,param2)
rownames(paramall)=c("GRM-N",0:2)


cat("\n",title,"\n")
Print(rmse_p,rmse_ii,rmse_iic, fmt=".4")
Print(paramall)



graph2pdf( outpdf=paste("RpgmORG/conv", title, ".pdf", sep="") )
plot_res()
graph2pdf( close=1 )





title="GRM-L -> GPCM   ncat=3"
title="G-P3"
itemparam=paramA1[6,]; ncat0=itemparam$ncat
res0=fitP2G_ls( itemparam, plot=4, print=1, npoints=51, method=0, init=1 )
rmse_p0=res0$rmse_p; rmse_ii0=res0$rmse_ii; rmse_iic0=res0$rmse_iic
param0=res0$paramNew[,3+(1:ncat0)]

res1=fitP2G_ls( itemparam, plot=4, print=1, npoints=51, method=1, init=1 )
rmse_p1=res1$rmse_p; rmse_ii1=res1$rmse_ii; rmse_iic1=res1$rmse_iic
param1=res1$paramNew[,3+(1:ncat0)]


res2=fitP2G_ls( itemparam, plot=4, print=1, npoints=51, method=2, init=1 )
rmse_p2=res2$rmse_p; rmse_ii2=res2$rmse_ii; rmse_iic2=res2$rmse_iic
param2=res2$paramNew[,3+(1:ncat0)]

rmse_p=rbind(rmse_p0,rmse_p1,rmse_p2)
rmse_ii=rbind(rmse_ii0,rmse_ii1,rmse_ii2)
rmse_iic=rbind(rmse_iic0,rmse_iic1,rmse_iic2)
rownames(rmse_p)=0:2; rownames(rmse_ii)=0:2; rownames(rmse_iic)=0:2

paramall=rbind(itemparam[1,3+(1:ncat0)],param0,param1,param2)
rownames(paramall)=c("GRM-N",0:2)


cat("\n",title,"\n")
Print(rmse_p,rmse_ii,rmse_iic, fmt=".4")
Print(paramall)



graph2pdf( outpdf=paste("RpgmORG/conv", title, ".pdf", sep="") )
plot_res()
graph2pdf( close=1 )




title="GRM-L -> GRM-N   ncat=3"
title="G-Gn3"
itemparam=paramA1[6,]; ncat0=itemparam$ncat
res0=fitGn2P_ls( itemparam, plot=4, print=1, npoints=51, method=0, init=1 )
rmse_p0=res0$rmse_p; rmse_ii0=res0$rmse_ii; rmse_iic0=res0$rmse_iic
param0=res0$paramNew[,3+(1:ncat0)]

res1=fitGn2P_ls( itemparam, plot=4, print=1, npoints=51, method=1, init=1 )
rmse_p1=res1$rmse_p; rmse_ii1=res1$rmse_ii; rmse_iic1=res1$rmse_iic
param1=res1$paramNew[,3+(1:ncat0)]


res2=fitGn2P_ls( itemparam, plot=4, print=1, npoints=51, method=2, init=1 )
rmse_p2=res2$rmse_p; rmse_ii2=res2$rmse_ii; rmse_iic2=res2$rmse_iic
param2=res2$paramNew[,3+(1:ncat0)]

rmse_p=rbind(rmse_p0,rmse_p1,rmse_p2)
rmse_ii=rbind(rmse_ii0,rmse_ii1,rmse_ii2)
rmse_iic=rbind(rmse_iic0,rmse_iic1,rmse_iic2)
rownames(rmse_p)=0:2; rownames(rmse_ii)=0:2; rownames(rmse_iic)=0:2

paramall=rbind(itemparam[1,3+(1:ncat0)],param0,param1,param2)
rownames(paramall)=c("GRM-N",0:2)


cat("\n",title,"\n")
Print(rmse_p,rmse_ii,rmse_iic, fmt=".4")
Print(paramall)



graph2pdf( outpdf=paste("RpgmORG/conv", title, ".pdf", sep="") )
plot_res()
graph2pdf( close=1 )
















#################################################################

title="GRM-N -> GPCM   ncat=4"
title="Gn-P4"
itemparam=paramA1[10,]; ncat0=itemparam$ncat
res0=fitP2G_ls( itemparam, plot=4, print=1, npoints=51, method=0, init=1 )
rmse_p0=res0$rmse_p; rmse_ii0=res0$rmse_ii; rmse_iic0=res0$rmse_iic
param0=res0$paramNew[,3+(1:ncat0)]

res1=fitP2G_ls( itemparam, plot=4, print=1, npoints=51, method=1, init=1 )
rmse_p1=res1$rmse_p; rmse_ii1=res1$rmse_ii; rmse_iic1=res1$rmse_iic
param1=res1$paramNew[,3+(1:ncat0)]


res2=fitP2G_ls( itemparam, plot=4, print=1, npoints=51, method=2, init=1 )
rmse_p2=res2$rmse_p; rmse_ii2=res2$rmse_ii; rmse_iic2=res2$rmse_iic
param2=res2$paramNew[,3+(1:ncat0)]

rmse_p=rbind(rmse_p0,rmse_p1,rmse_p2)
rmse_ii=rbind(rmse_ii0,rmse_ii1,rmse_ii2)
rmse_iic=rbind(rmse_iic0,rmse_iic1,rmse_iic2)
rownames(rmse_p)=0:2; rownames(rmse_ii)=0:2; rownames(rmse_iic)=0:2

paramall=rbind(itemparam[1,3+(1:ncat0)],param0,param1,param2)
rownames(paramall)=c("GRM-N",0:2)


cat("\n",title,"\n")
Print(rmse_p,rmse_ii,rmse_iic, fmt=".4")
Print(paramall)


graph2pdf( outpdf=paste("RpgmORG/conv", title, ".pdf", sep="") )
plot_res()
graph2pdf( close=1 )






title="GRM-L -> GPCM   ncat=4"
title="G-P4"
itemparam=paramA1[9,]; ncat0=itemparam$ncat
res0=fitP2G_ls( itemparam, plot=4, print=1, npoints=51, method=0, init=1 )
rmse_p0=res0$rmse_p; rmse_ii0=res0$rmse_ii; rmse_iic0=res0$rmse_iic
param0=res0$paramNew[,3+(1:ncat0)]

res1=fitP2G_ls( itemparam, plot=4, print=1, npoints=51, method=1, init=1 )
rmse_p1=res1$rmse_p; rmse_ii1=res1$rmse_ii; rmse_iic1=res1$rmse_iic
param1=res1$paramNew[,3+(1:ncat0)]


res2=fitP2G_ls( itemparam, plot=4, print=1, npoints=51, method=2, init=1 )
rmse_p2=res2$rmse_p; rmse_ii2=res2$rmse_ii; rmse_iic2=res2$rmse_iic
param2=res2$paramNew[,3+(1:ncat0)]

rmse_p=rbind(rmse_p0,rmse_p1,rmse_p2)
rmse_ii=rbind(rmse_ii0,rmse_ii1,rmse_ii2)
rmse_iic=rbind(rmse_iic0,rmse_iic1,rmse_iic2)
rownames(rmse_p)=0:2; rownames(rmse_ii)=0:2; rownames(rmse_iic)=0:2

paramall=rbind(itemparam[1,3+(1:ncat0)],param0,param1,param2)
rownames(paramall)=c("GRM-L",0:2)


cat("\n",title,"\n")
Print(rmse_p,rmse_ii,rmse_iic, fmt=".4")
Print(paramall)



graph2pdf( outpdf=paste("RpgmORG/conv", title, ".pdf", sep="") )
plot_res()
graph2pdf( close=1 )




title="GRM-L -> GRM-N   ncat=4"
title="G-Gn4"
itemparam=paramA1[9,]; ncat0=itemparam$ncat
res0=fitGn2P_ls( itemparam, plot=4, print=1, npoints=51, method=0, init=1 )
rmse_p0=res0$rmse_p; rmse_ii0=res0$rmse_ii; rmse_iic0=res0$rmse_iic
param0=res0$paramNew[,3+(1:ncat0)]

res1=fitGn2P_ls( itemparam, plot=4, print=1, npoints=51, method=1, init=1 )
rmse_p1=res1$rmse_p; rmse_ii1=res1$rmse_ii; rmse_iic1=res1$rmse_iic
param1=res1$paramNew[,3+(1:ncat0)]


res2=fitGn2P_ls( itemparam, plot=4, print=1, npoints=51, method=2, init=1 )
rmse_p2=res2$rmse_p; rmse_ii2=res2$rmse_ii; rmse_iic2=res2$rmse_iic
param2=res2$paramNew[,3+(1:ncat0)]

rmse_p=rbind(rmse_p0,rmse_p1,rmse_p2)
rmse_ii=rbind(rmse_ii0,rmse_ii1,rmse_ii2)
rmse_iic=rbind(rmse_iic0,rmse_iic1,rmse_iic2)
rownames(rmse_p)=0:2; rownames(rmse_ii)=0:2; rownames(rmse_iic)=0:2

paramall=rbind(itemparam[1,3+(1:ncat0)],param0,param1,param2)
rownames(paramall)=c("GRM-L",0:2)


cat("\n",title,"\n")
Print(rmse_p,rmse_ii,rmse_iic, fmt=".4")
Print(paramall)



graph2pdf( outpdf=paste("RpgmORG/conv", title, ".pdf", sep="") )
plot_res()
graph2pdf( close=1 )





################ back to the original ##########################






