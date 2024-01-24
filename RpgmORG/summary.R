Px_t=out_obscore2$Px_t
TRFx_t=out_obscore2$theta_stat$TRF
stdx_t=out_obscore2$theta_stat$stdx_t
cols=formatC(-4:4,digits=2,width=5,format="f")
Px_t1=Px_t[,cols]
scorex=colnames(Px_t1)
Px_t1r=Px_t1[nrow(Px_t1):1,]
Print(Px_t1r, fmt="6.2", fuzz=1e-3)
TRFx_t1=TRFx_t[loc]
stdx_t1=stdx_t[loc]
ES=rbind(TRFx_t1,stdx_t1)
rownames(ES)=c("mean","std")
colnames(ES)=colnames(Px_t1r)
Print(ES, fmt="6.2")

title="Plot of TRF and SEM of X"
title2=paste("score range:  X(",0,",",54,")", sep="")
matplot( colnames(Px_t), cbind(TRFx_t,stdx_t), type="l"
         , ylab="X", xlab="theta"
         , lty=1, col=1:2, cex=1, lwd=1
         , main=title, sub=title2 )
legend( 2,max(TRFx_t)-10, legend=c("E(X|theta)","SEM(theta)")
        , cex=1, lwd=1, col=1:2, merge=TRUE )

Py_t=res2$Py_t
TRFy_t=res2$TRFy_t
stdy_t=res2$stdy_t
cols=formatC(-4:4,digits=2,format="f")
loc=which(colnames(Py_t)%in%cols)
Py_t1=Py_t[,loc]
scorey=colnames(Py_t1)
Py_t1r=Py_t1[nrow(Py_t1):1,]
Print(Py_t1r, fmt="6.2", fuzz=1e-3)
TRFy_t1=TRFy_t[loc]
stdy_t1=stdy_t[loc]
ES=rbind(TRFy_t1,stdy_t1)
rownames(ES)=c("mean","std")
colnames(ES)=colnames(Py_t1r)
Print(ES, fmt="6.2")




Py_tc=res2c$Py_t
TRFy_tc=res2c$TRFy_t
stdy_tc=res2c$stdy_t
cols=formatC(-4:4,digits=2,format="f")
loc=which(colnames(Py_t)%in%cols)
Py_tc1=Py_tc[,cols]
scorey=colnames(Py_tc1)
Py_tc1r=Py_tc1[nrow(Py_tc1):1,]
Print(Py_tc1r, fmt="6.2", fuzz=1e-3)
TRFy_tc1=TRFy_tc[loc]
stdy_tc1=stdy_tc[loc]
ESc=rbind(TRFy_tc1,stdy_tc1)
rownames(ESc)=c("mean","std")
colnames(ESc)=colnames(Py_t1r)
Print(ESc, fmt="6.2")


