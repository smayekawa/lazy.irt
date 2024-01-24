
param=paramS1
param[c(1,4),"p1"]=1

locP=which(param$type=="P")
locG=which(param$type=="G")

thmin=-4
thmax=4
npoints=121
theta=seq(thmin,thmax,length=npoints)


#param[1,5:7]=param[1,5:7]+1
pres=fitG2P_ls( param[locP,], plot=1, thmin=thmin, thmax=thmax
                , npoints=npoints )
param[locG,]=pres$paramG
param[locG,"name"]=paste(param$name[locP],"g", sep="")

# param[1:2,"p1"]=param[4,"p1"]

# param=param[-3,]


param=paramS2[5:6,]
param$type="Gn"

resInfo=info_func( param, plot=4, print=3, theta=theta, numderiv=0, smallP=0 )




theta=resInfo$theta
trf=resInfo$TRF
trf_LO=resInfo$TRF_LO
icrf=resInfo$icrf
dicrf=resInfo$dicrf
fromP=resInfo$fromP
toP=resInfo$toP
v_LO=resInfo$v_LO

# tcc with best v from P'0/P0
tt=-dicrf[,fromP]/icrf[,fromP]
ttt=rowSums(tt)
maxt_1.7=apply(tt,2,max)/1.7
logP0=log(icrf[,fromP])

# tcc with best v
vp1=v_LO*icrf
vp=mapply( function(x,y){rowSums(vp1[,x:y])}, fromP,toP )
Print(max(abs(tt-vp)))

Print(theta,trf,trf_LO,ttt)
Print(maxt_1.7,fmt="6.3")

matplot(theta,cbind(trf,trf_LO,ttt), type="l", main="TRF and TRF_LO")

matplot(theta,tt, type="l", main="Best Item Score")

matplot(theta,logP0, type="l", main="log P0(theta)")



vmax=v_LO[length(theta),toP]
Print(vmax,max(trf_LO))
locG=which(param$type=="G")
tccmax=1.7*sum((param$ncat-1)[-locG]*param$p1[-locG]) +
 1.7*sum(param$p1[locG])
Print(tccmax,max(trf_LO))




comment(
 '
# v=1.7 k a
vp=v_LO[npoints,fromP[4]:toP[4]]
Print(vp/1.7/param[4,"p1"])

vp1=vp/vp[4]*3


paramx=param[3,]

res=icrfG( paramx, theta, print=0 )
str(res)
PP=res$PP
PQ=PP*(1-PP)

k=1
vk=(PQ[,k +1]-PQ[,k+1 +1])/(PP[,k +1]-PP[,k+1 +1]) + (PP[,1 +1])
plot(theta,vk,type="l", ylim=c(0,2), main=paste("k =",k,"  vk"))
k=2
vk=(PQ[,k +1]-PQ[,k+1 +1])/(PP[,k +1]-PP[,k+1 +1]) + (PP[,1 +1])
plot(theta,vk,type="l", ylim=c(0,2), main=paste("k =",k,"  vk"))
k=3
vk=(PQ[,k +1]-PQ[,k+1 +1])/(PP[,k +1]-PP[,k+1 +1]) + (PP[,1 +1])
plot(theta,vk,type="l", ylim=c(0,2), main=paste("k =",k,"  vk"))


k=1
vk=(PQ[,k +1]-PQ[,k+1 +1])/(PP[,k +1]-PP[,k+1 +1])
plot(theta,vk,type="l", ylim=c(-1,1), main=paste("k =",k,"  vk+PP1 "))
k=2
vk=(PQ[,k +1]-PQ[,k+1 +1])/(PP[,k +1]-PP[,k+1 +1])
plot(theta,vk,type="l", ylim=c(-1,1), main=paste("k =",k,"  vk+PP1 "))
k=3
vk=(PQ[,k +1]-PQ[,k+1 +1])/(PP[,k +1]-PP[,k+1 +1])
plot(theta,vk,type="l", ylim=c(-1,1), main=paste("k =",k,"  vk+PP1 "))


')



