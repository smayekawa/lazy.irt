# install packages
# install.packages("D:/RPGM/packages/lazy.irt_0.1.3.zip" , repos = NULL)
# install.packages("D:/RPGM/packages/lazy.tools_0.1.3.zip" , repos = NULL)
# install.packages("D:/RPGM/packages/lazy.mat_0.1.3.zip" , repos = NULL)
# install.packages("D:/RPGM/packages/lazy.accel_0.1.3.zip" , repos = NULL)

# use lazy.irt package
library(lazy.irt)

# discrete theta points
thmin=-4
thmax=4
npoints=121
theta=seq(thmin,thmax,length=npoints)

# parameter data frame
param=paramA1

# plot of icrf at theta
resirf=irf(param, theta, plot=1, print=0)

# calculate and plot the information functions and the best weights
resInfo=info_func( param, plot=4, print=0, theta=theta, numderiv=0, smallP=0 )

# store the result
theta=resInfo$theta
trf=resInfo$TRF
trf_LO=resInfo$TRF_LO
icrf=resInfo$icrf
dicrf=resInfo$dicrf
fromP=resInfo$fromP
toP=resInfo$toP
v_LO=resInfo$v_LO
w_LO=resInfo$w_LO
iname=param$name

# aveage best weights
abv=dnorm(theta)%*%v_LO/sum(dnorm(theta))
cn=colnames(abv)
abv=c(abv)
round(abv, digits=3)
Print(cn, abv,diff(abv))

# difference between the info with LISF and the info with the best weights
Print( max(abs(resInfo$info_LOW-resInfo$info_LO)) )


# irf with best v  = - P'0/P0
irfv=-dicrf[,fromP]/icrf[,fromP]
trfv=rowSums(irfv)
maxirfv=apply(irfv,2,max)
logP0=log(icrf[,fromP])

# irf with best v = v*icrf
vp1=v_LO*icrf
vp=mapply( function(x,y){rowSums(vp1[,x:y])}, fromP,toP )

# above two should be the same.
Print(max(abs(irfv-vp)))

# max value of irf at theta=4
Print(maxirfv)

# comparizon of trf: should be the same.
Print( max(abs(trfv-trf_LO)) )
# Print(theta,trf,trf_LO,trfv)
Print(maxirfv,fmt="6.3")

# plot of irf with the best weights
matplot( theta, irfv, type="l"
         , col=c("black", "red", "green", "blue"), lty=c(1:3), lwd=2
         , xlim=c(min(theta),max(theta)), ylim=c(0,max(irfv))
         , main="Item Response Function with Best Weight", ylab="score")
legend( range(theta)[1], max(irfv)-.1, iname
        , cex=0.8, col=c("black", "red", "green", "blue")
        , pch=NA, lty=c(1:3), title="legend" )

# max value of best item weight at theta=4
Print(apply(w_LO,2,max))


comment(
 '













param1=paramA1[c(5,8),]
param1$p2=-3

param1[2,c("p2","p3")]=c(-1,1)
temp=irf(param1,plot=1)


param2=fitGn2P_ls( param1, theta, maxiter=20, plot=1, wtype=1, method=0 )$paramG

param=rbind(param1,param2)


matplot( theta,irfv[,5:10], type="l"
         , col=c("black", "red", "green", "blue"), lty=c(1:3), lwd=2
         , xlim=c(min(theta),max(theta)), ylim=c(0,max(irfv))
         , main="Item Response Function with Best Weight", ylab="score")
 legend( range(theta)[1], max(irfv)-.1, iname[5:10]
         , cex=0.8, col=c("black", "red", "green", "blue")
         , pch=NA, lty=c(1:3), title="legend" )













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



