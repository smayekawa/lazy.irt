
param=paramB1[c(1:3,7:9,13:15),]

param=paramS2[paramS2$type=="G",]

param=paramS1


#param[1,5:7]=param[1,5:7]+1
pres=fitP2G_ls( param[3,], plot=1 )
param[4,]=pres$paramP
param[4,"name"]="Q31p"
#param[2,"type"]="P"
#param[2,"p1"]=0.802


#param=paramB1[2:3,]
#param[2,"type"]="G"
#param[2,"p2"]=0

resInfo=info_func( param, plot=3, print=3
                   , thmin=-3, thmax=3, npoints=151 )




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
