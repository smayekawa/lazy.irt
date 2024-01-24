temp="
 name type ncat p1 p2 p3
Q1 Qg   B3    2    2 0 0.25
Q2 Qh   B3    2    1 0 0
Q3 Qh2  B3    2    1 -0.5 0
";
paramSam=cards(temp,header=1)
paramSam$p1=paramSam$p1/1.7

param=paramSam[c(1,3),]

resInfo=info_func( param, plot=3, print=3
                   , thmin=-4, thmax=1, npoints=151, numderiv=1 )


Uc=matrix(c(1,0),5,2)

theta=resInfo$theta
trf=resInfo$TRF
trf_LO=resInfo$TRF_LO
icrf=resInfo$icrf
dicrf=resInfo$dicrf
fromP=resInfo$fromP
toP=resInfo$toP
v_LO=resInfo$v_LO

matplot(theta,cbind(trf_LO,v_LO[,2]), type="l")


stop()

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
