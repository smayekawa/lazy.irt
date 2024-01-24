

t=seq(-2,2,0.5)
u=10*logistic(t)


Print(t,u)
plot(t,u, type="b")


u05=seq(0,10,.5)
u05=seq(1,10)-0.5
tb=interpol( u,t, u05 )[,2]
tb=round(tb*100)/100

Print(u05,tb, fmt="8.4")

cc=cut( t, tb, include.lowest=TRUE)
levelname=attributes(cc)$levels
cc=as.numeric(cc)
ccn=levelname[cc]
ccu=u05[cc+1]-0.5

Print(t,cc,ccn,ccu)
