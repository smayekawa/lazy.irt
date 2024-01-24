
nitem=4
prob=rep(0.5,nitem)
pmat=rbind(prob,1-prob)
w1=rep(1,nitem)
w2=c(1,1,2,2)
w3=c(2,2,3,3)
w4=c(1,1,4,4)
w5=c(1,2,3,4)
sumsmnw( P=pmat, w=w1, print=9,plot=1 )
# dbinom(0:nitem,nitem,0.5)
#
sumsmnw( P=pmat, w=w2, print=9,plot=1 )
sumsmnw( P=pmat, w=w3, print=9,plot=1 )
sumsmnw( P=pmat, w=w4, print=9,plot=1 )


# simulate sumsmnw

U=gen01pat( nitem, sort=1 )
p=exp( U%*%log(prob) + (1-U)%*%log(1-prob) )

s1=U%*%w1
f1=table2freqdist(xtabs(p~s1))
tt=paste("prob=", paste(prob,collapse=", "), ",  w=", paste(w1,collapse=", "))
plot_freqdist(f1, title=tt)

s2=U%*%w2
f2=table2freqdist(xtabs(p~s2))
tt=paste("prob=", paste(prob,collapse=", "), ",  w=", paste(w2,collapse=", "))
plot_freqdist(f2, title=tt)

s3=U%*%w3
f3=table2freqdist(xtabs(p~s3))
tt=paste("prob=", paste(prob,collapse=", "), ",  w=", paste(w3,collapse=", "))
plot_freqdist(f3, title=tt)

s4=U%*%w4
f4=table2freqdist(xtabs(p~s4))
tt=paste("prob=", paste(prob,collapse=", "), ",  w=", paste(w4,collapse=", "))
plot_freqdist(f4, title=tt)

s5=U%*%w5
f5=table2freqdist(xtabs(p~s5))
tt=paste("prob=", paste(prob,collapse=", "), ",  w=", paste(w5,collapse=", "))
plot_freqdist(f5, title=tt)

score=cbind(s1,s2,s3,s4,s5)
Print(U,p,score)


ff=data.frame(score=(0:sum(w3)))
merge( ff, f3, by.x="score", by.y="s3", all=TRUE)
