set.seed(1701)

n=100
p1=0.6+runif(n)
p2=rnorm(n)
p3=0.5*runif(n)

paramB3=data.frame( name=paste("Q",1:n,sep=""), type="B3", ncat=2,
                    p1=p1, p2=p2, p3=p3, stringsAsFactors=0 )

res2 <- fit223_ls( paramB3, plot=1, print=0, wtype=1 )

mean(res2$rmse)
