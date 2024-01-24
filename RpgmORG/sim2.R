seed=1701

nitem=20
param=data.frame(name=paste("Q",1:nitem,sep=""), type="B", ncat=2
                 , p1=1, p2=0, p3=0, stringsAsFactors=0 )

nsim=2

# prepare storege for the fit stats

for( rep in 1:nsim ){

 # generate parameter data frame
 a=0.5+runif(nitem)
 b=rnorm(nitem, 0, 0.8)
 param[,c("p1","p2")]=cbind(a,b)

 Print(rep)
 Print(param)

 ################################################
 #  main body
 ################################################

 # accumulate fit statistics such as rmse_p maxad_t etc
 #

} # end of rep







