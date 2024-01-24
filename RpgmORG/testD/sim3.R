
library(lazy.irt)

seed=1701

nitem=20
param=data.frame(name=paste("Q",1:nitem,sep=""), type="B", ncat=2
                 , p1=1, p2=0, p3=0, stringsAsFactors=0 )

nsim=10
x <- matrix(ncol=2,nrow=nsim)
i <- 1

# prepare storege for the fit stats

for( rep in 1:nsim ){

 # generate parameter data frame
 a=0.5+runif(nitem,0,1.5)
 b=rnorm(nitem, 0, 0.8)
 param[,c("p1","p2")]=cbind(a,b)

 #  Print(rep)
 #  Print(param)

 ################################################
 nitem=nrow(param)

 thmin <- -2; thmax <- 2; npoint <- 5
 N <- 1000

 # discrete theta
 # theta0 <- seq(thmin,thmax,length=npoint)
 theta0 <- c(-2, -1, 0, 2, 3)
 theta <- unlist(lapply( theta0, rep, round(N/npoint) ))
 res2 <- gendataIRT( 1, param, theta=theta, compress=1 )
 Uc <- as.data.frame(res2$U)
 ncat <- res2$ncat
 type <- res2$type

 # lrt parameters
 nclass <- 5
 resm1 <- uLRT( Uc, nclass=nclass, estrho=1, monotone=1, alpha=20
                , maxiter=200, plot=1, print=1 )
 V1 <- resm1$V[,seq(2,2*resm1$nitems,2)]

 # conversion
 res <- fitI2L_ls( t(V1), print=1, plot=1 )
 param2=res$param
 theta2=res$theta


 ##################### equate result of fitI2L to the original

 ind=rbind(param,param2)
 ind[1:nitem,"gfid"]=1
 ind[(nitem+1):(2*nitem),"gfid"]=2
 ind[1:nitem,"giid"]=param$name
 ind[(nitem+1):(2*nitem),"giid"]=param$name

 resc=calr( ind, baseform=1 )

 #paramc=resc$param
 qr=resc$qr

 paramc2=param2
 paramc2$p1=paramc2$p1/qr[2,2]
 paramc2$p2=param2$p2*qr[2,2]+qr[2,1]
 theta2=theta2*qr[2,2]+qr[2,1]

 #  plot(theta0, theta2,type="b", main="original theta vs recovered theta")

 para1=param[,4:6]
 para2=paramc2[,4:6]
 dif_p=para1-para2
 maxad_p=max(abs(dif_p))
 rmse_p=sqrt(sum(dif_p^2)/nitem/2)
 dif_t=theta0-theta2
 maxad_t=max(abs(dif_t))
 rmse_t=sqrt(sum(dif_t^2)/nitem/2)
 #  Print(rmse_p, maxad_p, rmse_t, maxad_t)
 y <- c(rmse_p, rmse_t)
 x[i, ] <- y
 i <- i+1
 ################################################

 # accumulate fit statistics such as rmse_p maxad_t etc
 #

} # end of rep
