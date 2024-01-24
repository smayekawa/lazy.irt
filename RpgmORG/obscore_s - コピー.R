################# This is the core of this program #######################
# calculate cond dist of x given theta as the sum of scored multinomials #
##########################################################################

npoints=31; thmin=-4; thmax=4
theta=seq(thmin,thmax,length.out=npoints)
Pt=exp(-0.5*theta^2);
Pt=Pt/sum(Pt)
param1=paramS1; weight1=weightS11
ncat=param1$ncat
nitems=nrow(param1)
w=weight1$w
v=weight1[,5:8]

obscore_s <- function( param, weight
                       , thmin=-4, thmax=4, npoints=31, thdist=1 ){
 # Simplest verssion of lazy.irt::obscore.
 # Shin-ichi Mayekawa
 # 20230717
 #

 # generate thata and prior theta dist
 theta=seq(thmin,thmax,length.out=npoints)
 thname=format(theta,digits=2)
 if( thdist == 0 ) Pt=matrix(1/npoints,npoints)
 else if( thdist == 1 ){
  Pt=exp(-0.5*theta^2);
  Pt=Pt/sum(Pt)
 }

 # calculate icrf and trf=cond. mean of X given theta
 temp=irf( param, theta, weight, print=0, debug=0, plot=0 )
 icrf=temp$ICRF
 trf=temp$TRF
 fromP=temp$fromP
 toP=temp$toP
 maxscore_t=temp$maxscore_t
 rm(temp)

 # conditional dist of x given theta
 Px_t=matrix(0,maxscore_t+1,npoints)
 Px_t=NULL
 for( k in 1:npoints ){
  Pk=matrix(NA,max(ncat),nitems)
  for( j in 1:nitems ){
   Pk[1:ncat[j],j]=t( icrf[k,fromP[j]:toP[j]] )
  }
  Px_t=cbind( Px_t
              , sumsmnw( Pk, t(v), w, compress=0
                         , print=0, plot=0, debug=0 )[,2] )
 }
 rownames(Px_t)=0:maxscore_t; colnames(Px_t)=thname

 # joint
 Pxt=Px_t%*%Diag(Pt)
 Pxt[is.nan(Pxt)]=NA

 # marginal x
 Px=matrix(rowSums(Pxt),,1)

 return( Px )

} # end of obscore_s

obscore_s( param1, weight1, thmin=-4, thmax=4, npoints=31, thdist=1 )


