gendataIRT <- function( Ntot, param, DinP=1
                      , theta=NULL, npoints=31, thmin=-3, thmax=3
                      , thd=NULL, thdist="NORMAL", thmean=0, thstd=1
                      , nomiss=0, compress=0, sort=1 ){
 # generation of Item response data
 # Shin-ichi Mayekawa
 # 121017,18,24
 # DinP: 121110(London)
 # renamed as gendataIRT: 121110(London)
 # make the data integer when possible: 121126
 # bugfix: 20150614
 # nomiss: 20150615
 #
 #
 # Args:
 #
 #   npoints  # of theta points:  or the length of theta
 #
 #   Ntot     # of total observations to be generated
 #            If Ntot == 1, one observation per each theta point,
 #            else the # of observations per each theta point is
 #            proportional to round(Ntot*thd*) where thd is
 #            the npoint x 1 vector or theta distribution.
 #            Therefore, Ntot must be large enough.
 #   param    parameter data frame
 #   theta    theta points
 #            If null, npoints, thmin and thmax will be used to generate it
 #   thd      theta distribution probability vector
 #   thdist   = "NORMAL" or "UNIFORM" or "rnorm"  or  "runif"
 #            When "NORMAL" or "rnorm",
 #                thmean and thstd will be used to generate thd.
 #            When "UNIFORM" or "runif",
 #                thmin and thmax will be used to generate thd.
 #            When "rnorm" or "runif", theta will be generated using
 #                npoints random numbers and thd and Ntot is set equal to 1.
 #   nomiss   = 1 to omit those theta points with N=0.
 #   compress = 1 to compress the data when Ntot = 1.
 #            See the Values section.
 #
 #
 # Values:
 #            list of U, fromP, toP, N, theta, thd
 #
 #            If compress = 0 and Ntot > 1,
 #              U    is   npoints x sum of ncat[j]
 #               U[,fromP[j]:toP[j]]  is npoints x ncat[j]
 #               U[,fromP[j]:toP[j]][i,k] contains # of responses to
 #               the k-th category (k=0,1,...,ncat[j]) at theta[i]
 #               to which N[i]=round(Ntot*thd[i]) observations belong.
 #
 #
 #            If compress = 1 and Ntot = 1,
 #              U    is   npoints x nitems
 #               U[i,j] contains the response to the j-th item at theta[i]
 #               where 0 <= U[i,j] <= ncat[j]-1.
 #              fromP and toP are not compressed.
 #
 #            (theta, thd, N)   are   npoints x 1
 #            N[i]  is the # of observations at theta[i]
 #
 #            Response probability out of Ntot can be obtained as
 #              P = U/N
 #
 #            U will be of integer type.
 #

 # const
 nitems=nrow(param)

 # generate theta
 if( is.null(theta) ){
  theta=seq(thmin,thmax,length.out=npoints)
 }
 else{
  npoints=length(as.vector(theta))
 }

 # error check
 if( toupper(thdist) == "NORMAL"  ||  toupper(thdist) == "UNIFORM" )
  if( Ntot < 3*npoints ){
   cat("\n\nwarning1(gendataIRT): Ntot may be too small. \n\n")
  }

 if( Ntot > 1  &  compress == 1 ){
  cat("\n\nwarning1(gendataIRT): compress=1 cannot be used when
      Ntot>1. \n\n")
 }

 # theta distribution
 if( is.null(thd) ){
  if( toupper(thdist) == "NORMAL" ){
   thd=exp(-0.5*( (theta-thmean)/thstd )^2); thd=thd/sum(thd)
  }
  else if( toupper(thdist) == "UNIFORM" ){
   thd=1/npoints
  }
  else if( toupper(thdist) == "RNORM" ){
   thd=1; Ntot=1;
   theta=rnorm(npoints,thmean,thstd)
  }
  else if( toupper(thdist) == "RUNIF" ){
   thd=1; Ntot=1;
   theta=runif(npoints,thmin,thmax)
  }
  else{
   thd=1;
   Ntot=1
  }
 }

 if( sort == 1 ) theta=sort(theta)

 # # of obs. per theta point
 if( Ntot > 1 ){
  N=round(Ntot*thd)
  if( sum(N) < Ntot ){
   N[round(npoints/2):(round(npoints/2)-(Ntot-sum(N))+1)]=
                N[round(npoints/2):(round(npoints/2)-(Ntot-sum(N))+1)]+1
  }
 }
 else N=rep(1,npoints)

 # generate irf
 res=irf( param, theta, DinP=DinP, zero=1, print=0, plot=0 )
 icrf=res$ICRF; fromP=res$fromP; toP=res$toP; ncat=res$maxscore_i+1
 rm(res)
 # Print(icrf, fromP, toP, ncat)

 # generate dummy expanded response according to icrf
 Y=matrix(NA,nrow(icrf),ncol(icrf))
 dimnames(Y)=dimnames(icrf)
 for( j in 1:nitems ){
  domj=c(0:(ncat[j]-1))
  Pj=icrf[,fromP[j]:toP[j]]
  # Print(ncat[j],Pj,digits=3)
  for( i in 1:npoints ){
   Y[i,fromP[j]:toP[j]]=t( rmultinom( 1, N[i], Pj[i,] ) )
  }
 }

 if( nomiss == 1 ){
  # delete rows with N=0
  loc=which( N > 0 )
  N=N[loc]; theta=theta[loc]; thd=thd[loc]
  icrf=icrf[loc,]; Y=Y[loc,]
  npoints=length(loc)
 }

 # compress if requested
 if( Ntot == 1 && compress == 1 ){
  U=matrix(NA,npoints,nitems)
  for( j in 1:nitems ){
   loc=apply( Y[,fromP[j]:toP[j]], 1, function(y) which( y == 1) )
   U[,j]=loc-1
  }
  Y=matrix(as.integer(U),nrow(U),ncol(U))
  rownames(Y)=rownames(icrf)
  colnames(Y)=param$name
 }

 return( list( U=Y, N=N, npoints=npoints, theta=theta, thd=thd
             , Ntot=Ntot, fromP=fromP, toP=toP
             , thmean=thmean, thstd=thstd
             , compress=compress, sort=sort ) )

} # end of gendataIRT




dummy_expand <- function( Uc, ncat=NULL, type=NULL, zero=1 ){
 # dummy expand the categorical valued variable u
 # Shin-ichi Mayekawa
 # 20121022
 # matrix input: 121024
 # missing value: 121025
 # make the data integer when possible: 121126
 # vectorize fromP: 20151113

 nitems=ncol(Uc)
 npoints=nrow(Uc)
 iname=colnames(Uc)

 if( is.null(ncat) ){
  ncat=matrix(apply(Uc,2,max,na.rm=1),,1)+1
 }
 toP=cumsum(ncat)
 fromP=c(toP-ncat+1)
 # Print(nitems,npoints,ncat,toP,fromP, iname)
 U=matrix(as.integer(0),npoints,toP[nitems])
 for( j in 1:nitems ){
  uj=Uc[,j]
  locobs=!is.na(uj)
  min=0; max=ncat[j]-1
  factuj=factor(uj, levels=min:max)
  Uj=model.matrix(~factuj-1)
  #Print(uj,factuj,Uj)
  #  attributes(Uj)$assign=NULL
  #  attributes(Uj)$contrasts=NULL
  #  colnames(Uj)=min:max
  U[locobs,fromP[j]:toP[j]]=as.integer(Uj)
 }
 if( is.null(iname) )
  iname=paste("Q",1:nitems,sep="")
 Pname=character(sum(ncat))
 catname=outer(paste(iname,"_",sep=""),0:max(ncat),paste,sep="")
 for( j in 1:nitems ){
  Pname[fromP[j]:toP[j]]=catname[j,1:ncat[j]]
 }
 colnames(U)=Pname
 rownames(U)=rownames(Uc)

 if( zero != 1 ){
  # remove the 0-th category
  U=U[,-fromP]
  fromP=c(fromP-0:(nitems-1))
  toP=toP-1:nitems
 }

 return( list(U=U, ncat=ncat, type=type, zero=zero, fromP=fromP, toP=toP) )

} # end of dummy_expand





set.seed(1701)
# 100 observations tabulated at seven theta points in [-3,3]
res <- gendataIRT( 1000, paramS1, npoints=7, thmin=-3, thmax=3 )
Print(res$N,res$U,sum(res$N))

# seven theta points in [-3,3]
set.seed(1701)
res1 <- gendataIRT( 1, paramS1, npoints=7, thmin=-3, thmax=3 )
Print(res1$N,res1$U,sum(res1$N))

# seven theta points in [-3,3] and compressed
set.seed(1701)
res1 <- gendataIRT( 1, paramS1, npoints=7, thmin=-3, thmax=3, compress=1 )
Print(res1$N,res1$U,sum(res1$N))

# 10 observations with theta from N( 0, 1 )
set.seed(1701)
res1 <- gendataIRT( 1, npoints=10, paramS1, thdist="rnorm", compress=1 )
Print(res1$N,res1$U,sum(res1$N))

comments('
res=gendataIRT( 10, paramS1, thd=NULL
             , npoints=31, thmin=-3, thmax=3, compress=0 )
U=res$U
N=res$U
P=U/res$N


check1=mapply( function( from,to) apply(P[,from:to],1,sum), fromP,toP )
check2=mapply( function( from,to) apply(U[,from:to],1,sum), fromP,toP )
Print(check1,check2)


res=genndataIRT( 1, param, thdist=""
             , npoints=31, thmin=-3, thmax=3, compress=0 )
U=res$U
P=U/res$N
check1=mapply( function( from,to) apply(P[,from:to],1,sum), fromP,toP )
check2=mapply( function( from,to) apply(U[,from:to],1,sum), fromP,toP )
Print(check1,check2)


res=genndataIRT( 1, param
                 , npoints=31, thmin=-3, thmax=3, compress=1 )



res=genndataIRT( 1000, param[19,]
                 , npoints=11, thmin=-3, thmax=3, compress=0 )
U=res$U
P=U/res$N


res=genndataIRT( 1, param[19,], thdist="rnorm", thd=NULL
                 , npoints=110, thmin=-3, thmax=3, compress=1 )

')



comments('
paramtest1=data.frame( name="Q1", type="N", ncat=4
                       , p1=1.19, p2=2*1.19, p3=3*1.19
                       , p4=1.19, p5=1.19, p6=0, p7=NA, p8=NA
                       , stringsAsFactors=0)
paramtest2=data.frame( name="Q2", type="N", ncat=5
                       , p1=1.2, p2=2*1.2, p3=3*1.2, p4=4*1.2
                       , p5=1.785, p6=2.975, p7=1.785, p8=0
                       , stringsAsFactors=0)
paramtest=rbind(paramtest1,paramtest2)

theta=seq(-4,4,length.out=11)

# temp=irf(paramtest,theta,plot=1,print=0)

resg=gendataIRT( 1, paramtest, theta=theta, thdist="NORMAL", thd=NULL
                 , thmean=0, thstd=1, compress=1 )
testdata=resg$U

Print(typeof(testdata),"/",testdata)

testdata1=dummy_expand( testdata )
Print(typeof(testdata1$U),"/",testdata1$U)

')
