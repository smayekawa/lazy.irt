convP2PN <- function( param, print=0, DinP=1 ){
 # convert the original parametrization ( no 1.7 and step parameters)
 # to x-intercept parameters with 1.7
 # Shin-ichi Mayekawa
 # 20121015
 # DinP option:121109(London)
 # 121113
 # 
 #    ICRF of category k of item j at theta, namely, p_{jk}(theta)
 #    is defined as
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0,ncat[j]-1} Ez_{jk}(theta)
 #      where 
 #
 #    if type = "PN",
 #     Ez_{jk}(theta) = exp( 1.7 a_j k (theta - b_{jk}) )
 #                                              , k=0,1, ..., ncat[j]-1
 #    with b_{j0}=0.
 #
 #    or
 # 
 #    if type = "P",
 #     Ez_{jk}(theta) = exp( 1.7^DinP a*_j k (theta - sum_{h=0}^k b*_{jh}) )
 #                    = exp( 1.7^DinP a*_j k (theta - sum_{h=1}^k b*_{jh}) )
 #                                              , k=0,1, ..., ncat[j]-1
 #    with b*_{j0}=0.
 #
 #    Note that if DinP=1, 1.7 will be used.
 #
 #    This b*_{jk} is the original step parameter and it is 
 #    the value of theta where P_{jk-1} and P_{jk} intersect.
 #
 #
 locP=which( param$type == "P" )
 if( length(locP) == 0 )
  return( param )
 
 param0=param[locP,,drop=0]
 ncat=param0$ncat
 locparam=grep("^p[[:digit:]]$",colnames(param))
 if( length(locparam) >= 1 ) locparam=locparam[1]
 else locparam=grep("^a$",colnames(param))
 bs=param0[,(locparam+1):(locparam+max(ncat)-1),drop=0]
 b=bs
 b=t( apply( bs,1,cumsum)*(1/(1:ncol(bs))) )
 if( DinP == 0 ) param0[,locparam]=param0[,locparam]/1.7
 param0[,(locparam+1):(locparam+max(ncat)-1)]=b
 param0$type="PN"
 param[locP,]=param0
 
 if( print > 0 ){
  cat("\n item parameters with D in P =", DinP, "\n")
  Print(param0)
 }
 
 return( param )
 
} # end of convP2PN



convPN2P <- function( param, print=0, DinP=1 ){
 # convert the x-intercept parametrization ( with 1.7 and b-type parameters)
 # to the original step parameters w/o 1.7
 # Shin-ichi Mayekawa
 # 20121015
 # DinP option:121109(London)
 # 
 #    ICRF of category k of item j at theta, namely, p_{jk}(theta)
 #    is defined as
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0,ncat[j]-1} Ez_{jk}(theta)
 #      where 
 #
 #    if type = "PN",
 #     Ez_{jk}(theta) = exp( 1.7 a_j k (theta - b_{jk}) )
 #                                              , k=0,1, ..., ncat[j]-1
 #    with b_{j0}=0.
 #
 #    or
 # 
 #    if type = "P",
 #     Ez_{jk}(theta) = exp( 1.7^DinP a*_j k (theta - sum_{h=0}^k b*_{jh}) )
 #                    = exp( 1.7^DinP a*_j k (theta - sum_{h=1}^k b*_{jh}) )
 #                                              , k=0,1, ..., ncat[j]-1
 #    with b*_{j0}=0.
 #
 #    Note that if DinP=1, 1.7 will be used.
 #
 #    This b*_{jk} is the original step parameter and it is 
 #    the value of theta where P_{jk-1} and P_{jk} intersect.
 #
 #

 locPN=which( param$type == "PN" )
 if( length(locPN) == 0 )
  return( param )
 
 param0=param[locPN,,drop=0]
 ncat=param0$ncat
 locparam=grep("^p[[:digit:]]$",colnames(param))
 if( length(locparam) >= 1 ) locparam=locparam[1]
 else locparam=grep("^a$",colnames(param))
 b=param0[,(locparam+1):(locparam+max(ncat)-1)]
 bs=b
 for( j in 1:nrow(param0) ){
  for( k in 2:(ncat[j]-1) ){
   bs[j,k]=k*b[j,k]-(k-1)*b[j,k-1]
  }
 }
 
 if( DinP == 0 ) param0[,locparam]=1.7*param0[,locparam]
 param0[,(locparam+1):(locparam+max(ncat)-1)]=bs
 param0$type="P"
 param[locPN,]=param0
 
 if( print > 0 ){
  cat("\n item parameters with D in P =", DinP, "\n")
  Print(param0)
 }
  
 return( param )
 
} # end of convPN2P



convP2PN0 <- function( param, print=0, DinP=1 ){
 # convert the original parametrization ( no 1.7 and step parameters)
 # to a_j theta + c_jk type w/o 1.7
 # Shin-ichi Mayekawa
 # 20121122
 # 
 #    ICRF of category k of item j at theta, namely, p_{jk}(theta)
 #    is defined as
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0,ncat[j]-1} Ez_{jk}(theta)
 #      where 
 #
 #    if type = "PN",
 #     Ez_{jk}(theta) = exp( 1.7 a_j k (theta - b_{jk}) )
 #                                              , k=0,1, ..., ncat[j]-1
 #    with b_{j0}=0.
 #
 #    or
 # 
 #    if type = "P",
 #     Ez_{jk}(theta) = exp( 1.7^DinP a*_j k (theta - sum_{h=0}^k b*_{jh}) )
 #                    = exp( 1.7^DinP a*_j k (theta - sum_{h=1}^k b*_{jh}) )
 #                                              , k=0,1, ..., ncat[j]-1
 #    with b*_{j0}=0.
 #
 #    Note that if DinP=1, 1.7 will be used.
 #
 #    This b*_{jk} is the original step parameter and it is 
 #    the value of theta where P_{jk-1} and P_{jk} intersect.
 #
 #
 locP=which( param$type == "P" )
 if( length(locP) == 0 )
  return( param )
 
 param0=param[locP,,drop=0]
 ncat=param0$ncat
 locparam=grep("^p[[:digit:]]$",colnames(param))
 if( length(locparam) >= 1 ) locparam=locparam[1]
 else locparam=grep("^a$",colnames(param))
 bs=param0[,(locparam+1):(locparam+max(ncat)-1),drop=0]
 
 c=-t( apply( bs,1,cumsum) )*param0[,locparam]
 param0[,(locparam+1):(locparam+max(ncat)-1)]=c
 if( DinP == 1 ) 
  param0[,locparam:(locparam+max(ncat)-1)]=
      param0[,locparam:(locparam+max(ncat)-1)]*1.7
 param0$type="PN0"
 param0$type="PN"
 param[locP,]=param0
 
 if( print > 0 ){
  cat("\n item parameters with D in P =", DinP, "\n")
  Print(param0)
 }
 
 return( param )
 
} # end of convP2PN0



convP2N <- function( param, print=0, DinP=1 ){
 # convert the original parametrization ( no 1.7 and step parameters)
 # to a_jk theta + c_jk type nominal response model w/o 1.7
 # where a_jk = a_k * k
 # Shin-ichi Mayekawa
 # 20121122,1202
 # 
 #    ICRF of category k of item j at theta, namely, p_{jk}(theta)
 #    is defined as
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0,ncat[j]-1} Ez_{jk}(theta)
 #      where 
 #
 # 
 #    if type = "P",
 #     Ez_{jk}(theta) = exp( 1.7^DinP a*_j k (theta - sum_{h=0}^k b*_{jh}) )
 #                    = exp( 1.7^DinP a*_j k (theta - sum_{h=1}^k b*_{jh}) )
 #                                              , k=0,1, ..., ncat[j]-1
 #    with b*_{j0}=0.
 #
 #    Note that if DinP=1, 1.7 will be used.
 #
 #    This b*_{jk} is the original step parameter and it is 
 #    the value of theta where P_{jk-1} and P_{jk} intersect.
 #
 #
 
 locP=which( param$type == "P" )
 if( length(locP) == 0 )
  return( param )
 
 nitems=nrow(param)
 param0=param[locP,,drop=0]
 ncat=param0$ncat
 maxncat=max(ncat)
 locparam=grep("^p[[:digit:]]$",colnames(param))
 if( length(locparam) >= 1 ) locparam=locparam[1]
 else locparam=grep("^a$",colnames(param))
 bs=param0[,(locparam+1):(locparam+max(ncat)-1),drop=0]
 
 c=-t( apply( bs,1,cumsum) )*param0[,locparam]
 param0[,(locparam+1):(locparam+max(ncat)-1)]=c
 if( DinP == 1 ) 
  param0[,locparam:(locparam+max(ncat)-1)]=
     param0[,locparam:(locparam+max(ncat)-1)]*1.7
 a=param0[,locparam]
 param0$type="PN0"
 param[locP,]=param0
 
 np=max(2*(maxncat-1),max(param$ncat))
 param0=data.frame( param[,1:3]
        ,matrix(NA,nitems,np,dimnames=list(NULL,paste("p",1:np,sep="")))
                   )
 param0[,1:ncol(param)]=param
 for( j in 1:length(locP) ){
  param0[locP[j],4:(4+ncat[j]-2)]=a[j]*(1:(ncat[j]-1))
  param0[locP[j],(4+ncat[j]-1):(4+2*(ncat[j]-1)-1)]=
    param[locP[j],5:(5+ncat[j]-2)]
 }
 param0$type[locP]="N"
 param=param0
 
 if( print > 0 ){
  cat("\n item parameters with D in P =", DinP, "\n")
  Print(param0)
 }
 
 return( param )
 
} # end of convP2N





paramtest1=data.frame(
 name="Q1", type="B3", ncat=2,
 p1=1, p2=0, p3=0.3, p4=NA, p5=NA, p6=NA
 , stringsAsFactors=0 )

paramtest2=data.frame(
 name="Q2", type="P", ncat=2,
 p1=1, p2=0, p3=0, p4=NA, p5=NA, p6=NA
 , stringsAsFactors=0 )

paramtest3=data.frame(
 name="Q3", type="G", ncat=4,
 p1=1, p2=-1, p3=0, p4=1, p5=NA, p6=NA
 , stringsAsFactors=0 )

paramtest4=data.frame(
 name="Q4", type="P", ncat=5,
 p1=.7, p2=-1.5, p3=-1, p4=1, p5=1.5, p6=NA
 , stringsAsFactors=0 )

paramtest5=data.frame(
 name="Q5", type="P", ncat=6,
 p1=.7, p2=-1.5, p3=-1, p4=0, p5=1, p6=1.5
 , stringsAsFactors=0 )


paramtest=rbind(paramtest1, paramtest2, paramtest3, paramtest4, paramtest5)



theta=seq(-4,4,length.out=21)

comments('
temp=irf(paramtest4,theta,plot=1,print=0)

temp=icrfPN0( convP2PN0(paramtest4), theta, print=0, plot=1)

temp=irf(paramtest4,theta,plot=1,print=0)
')

temp=irf(paramtest2,theta,plot=1,print=0)
temp=irf(convP2N(paramtest2),theta,plot=1,print=0)


