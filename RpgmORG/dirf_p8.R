#' Calculation of the Derivatives of the Item Category Response Function
#' with respect to Item Parameters
#'
#' @param paramj item parameters data frame for ONE item
#' @param theta Discrete theta values
#' @param weight Weight data frame:  NOT used.
#' @param smallP Minimum value of probability
#' @param DinP = 1 to include D=1.7 in logistic function
#' @param thmin Minimum value of discrete thata value
#' @param thmax Maximum value of discrete thata value
#' @param npoints # of discrete points for theta
#' @param Pj icrf:   npoints x (ncatj-1)  (no zero category)  or NULL
#' @param PPj icbrf of the Graded Response Model   or NULL
#' @param zero = 1 to include the zero-th category in output
#' @param cat.first  = 1 to chage the category fist in the rows of Jack.
#' @param log = 1 to obtain the log Jacobian: d log(ICRF) d param
#' @param print = 1 to print result
#'
#' @return
#'   list of (Jack, Pj, PPj) \cr
#'   Jack     (length(theta) x (ncatj-1)) x ncatj \cr
#'             derivative of vec(Pj) with respect to (a, b1, b2, ...) \cr
#'           If cat.first = 0 \cr
#'             theta changes first, then k changes from 1 to ncatj \cr
#'           If cat.first = 1, category(k) changes first, then theta. \cr
#'   Pj Item Category Response Function (icrf)  (length(theta) x (ncatj-1))\cr
#'   PPj will be output when type="G" or "Gn".
#'
#'
#' @examples
#' res=dirf_p( paramS1[3,], npoints=5, print=1 )
#' res=dirf_p( paramS1[3,], npoints=5, print=1, zero=1 )
#' res=dirf_p( paramS1[3,], npoints=5, print=1, cat.first=1 )
#'
#'
#' @export
#'

dirf_p <- function( paramj, theta=NULL, weight=NULL, smallP=0, DinP=1
                  , thmin=-4, thmax=4, npoints=21
                  , Pj=NULL, PPj=NULL
                  , zero=0, cat.first=0, log=0, print=0 ){
 # derivative of ICRF with respect to item paramters
 # Shin-ichi Mayekawa
 # 121013
 # zero=1: 121018
 # cat.first: 121022
 # log: 121022,24
 # smallP: 121022
 # calculation of PPj from Pj for Graded: 121105
 # type = "P": 121106,08(London)
 # DinP: 121109(London)
 # Check Pj in denominator by smallP: 121115
 # Check Pj in denominator by 1e-7: 121119
 # P with ncat=2: 121122
 # nominal response model: 121123
 # bugfix for type "B": 20150614
 # error check: 201506014
 # rownames for Jack: 20150615
 # check if c is na or null: 20150615
 # Bn and Gn: 20180211
 # roxgen2 comment: 20180212
 #
 #
 #
 # Args:
 #  paramj   item parameters data frame for one item
 #  theta    lenth x 1 theta values
 #  zero     = 1 to include the 0-th category
 #  Pj       icrf:   npoints x (ncatj-1)  (no zero category)
 #  PPj      icbrf for the Graded Response Model
 #  cat.first = 1 to chage the category fist in the rows of Jack.
 #           ( 10000 theta points or less.)
 #  log      = 1 to obtain the log Jacobian: d log(ICRF) d param
 #  smallp   = small value for icrf etc
 #
 # Values:
 #   list of (Jack, Pj, PPj)
 #   Jack     (length(theta) x (ncatj-1)) x ncatj
 #             derivative of vec(Pj) with respect to (a, b1, b2, ...)
 #           If cat.first = 0
 #             theta changes first, then k changes from 1 to ncatj
 #           If cat.first = 1, category(k) changes first, then theta.
 #
 # Needs:
 #   icrfB, icrfG, icrfPN
 #


 JacobianMat <- function( parvec, infcn, eps = 1e-06 ){
  # Function to calculate the difference-quotient approx gradient matrix
  # of an arbitrary input function infcn.
  # Now recoded to use central differences !
  #
  # Original version was Gradmat in Sec4NotF09.pdf
  #

  dd = length(parvec)
  aa = length(infcn(parvec))
  epsmat = (diag(dd) * eps)/2
  gmat = array(0, dim = c(aa, dd))
  for(i in 1:dd)
   gmat[, i]=( infcn(parvec + epsmat[, i] )
               - infcn(parvec - epsmat[, i] ) )/eps

  # return a vector if infcn is unidimensional
  if( aa == 1 ) gmat=c(gmat)
  return( gmat )

 } # end of JacobianMat


 icrfP00 <- function( paramnum ){
  # This is to be used in conjunction with JacobianMat
  paramj0=paramj
  paramj0[,4:(4+length(paramnum)-1)]=paramnum
  res=irf( paramj0, theta, zero=0, print=0 )$ICRF
  res=matrix(res,,1)
  # vertically stack all categories
  return(  matrix( res, , 1 )  )
 } # end of icrfP00


 # error check
 if( nrow(paramj) > 1 ){
  cat("\n\nerror1:(dirf_p)  paramj data frame contains more than 2 items.\n ")
  return()
 }

 # const
 ncatj=paramj$ncat
 ncatj1=ncatj-1
 typej=paramj$type

 # # of parameters to be estimated
 if( typej %in% c("B3","Bn3") ) np=ncatj+1
 else if( typej == "N" ) np=2*ncatj1
 else np=ncatj

 if( is.null(theta) ){
  theta=seq(thmin,thmax,length.out=npoints)
 }
 else{
  if( is.matrix(theta) ) theta=as.vector(theta)
  thmin=theta[1]; thmax=theta[length(theta)]
 }
 lenth=length(theta)

 if( DinP == 1 ) DD=1.7
 else DD=1

 # item param for type Bx, G, P, PN
 locparam=grep("^p[[:digit:]]$",colnames(paramj))
 if( length(locparam) >= 1 ) locparam=locparam[1]
 else locparam=grep("^a$",colnames(paramj))
 a=paramj$p1
 b=as.matrix(paramj[1,(locparam+1):(locparam+ncatj1)],1)
 if( typej %in% c("G","Gn") ) b=cbind(b,9999)
 c=paramj$p3
 # 20150615
 if( is.null(c) ) c=0; if( is.na(c) ) c=0

 # storage
 Jack=matrix(0,ncatj1*lenth,np)
 colnames(Jack)=paste("p",1:np,sep="")
 # thname=paste( "t", 1:lenth, sep="" )
 thname=format(theta,digits=3)
 catname=paste( "c", 1:ncatj1, sep="" )
 rnJ=apply( expand.grid( thname, catname ), 1, paste, collapse="-" )
 rownames(Jack)=rnJ

  # icrf at current paramnum
 if( length(grep("^B[[:digit:]]*$", typej)) > 0 ){

  if( is.null(Pj) ){
   # icrf for k=1
   Pj=icrfB( paramj, theta, smallP=smallP, zero=0 )
  }
  # Jacobian
  Jack[,1]=1.7/(1-c)*(theta-b[1])*(Pj-c)*(1-Pj)
  Jack[,2]=-1.7/(1-c)*a*(Pj-c)*(1-Pj)
  if( typej == "B3" ) Jack[,3]=(1-Pj)/(1-c)

 } # end of type = "B"
 else if( length(grep("^Bn[[:digit:]]*$", typej)) > 0 ){

  dd=dnorm(a*(theta-b[1]))
  if( is.null(Pj) ){
   # icrf for k=1
   Pj=icrfBn( paramj, theta, smallP=smallP, zero=0 )
  }
  # Jacobian
  Jack[,1]=dd*(theta-b[1])*(1-c)
  Jack[,2]=-dd*a*(1-c)
  if( typej == "Bn3" ) Jack[,3]=(1-Pj)/(1-c)

 } # end of type = "Bn"
 else if( typej == "G" ){

  if( is.null(Pj) ){
   # icrf : no zero-th category
   temp=icrfG( paramj, theta, smallP=smallP )
   Pj=temp$P[,2:ncatj,drop=0];
   PPj=temp$PP[,2:ncatj,drop=0]; PPj=cbind(PPj,0)
  }
  else if( is.null(PPj) ){
   # calculation of icbrf from icrf
   PPj=matrix( 0,lenth, ncatj)
   for( k in ncatj1:1 ){
    PPj[,k]=Pj[,k] + PPj[,k+1]
   }
  }
  # Jacobian
  for( k in 1:ncatj1 ){
   temp=1.7*(
    (theta-b[k])*PPj[,k]*(1-PPj[,k])-(theta-b[k+1])*PPj[,k+1]*(1-PPj[,k+1]) )
   Jack[((k-1)*lenth+1):(k*lenth),1]=temp
   for( z in 1:ncatj1 ){
    temp=0
    if( k == z ) temp=-1.7*a*PPj[,z]*(1-PPj[,z])
    else if( (k+1) == z ) temp=1.7*a*PPj[,z]*(1-PPj[,z])
    Jack[((k-1)*lenth+1):(k*lenth),1+z]=temp
   } # z
  } # k

 } # end of type = "G"
 else if( typej == "Gn" ){

  maxZ=700
  Z=-(a*outer(c(b),theta,"-"))
  Z[which(Z > maxZ)]=maxZ
  Z[which(Z < -maxZ)]=-maxZ
  dPj=t( dnorm(Z) )

  if( is.null(Pj) ){
   # icrf : no zero-th category
   temp=icrfGn( paramj, theta, smallP=smallP )
   Pj=temp$P[,2:ncatj,drop=0];
   PPj=temp$PP[,2:ncatj,drop=0]; PPj=cbind(PPj,0)
  }
  else if( is.null(PPj) ){
   # calculation of icbrf from icrf
   PPj=matrix( 0,lenth, ncatj)
   for( k in ncatj1:1 ){
    PPj[,k]=Pj[,k] + PPj[,k+1]
   }
  }
  # Jacobian
  for( k in 1:ncatj1 ){
   temp=( (theta-b[k])*dPj[,k]-(theta-b[k+1])*dPj[,k+1] )
   Jack[((k-1)*lenth+1):(k*lenth),1]=temp
   for( z in 1:ncatj1 ){
    temp=0
    if( k == z ) temp=-a*dPj[,z]
    else if( (k+1) == z ) temp=a*dPj[,z]
    Jack[((k-1)*lenth+1):(k*lenth),1+z]=temp
   } # z
  } # k

 } # end of type = "Gn"
 else if( typej == "PN" ){

  if( is.null(Pj) ){
   # icrf : no zero-th category
   Pj=icrfPN( paramj, theta, smallP=smallP )$P[,2:ncatj,drop=0]
  }
  # Jacobian
  thetab=matrix(theta,lenth,ncatj1)-matrix(1,lenth,1)%*%b
  Pqthetab=rowSums( (Pj%*%diag(1:ncatj1))*thetab )
  for( k in 1:ncatj1 ){
   temp=1.7*Pj[,k]*(k*thetab[,k]-Pqthetab)
   Jack[((k-1)*lenth+1):(k*lenth),1]=temp
   for( z in 1:ncatj1 ){
    if( k == z ) temp=-1.7*a*z*Pj[,z]*(1-Pj[,z])
    else  temp=1.7*a*z*Pj[,k]*Pj[,z]
    Jack[((k-1)*lenth+1):(k*lenth),1+z]=temp
   } # z
  } # k

 } # end of type = "PN"
 else if( typej == "P" ){

  if( is.null(Pj) ){
   # icrf : no zero-th category
   Pj=icrfP( paramj, theta, smallP=smallP, DinP=DinP )$P[,2:ncatj,drop=0]
  }
  # Jacobian
  thetab=matrix(theta,lenth,ncatj1)-matrix(1,lenth,1)%*%b
  thetabc=matrix( t( apply(thetab,1,cumsum) ), lenth )
  Pqthetabc=rowSums( Pj*thetabc )
  # reverse cum sum
  Pqac=DD*a*matrix( t(apply(Pj[,ncatj-(1:ncatj1),drop=0],1,cumsum)) ,lenth)
  Pqac=Pqac[,ncatj-(1:ncatj1),drop=0] # reverse again
  for( k in 1:ncatj1 ){
   temp=DD*Pj[,k]*(thetabc[,k]-Pqthetabc)
   Jack[((k-1)*lenth+1):(k*lenth),1]=temp
   for( z in 1:ncatj1 ){
    if( k >= z ) temp=-DD*a + Pqac[,z]
    else  temp=Pqac[,z]
    Jack[((k-1)*lenth+1):(k*lenth),1+z]=Pj[,k]*temp
   } # z
  } # k
  #  Jack=JacobianMat( c(a,b), icrfP00, eps = 1e-06 )
  #  Print(Jack)

 } # end of type = "P"
 else if( typej == "N" ){

  if( is.null(Pj) ){
   # icrf : no zero-th category
   Pj=icrfN( paramj, theta, smallP=smallP, DinP=DinP )$P[,2:ncatj,drop=0]
  }

  # Jacobian
  for( k in 1:ncatj1 ){
   for( z in 1:ncatj1 ){
    if( k == z ){
     Jack[((k-1)*lenth+1):(k*lenth),k]=theta*Pj[,k]*(1-Pj[,k])
     Jack[((k-1)*lenth+1):(k*lenth),ncatj1+k]=Pj[,k]*(1-Pj[,k])
    }
    else{
     Jack[((k-1)*lenth+1):(k*lenth),z]=-theta*Pj[,k]*Pj[,z]
     Jack[((k-1)*lenth+1):(k*lenth),ncatj1+z]=-Pj[,k]*Pj[,z]
    }
   }
  }

 } # end of type = "N"



 if( zero == 1 ){
  # dP0dz = d(1-sum_k Pj)dz = 0 - sum_k dPjdz
  temp=0
  rntemp=paste(thname,"-c0",sep="")
  for( k in 1:ncatj1 ){
   temp=temp+Jack[((k-1)*lenth+1):(k*lenth),]
  } # k
  Jack=rbind(-temp, Jack)
  rownames(Jack)=c(rntemp,rnJ)
  Pj=cbind(1-apply(Pj,1,sum),Pj)
  colnames(Pj)[1]="0"
 }

 if( log == 1 ){
  # Pj[which(Pj < smallP)]=smallP
  # Pj[which(Pj < 1e-307)]=1e-307
  # Pj[which(Pj < 1e-30)]=1e-30
  # This 1e-7 seems to be the best.
  Pj[which(Pj < 1e-7)]=1e-7
  Jack=Jack/as.vector(matrix(Pj,,1))
 }

 if( cat.first ){
  # sort the rows of Jack to make the category change first.
  ncat1=ncatj-1; if( zero == 1 ) ncat1=ncatj
  index=10000*rep(1:lenth,ncat1) + as.vector(
                  sapply(1:ncat1, function(x) rep(x,lenth)) )
  Jack=Jack[order(index),,drop=0]
 }

 if( print > 0 ){
  Print(lenth, zero, log, cat.first, smallP,paramj)
  Print(Pj,Jack, fmt="5.2")
 }

 return( list(Jack=Jack, Pj=Pj, PPj=PPj) )

} # end of dirf_p






comment(
'

paramtest1=data.frame( name="Q1", type="N", ncat=4
                       , p1=1.19, p2=2*1.19, p3=3*1.19
                       , p4=1.19, p5=1.19, p6=0, p7=NA, p8=NA
                       , stringsAsFactors=0)
paramtest2=data.frame( name="Q1", type="N", ncat=5
                       , p1=1.2, p2=2*1.2, p3=3*1.2, p4=4*1.2
                       , p5=1.785, p6=2.975, p7=1.785, p8=0
                       , stringsAsFactors=0)

res=dirf_p( paramtest2, thmin=-4, thmax=4, npoints=5, print=1, zero=0 )





####  testing type "Bn"
thmax=3
thmin=-thmax
theta=seq(thmin,thmax,length=5)

param=paramS1[2,]
param$p1=1
param$p3=0.01
param$type="Bn"

res=dirf_p( param, theta, print=1, zero=0 )




#### testing type "Gn"
param=paramS2[5:6,]
param$type="Gn"

res=dirf_p( param[2,], theta, print=1, zero=0 )


')

