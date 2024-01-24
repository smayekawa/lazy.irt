#' Generation of Simulated Item Response Data
#'
#' @param Ntot # of total observations to be generated
#'            If Ntot == 1, one observation per each theta point
#'            will be generated, \cr
#'            otherwize, the # of observations per each theta point is \cr
#'            proportional to round(Ntot*thd) where thd is \cr
#'            the npoint x 1 vector of theta distribution. \cr
#'            Therefore, Ntot must be large enough. \cr
#' @param param Item Parameter Data Frame
#' @param DinP = 1 to include D=1.7 in logistic function
#' @param Nmat npoints x nitems matrix consisting of # of trials, or NULL.
#' This will be used only when Ntot=1 and compress=0.
#' @param zero = 0 to exclude zero-th category when not compressed.
#' @param theta Discrete theta values
#' @param npoints # of discrete points for theta, or # of random variables
#' @param thmin Minimum value of discrete thata value
#' @param thmax Maximum value of discrete thata value
#' @param thd theta distribution probability vector
#' @param thdist = "NORMAL" or "UNIFORM" or "rnorm"  or  "runif" \cr
#'            When "NORMAL" or "rnorm", \cr
#'                thmean and thstd will be used to generate thd. \cr
#'            When "UNIFORM" or "runif", \cr
#'                thmin and thmax will be used to generate thd. \cr
#'            When "rnorm" or "runif", theta will be generated using \cr
#'                npoints random numbers and thd  and Ntot is set equal to 1.
#' @param thmean mean of theta distribution
#' @param thstd std of theta distribution
#' @param nomiss = 1 to exclude those theta points with N = 0.
#' @param compress = 1 to compress the data when Ntot = 1.
#' @param sort = 1 to sort the result according to theta.
#' @param thetaname = Name of the variable which will contain theta value
#' or NULL.
#'
#' @return A list of \cr
#'  U, N, npoints, theta, thd, Ntot, fromP, toP, type, ncat, thmean, thstd
#'  , compress, sort
#'
#' \cr
#' where
#' \cr
#' \cr
#'            If compress = 0 and Ntot is large, \cr
#'              U    is   npoints x sum of ncat[j] \cr
#'               U[,fromP[j]:toP[j]]  is npoints x ncat[j] \cr
#'               U[,fromP[j]:toP[j]][i,k] contains # of responses to \cr
#'               the k-th category (k=0,1,...,ncat[j]) at theta[i] \cr
#'               to which N[i]=round(Ntot*thd[i]) observations belong. \cr
#' \cr
#' \cr
#'            If compress = 1 and Ntot = 1, \cr
#'              U    is   npoints x nitems \cr
#'               U[i,j] contains the response to the j-th item at theta[i] \cr
#'               where 0 <= U[i,j] <= ncat[j]-1. \cr
#'              In this case, fromP and toP are not compressed. \cr
#' \cr
#'            (theta, thd, N)   are   npoints x 1 \cr
#'            N[i]  is the # of observations at theta[i] \cr
#'
#' @examples
#' # 10 observations with theta from N( 0, 1 ): not compressed
#' set.seed(1702)
#' res1 <- gendataIRT( 1, paramS1, npoints=10, thdist="rnorm" )
#' Print(res1$N,res1$U,res1$theta)
#' #
#' # 10 observations with theta from N( 0, 1 ): compressed
#' set.seed(1702)
#' res1 <- gendataIRT( 1, paramS1, npoints=10, thdist="rnorm", compress=1 )
#' Print(res1$N,res1$U,res1$theta)
#' #
#' # 10 observations with theta from N( 0, 1 ): no zero category
#' set.seed(1702)
#' res1 <- gendataIRT( 1, paramS1, npoints=10, thdist="rnorm", zero=0 )
#' Print(res1$N,res1$U,res1$theta)
#' #
#' #
#' #
#' # seven theta points in [-3,3]
#' set.seed(1701)
#' res1 <- gendataIRT( 1, paramS1, npoints=7, thmin=-3, thmax=3 )
#' Print(res1$N,res1$U,res1$theta)
#' #
#' # seven theta points in [-3,3] and compressed
#' set.seed(1701)
#' res1 <- gendataIRT( 1, paramS1, npoints=7, thmin=-3, thmax=3, compress=1 )
#' Print(res1$N,res1$U,res1$theta)
#' #
#' set.seed(1701)
#' # 100 observations tabulated at seven theta points in [-3,3]
#' res1 <- gendataIRT( 1000, paramS1, npoints=7, thmin=-3, thmax=3 )
#' Print(res1$N,res1$U,res1$theta)
#' #
#' #
#' # Using Nmat argument
#' set.seed(1701)
#' npoints=20
#' param=paramS1
#' nitems=nrow(param)
#' Nmat=matrix(100,npoints,nitems)
#' res1 <- gendataIRT( 1, param, Nmat=Nmat, npoints=npoints
#'                   , thdist="rnorm", compress=1 )
#' Print(res1$N,res1$U,res1$theta)
#'
#'
#' @export
#'

gendataIRT <- function( Ntot, param, DinP=1, Nmat=NULL, zero=1
                      , theta=NULL, npoints=31, thmin=-3, thmax=3
                      , thd=NULL, thdist="NORMAL", thmean=0, thstd=1
                      , nomiss=0, compress=0, sort=1, thetaname=NULL ){
 # generation of Item response data
 # Shin-ichi Mayekawa
 # 121017,18,24
 # DinP: 121110(London)
 # renamed as gendataIRT: 121110(London)
 # make the data integer when possible: 121126
 # bugfix: 20150614
 # nomiss: 20150615
 # theta as rownames abandoned: 20161103
 # Nmat: 20161103
 # zero: 20170303
 # bugfix: 20170307
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
 #   zero     = 0 to exclude zero-th category
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
 type=param$type
 ncat=param$ncat

 # generate theta
 if( is.null(theta) ){
  theta=seq(thmin,thmax,length.out=npoints)
 }
 else{
  npoints=length(as.vector(theta))
 }

 # error check
 if( toupper(thdist) == "NORMAL"  ||  toupper(thdist) == "UNIFORM" )
  if( Ntot > 1  &&  Ntot < 3*npoints ){
   cat("\n\nwarning1(gendataIRT): Ntot may be too small. \n\n")
  }

 if( Ntot > 1  &  compress == 1 ){
  cat("\n\nwarning1(gendataIRT): compress=1 cannot be used when
      Ntot > 1 or Nmat is given. \n\n")
 }

 if( !is.null(Nmat) ){
  Nmatin=1
  if( Ntot != 1 ){
   cat("\n\nerror1(gendataIRT): Ntot must be 1 when Nmat is given.\n\n")
   return()
  }
  if( npoints != nrow(Nmat) || nrow(param) != ncol(Nmat) ){
   cat("\n\nerror1(gendataIRT): Nmat does not comform to param.\n\n")
   return()
  }
 }
 else Nmatin=0

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

 if( !( Ntot == 1  && !is.null(Nmat) ) ){
  Nmat=matrix(N,length(N),nitems)
 }

 # generate irf
 res=irf( param, theta, DinP=DinP, zero=1, print=0, plot=0 )
 icrf=res$ICRF; fromP=res$fromP; toP=res$toP; ncat=res$maxscore_i+1
 thname=rownames(icrf)
 rm(res)
 # Print(icrf, fromP, toP, ncat)

 # generate dummy expanded response according to icrf
 Y=matrix(NA,nrow(icrf),ncol(icrf))
 colnames(Y)=colnames(icrf)
 rownames(Y)=1:nrow(Y)
 for( j in 1:nitems ){
  domj=c(0:(ncat[j]-1))
  Pj=icrf[,fromP[j]:toP[j]]
  # Print(ncat[j],Pj,digits=3)
  for( i in 1:npoints ){
   Y[i,fromP[j]:toP[j]]=t( rmultinom( 1, Nmat[i,j], Pj[i,] ) )
  }
 }

 if( nomiss == 1 ){
  # delete rows with N=0
  locobs=which( rowSums(Nmat) > 0 )
  N=N[locobs]; theta=theta[locobs]; thd=thd[locobs]
  icrf=icrf[locobs,]; Y=Y[locobs,]
  npoints=length(locobs)
  Nmat=Nmat[locobs,]
 }

 # compress if requested
 if( Ntot == 1 && compress == 1 && all( Nmat==1 ) ){
  U=matrix(NA,npoints,nitems)
  for( j in 1:nitems ){
   loc=apply( Y[,fromP[j]:toP[j]], 1, function(y) which( y == 1) )
   U[,j]=loc-1
  }
  U=matrix(as.integer(U),nrow(U),ncol(U))
  colnames(U)=param$name
 }
 else U=Y
 rm(Y)

 rownames(U)=1:nrow(U)

 if( !is.null(thetaname) ){
  U=cbind(U, theta)
  colnames(U)[ncol(U)]=thetaname
 }

 # zero
 if( zero == 0 ){
  loc=setdiff(1:ncol(U),grep("_0",colnames(U)))
  U=U[,loc]
 }

 res=list( U=U, N=N, npoints=npoints, theta=theta, thd=thd, icrf=icrf
           , Ntot=Ntot, fromP=fromP, toP=toP, type=type, ncat=ncat
           , thmean=thmean, thstd=thstd
           , compress=compress, sort=sort, zero=zero )
 if( Nmatin ) res=append(res,list(Nmat=Nmat))
 return( res )

} # end of gendataIRT


