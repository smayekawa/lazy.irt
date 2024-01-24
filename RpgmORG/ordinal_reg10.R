#' Oridinal Regression
#'
#' This program regresses R on x where
#' R is the length(x) x ncat matrix.
#'
#' @param R A matrix consisting of frequency of criterion variable:
#'   lenght(x) x ncat \cr
#'   R[,1] = Frequency of the zero-th category at x,
#'   R[,2] = Frequency of the first category at x, up to R[,ncat]
#' @param x A vector or matrix of unidimensional regressor variable
#' @param type Item type  "B", "B3", "G", "P", "PN", "N"
#' @param param A data frame containing the initial value of parameters,
#'  or  NULL
#' @param maxiter Max # of iterations
#' @param eps Criterion for convergence in llh
#' @param epsd Criterion for convergence in absolute value of parameter value
#' @param smallP Smallest value of probability
#' @param DinP = 0 to omit 1.7
#' @param minp1 Min value of the slope parameter for type != "N" items.
#' @param maxabsparam Max value of parameters
#' @param print = 0 to surpress output
#'
#' @details
#' Let R be n x ncat matrix, and x, n x 1.  \cr
#' This function maximizes the following log likelihood
#' w.r.t. item parameters. \cr
#' \eqn{ llh = \sum_{i=1}^{n} \sum_{j=1}^{ncat}  R[i,j] * log( P_j(x[i]) }
#' \cr
#' \cr
#' When \eqn{P_j} is 3PLM, \eqn{c} parameters will be transformed to
#' \eqn{y=logit(c)} and \eqn{llh} will be minimized w.r.t. \eqn{y}. \cr
#' Note that \cr
#' \eqn{d llh / d y = d llh / d c \times c(1-c)}, \cr
#' where \eqn{d logistic(y) / d y = c(1-c)}
#'
#'
#' @return A list of: \cr
#' param Parameter data frame \cr
#' llh Maximized log likelihood \cr
#' P Probability at x \cr
#' x input \cr
#' g gradient \cr
#' R input \cr
#' maxag Maximum of absolute gradient values
#'
#' @examples
#' # binary data
#' set.seed(1701)
#' theta=seq(-4,4,length.out=51)
#' resg=gendataIRT( 500, paramS1[1,], theta=theta, thdist="NORMAL", thd=NULL
#'                , thmean=0, thstd=1, compress=0 )
#' testdata=as.matrix( resg$U )
#'
#' # 3PLM
#' res=ordinal_reg( testdata, theta, type="B3", param=NULL, print=1 )
#' # 2PLM
#' res=ordinal_reg( testdata, theta, type="B", param=NULL, print=1 )
#' res=ordinal_reg( testdata, theta, type="P", param=NULL, print=1 )
#' # 2PLM with nonzero c parameter
#' initp=paramS1[2,]; initp$type="B"
#' res=ordinal_reg( testdata, theta, type="B", param=initp, print=1 )
#'
#'
#' # polytomous data
#' set.seed(1701)
#' theta=seq(-4,4,length.out=51)
#' resg=gendataIRT( 500, paramS1[3,], theta=theta, thdist="NORMAL", thd=NULL
#'                  , thmean=0, thstd=1, compress=0 )
#' testdata=as.matrix( resg$U )
#'
#' # graded response model
#' res=ordinal_reg( testdata, theta, type="G", param=NULL, print=1 )
#' # partial credit model
#' res=ordinal_reg( testdata, theta, type="P", param=NULL, print=1 )
#' res=ordinal_reg( testdata, theta, type="PN", param=NULL, print=1 )
#'
#'
#' @export
#'

ordinal_reg <- function( R, x, type="P", param=NULL
                         , maxiter=100, eps=1e-8, epsd=1e-5, smallP=0, DinP=1
                         , minp1=1e-2, maxabsparam=20
                         , print=1 ){
 # ordinal regression of (R,N) on scalar x
 # Shin-ichi Mayekawa
 # 20121022,23,24
 # valid parameter value: 121025
 # bugfix: 121104
 # pass Pj to dirf_p: 121105
 # DinP: 121110(London)
 # check too large values of abs(paramnum): 121115
 # bigfix: 121118
 # dc=logit(c): 121118
 # initial p3=c=0.3 for type="B3" items: 121118
 # nominal model: 121123
 # matSwp: 130204 ????
 # maxabsparam: 130204
 # minp1: 161115
 # solve -> matSwp: 20161116
 # avoid use of b_diag: 20161202
 # comment of loggit(c): 20180211
 #
 #  Args:
 #
 #    R       npoints x ncat   response data matrix
 #    x       npoints x 1      regressor matrix
 #    type    = B | B2 | B3 | G | P | PN | N   item type
 #    param  initial parameter data.frame: name, type, ncat, p1, p2, ...
 #
 #
 #  Values:
 #   list of param, P etc
 #
 # Needs:
 #  irf, dirf_p
 #
 # This program finds the parameter which maximizes the following likelihood
 #    llh = sum( R*log(P(x | param)) )
 #  where R is the npoints x ncat dummy expanded responce matrix,
 #  x is the npoints x 1 regressor vector, and
 #  P(x | param) is the npoints x ncat probability vector calculated as
 #  the IRT binary or polytomous item category response function at theta=x.
 #  The Fisher Scoring method is used to estimate the parameters.
 #

 # Print("**** Top of Ordinal Reg")


 # const
 if( is.vector(x) ){
  nobs=length(x); nq=1
 }
 else{
  nobs=nrow(x); nq=ncol(x)
 }
 ncat=ncol(R)  # This is the # of categories.  Tha max score is m=ncat-1.
 ncat1=ncat-1
 N=matrix(rowSums(R),,1)
 ntot=sum(N); nmiss=sum(N==0)
 nitems=1


 # error
 if( nq > 1 ){
  cat("\n\nerror:(ordinal_reg) *** Only one regressor is allowed. ******\n\n")
  return( NULL )
 }
 if( !substr(type,1,1) %in% c("B","G","P","N") ){
  cat("\n\nerror:(ordinal_reg) *** Incorrect type spec. **********\n\n")
  Print(type)
  return( NULL )
 }
 if( length(grep("^B[[:digit:]]*$",type)) > 0  && ncat >= 3 ){
  cat("\n\nerror:(ordinal_reg) *** Incorrect type spec. **********\n\n")
  Print(type, ncat)
  return( NULL )
 }

 # of parameters of the model:
 if( type == "B3" ){
  np=3; ncatx1=3
 }
 else if( type == "N" ){
  np=2*ncat1; ncatx1=np
 }
 else{
  np=ncat; ncatx1=ncat
 }

 # initial values
 if( is.null(param) ){
  # not given
  paramhead=as.data.frame(
   cbind( name="Q1" , type=type )
   , stringsAsFactors=0 )
  paramnum=matrix(0,1,np+1
                  , dimnames=list(NULL,c("ncat", paste("p",1:ncatx1,sep=""))))
  paramnum[1,1]=ncat
  paramnum[1,2]=1
  # increasing b-parameters with mean=0
  if( length(grep("^B[[:digit:]]*$",type)) == 0 ){
   paramnum[,3:(np+1)]=(1:(np-1))-np*(np-1)/2/(np-1)
  }
  else if( type == "B3" ){
   paramnum[,3]=0
   paramnum[,4]=0.3
  }
  else{
   paramnum=cbind(paramnum,0)
   colnames(paramnum)[4]="p3"
  }
  param=cbind(paramhead, paramnum)
 }
 else{
  # initial given
  # delete exessivecolumns
  locnotna=which( !apply( param, 2, function(x) all(is.na(x)) ) )
  param=param[,locnotna]
  param$type=type
 }

 if( is.null(param)  &&  type == "N" ){
  paramhead=as.data.frame(
   cbind( name="Q1" , type=type )
   , stringsAsFactors=0 )
  paramnum=matrix(0,1,np+1
                  , dimnames=list(NULL,c("ncat", paste("p",1:ncatx1,sep=""))))
  paramnum[1,1:(np+1)]=c(ncat, 1:ncat1, (1:(ncat-1))-ncat*(ncat-1)/2/(ncat-1))
  paramnum[1,1:(np+1)]=c(ncat, rep(1,ncat1), rep(0,ncat1))
  param=cbind(paramhead, paramnum)
  Print(param)
 }

 if( print > 0 ){
  cat("\nOrdinal Regression\n")
  cat(" # of observations =", nobs, ",  # of regressors =", nq,"\n")
  cat(" total # of responses =", ntot, " with missing =", nmiss,"\n")
  cat(" Model Type =", type, " with the # of categories =",ncat,"\n")
  cat(" # of parameters to be estimated =", np, "\n")
  cat(" max # of iterations =", maxiter," with eps =",eps
      ," and epsd =", epsd,"\n")
  cat(" minimum value of the the slope parameter for type != 'N' items ="
      , minp1, "\n")
  cat(" maximum value of the absolute value of the parameters ="
      , maxabsparam, "\n")
  cat(" print level =", print,"\n")
  cat(" inital parameter value \n")
  Print(param)
  if( print >= 3 ){
   Print(x,R,N)
  }
 }


 # initial llh
 # P is npoints x ncat irf
 P=irf( param, x, smallP=smallP, DinP=DinP, print=0 )$ICRF
 P[P<10e-324]=10e-324
 llh=sum(R*log(P))
 #Print("INITIAL", llh, R,N,P)

 # main iteration loop
 llhp=llh
 paramnump=param[,4:(4+np-1)]
 converged=0
 P1=NULL


 for( llll in 1:maxiter ){

  # log-Jacobian matrix by expanding the category first
  Jac=dirf_p( param, x, smallP=smallP, DinP=DinP, Pj=P1
              , log=1, cat.first=1, zero=1, print=0 )$Jack

  # simple d llh / d param
  #  gx=t( t(matrix(t(R),,1))%*%Jac )

  paramnum=param[,4:(4+np-1)]

  # expectation of  t(R[i,])%*%R[i,]
  # Dr=NULL
  # for( q in 1:nobs ){
  #  Drk=N[q]*( Diag(P[q,])-t(P[q,,drop=0])%*%P[q,,drop=0] )
  #  Dr=b_diag( Dr, Drk )
  # }

  Dr=matrix(0,nrow(Jac),nrow(Jac))
  for( q in 1:nobs ){
   Dr[((q-1)*ncat+1):(q*ncat),((q-1)*ncat+1):(q*ncat)]=
    N[q]*( Diag(P[q,])-t(P[q,,drop=0])%*%P[q,,drop=0] )
  }


  # type B3 or type PO or others
  if( type == "B3" ){

   # logit transformation of c=p3:    d logit(c) / d x = c*(1-c)
   locc=which(colnames(param)=="p3")
   c=param[1,locc]
   if( c == 0 ) c=1e-323
   Jac[,3]=Jac[,3]*c*(1-c)
   if( c > 0 ) dc=log(c/(1-c))
   else dc=-709
   paramnum[3]=dc

   # gradient etc
   delta=matrix(t((R-as.vector(N)*P)),,1)
   g=t(Jac)%*%delta
   H=t(Jac)%*%Dr%*%Jac

   # use sweep
   td=t( matSwp(H)%*%g )

   paramnew=param
   llhnew=llhp

   step=1; ok=0; llhnew=llhp
   for( lllll in 1:20 ){
    paramnumnew=paramnum + step*td
    if( paramnumnew[1] <=  0.1 ) paramnumnew[1]=0.1
    #    if( paramnumnew[3] < -709 ) paramnumnew[3]=-709
    if( paramnumnew[3] < -20 ) paramnumnew[3]=-20
    paramnumnew[3]=1/(1+exp(-paramnumnew[3]))
    paramnew[,4:(4+np-1)]=paramnumnew
    P=irf( paramnew, x, smallP=smallP, print=0 )$ICRF
    P[P<10e-324]=10e-324
    llhnew=sum( R*log(P) )
    if( llhnew >= llhp ){
     ok=1
     break
    }
    step=step/2
   } # end of lllll

  } # end of B3
  else if( type == "PO" ){

   comments('***')

  } # end of PO
  else{

   # Nominal, 2PLM, Graded, or GPCM w/o restrictions

   # gradient etc
   delta=matrix(t((R-as.vector(N)*P)),,1)
   g=t(Jac)%*%delta
   H=t(Jac)%*%Dr%*%Jac

   td=t( matSwp(H)%*%g )

   paramnew=param
   llhnew=llhp

   step=1; ok=0; llhnew=llhp
   for( lllll in 1:20 ){
    paramnumnew=paramnum + step*td
    paramnumnew[paramnumnew > maxabsparam]=maxabsparam
    paramnumnew[paramnumnew < -maxabsparam]=-maxabsparam
    if( type != "N"  &&  paramnumnew[1] <=  minp1 ) paramnumnew[1]=minp1
    paramnew[,4:(4+np-1)]=paramnumnew
    P=irf( paramnew, x, smallP=smallP, print=0 )$ICRF
    P[P<10e-324]=10e-324
    llhnew=sum( R*log(P) )
    if( llhnew >= llhp ){
     ok=1
     break
    }
    step=step/2
   } # end of lllll


  } # end of others


  if( ok ){
   # update parameters
   param=paramnew
   llh=llhnew
  }
  else{
   if( print >= 0 )
    cat("\n\n(ordinal_reg)**** Iteration did not converge. ****\n\n")
   llhimpr=0
   maxag=9999
   break
  }

  # convergence
  maxag=max(abs(g))
  maxadp=max(abs(paramnump-paramnumnew))
  llhimpr=(llh-llhp)/abs(llhp)
  if( print >= 2 ) Print(llll,llh,llhp,llhimpr, maxadp, maxag, digits=3)
  if( llhimpr <= eps  &&  maxadp <= epsd ){
   converged=1;
   break
  }

  # next iteration
  llhp=llh
  paramnump=paramnumnew
  P1=P[,2:ncat,drop=F]

 } # end of llll loop



 iter=llll
 if( print >= 1 && converged )
  cat("\nIteration converged with ", iter, " iterations.  eps= ",eps,"\n")

 if( print > 0 ){
  cat("\nOrdinal Regression\n")
  cat(" # of observations =", nobs, ",  # of regressors =", nq,"\n")
  cat(" total # of responses =", ntot, " with missing =", nmiss,"\n")
  cat(" Model Type =", type, "\n")
  cat(" # of parameters to be estimated =", np, "\n")
  cat(" # of iterations before termination =",iter)
  cat(" with eps =", eps, " and epsd =", epsd,"\n")
  cat(" final log likelihood value =", llh," with relative imp. ="
      , llhimpr," \n")
  cat(" minimum value of the the slope parameter for type != 'N' items ="
      , minp1, "\n")
  cat(" max value of the gradient =", maxag, "\n")
  cat(" print level =", print,"\n")
  cat("\n estimated parameter value  (DinP =",DinP,")\n")
  print(param)
  cat("\n final gradient value\n")
  print(t(g))
 }

 return( list(param=param, llh=llh, P=P, x=x, g=g, R=R, maxag=maxag) )

} # end of ordinal_reg






comments('


seed=34955
set.seed(seed)


paramtest=data.frame(
 name="Q1", type="B3", ncat=2,
 p1=1, p2=0, p3=0.3
 )

paramtest=data.frame(
 name="Q1", type="P", ncat=4,
 p1=.7, p2=-1, p3=0, p4=1
 )

paramtest=data.frame(
 name="Q1", type="P", ncat=5,
 p1=.7, p2=-1.5, p3=-1, p4=1, p5=1.5
)

paramtest=data.frame(
 name="Q1", type="P", ncat=6,
 p1=.7, p2=-1.5, p3=-1, p4=0, p5=1, p6=1.5
)

theta=seq(-4,4,length.out=51)

# temp=irf(paramtest,theta,plot=1,print=0)

resg=gendataIRT( 500, paramS1[1,], theta=theta, thdist="NORMAL", thd=NULL
               , thmean=0, thstd=1, compress=0 )
testdata=resg$U

# testdata=testdata[,c(2,1,5,4,3)]

# testdata=testdata[,c(2,1,3,4,5,6)]
# testdata=testdata[,c(2,1,5,4,3,6)]




res=ordinal_reg( testdata, theta, type="B", param=NULL
               , maxiter=30, eps=1e-9, print=2 )

temp=irf(res$param,theta,plot=1,print=0)



res=ordinal_reg( testdata, theta, type="G", param=NULL
                 , maxiter=30, eps=1e-9, print=2 )





paramtest1=data.frame(
 name="Q2", type="P", ncat=2,
 p1=1, p2=0, p3=NA, p4=NA, p5=NA, p6=NA
 , stringsAsFactors=0 )

paramtest2=data.frame(
 name="Q2", type="P", ncat=3,
 p1=1, p2=-1, p3=0, p4=NA, p5=NA, p6=NA
 , stringsAsFactors=0 )

paramtest3=data.frame(
 name="Q3", type="P", ncat=4,
 p1=1, p2=-1, p3=0, p4=2, p5=NA, p6=NA
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

paramtest0=convP2N( paramtest )


resg=gendataIRT( 500, paramtest0[3,], theta=theta, thdist="NORMAL", thd=NULL
                 , thmean=0, thstd=1, compress=0 )
testdata=resg$U




res=ordinal_reg( testdata, theta, type="N", param=NULL
                 , maxiter=30, eps=1e-9, epsd=1e-9, print=3 )

temp=irf( res$param, theta, print=0, plot=1 )




testdata_rev=testdata[,c(4,3,2,1)]

res=ordinal_reg( testdata_rev, theta, type="N", param=NULL
                 , maxiter=30, eps=1e-9, epsd=1e-9, print=2 )
temp=irf( res$param, theta, print=0, plot=1 )


res=ordinal_reg( testdata_rev, theta, type="P", param=NULL
                 , maxiter=30, eps=1e-9, epsd=1e-9, print=2 )
temp=irf( res$param, theta, print=0, plot=1 )



testdata_rev=testdata[,c(1,2,4,3)]

res=ordinal_reg( testdata_rev, theta, type="N", param=NULL
                 , maxiter=30, eps=1e-9, epsd=1e-9, print=2 )
temp=irf( res$param, theta, print=0, plot=1 )


res=ordinal_reg( testdata_rev, theta, type="P", param=NULL
                 , maxiter=30, eps=1e-9, epsd=1e-9, print=2 )
temp=irf( res$param, theta, print=0, plot=1 )

res=ordinal_reg( testdata_rev, theta, type="G", param=NULL
                 , maxiter=30, eps=1e-9, epsd=1e-9, print=2 )
temp=irf( res$param, theta, print=0, plot=1 )


')

set.seed(1701)
theta=seq(-4,4,length.out=51)
resg=gendataIRT( 500, paramS1[1,], theta=theta, thdist="NORMAL", thd=NULL
                 , thmean=0, thstd=1, compress=0 )
testdata=as.matrix( resg$U )
res=ordinal_reg( testdata, theta, type="G", param=NULL, print=2 )

comments('
set.seed(1701)
theta=seq(-4,4,length.out=51)
resg=gendataIRT( 500, paramS1[1,], theta=theta, thdist="NORMAL", thd=NULL
               , thmean=0, thstd=1, compress=0 )
testdata=as.matrix( resg$U )

res=ordinal_reg( testdata, theta, type="B3", param=NULL, print=2 )
res=ordinal_reg( testdata, theta, type="B", param=NULL, print=2 )
res=ordinal_reg( testdata, theta, type="P", param=NULL, print=2 )


set.seed(1701)
theta=seq(-4,4,length.out=51)
resg=gendataIRT( 500, paramS1[3,], theta=theta, thdist="NORMAL", thd=NULL
                 , thmean=0, thstd=1, compress=0 )
testdata=as.matrix( resg$U )

res=ordinal_reg( testdata, theta, type="G", param=NULL, print=1 )
res=ordinal_reg( testdata, theta, type="P", param=NULL, print=1 )
res=ordinal_reg( testdata, theta, type="PN", param=NULL, print=1 )

')
