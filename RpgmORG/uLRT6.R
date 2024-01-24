#' Item Parameter Estimation of Unidimensional LRT
#'
#' Marginal likelihood will be maximized by the EM algorithm
#' using ordinal_reg funcion and smn function in the M-step.
#'
#' @param Uc   n x nitems+1 or +2    compressed item resopnse data frame
#'          in BILOG-MG's expanded format. (idvar and groupvar)
#' @param U n x sum(ncat)+1 or +2 uncompressed irem response data frame
#' @param   groupvar name of the grouping varible contained as a column of Uc
#'  or NULL.
#' @param   idvar    name of the id varible contained as a column of Uc
#'  or NULL
#' @param   ncat     nitems x 1 #  of categories for each item or fNULL.
#' @param   type     nitems x 1 vector of item types consisting of:
#'          "Bn" | "G" | "PN" | "P" \cr
#'          If NULL, P will be used. \cr
#'          *** Currently, polytomous items are not allowed.
#' @param itemname nitems x 1 vector of item names or NULL.
#' @param nclass # of classes
#' @param V  npoitns x sum(ncat)  initial value for item param\cr
#'         Only the subset of the items can be given.\cr
#'         param$name and the colname(Uc) will be used to match the items.
#' @param   fixeditems list of items whose item parameters are to be fixed\cr
#'   to the values given in the param data frame.\cr
#'  Item numbers or item names\cr
#'  *** Currently, this option is not available.
#' @param  baseform base form number whose theta distribution is
#' fixed at the values given in msn[baseform,] \cr
#'  If there are fixed items, the values of the item paramters
#'  given in the param data frame must be on the baseform scale.\cr
#' @param   rho initial value of probability vector of each latent class.
#' @param   estrho    = 1  to estimate multigroup theta means
#' @param   alpha = vector of Dirichle parameters of length nclass \cr
#' or a number A to set alpha=rep(A,nclass) \cr
#' or a negative number -A to set alpha=A*normal pdf
#' @param   vmin minimum value of V
#' @param   vmax maximum value of V
#' @param   maxiter  max # of iterations
#' @param   eps      eps for the relative improvement of lmlh
#' @param   epsd     eps for the max. abs. diff. of msn\cr
#'     Seems this is important.
#' @param   maxiter2 max # of iterations for ordinal_reg and smn
#' when llll <= nstrict.
#' @param   eps2     eps for the relative improvement of
#' llh in ordinal_reg and smn.
#' @param   nstrict  maxiter2 will be reduced to maxiter22
#' after nstrict iterations.
#' @param   maxiter22  max # of iterations for ordinal_reg
#' and smn when llll > nstrict.
#' @param minp1 Minimum value of p1 parameter for type != "N" items.
#' @param maxabsparam Maximum absolute value of parameters.
#'
#' @param SQUAREM = 3 :   See the help of iSQUAREM in lazy.accel package.
#' @param nSQUAREM when to star iSQUAREM
#' @param minalpha = -999 :   See the help of iSQUAREM in lazy.accel package.
#' @param  maxalpha = -1 :   See the help of iSQUAREM in lazy.accel package.
#' @param  always = 0 :   See the help of iSQUAREM in lazy.accel package.
#' @param  reset1 = 1 :   See the help of iSQUAREM in lazy.accel package.
#' @param  reset2 = 2 :   See the help of iSQUAREM in lazy.accel package.
#'
#' @param print = 1 to print the result
#' @param plot = 1 to plot the estimated theta distributions
#' @param smallP = Minimum value of P
#' @param debug = 1 to print intermediate result
#'
#' @references
#' Varadhan, R. and Roland, C.(2007) Simple and Globally Convergent Methods
#' for Accelerating the Convergence of Any EM Algorithm.
#' Scandinavian Journal of Statistics, Vol. 35: 335-353. \cr
#' Shojima, K. (2008a) Neural test theory.
#' In K. Shigemasu, A. Okada, T. Imaizumi, & T. Hoshino (Eds.)
#' New trends in psychometrics. Universal Academy Press, Inc
#'
#' @return A list of:
#' V Estimated Item Parameters in a data frame \cr
#' rho Class distribution \cr
#' converged = 1 if converged, = 0 otherwise. \cr
#' lmlh Log Marginal LIkelihood maximized \cr
#' iter_hist A list of iteration history \cr
#' H Posterior distribution of theta given U.  \cr
#' id Id variable \cr
#'
#' @details
#' *** Currently, polytomous items are not allowed. ****** \cr
#' *** Currently, multi-group analysis is not available. \cr
#' \cr
#'  Uc is the compressed data if it is n x nitems. \cr
#'  U is the uncompressed data if it is n x sum(ncat). \cr
#'  \cr
#'  When U, instead of Uc, is given: \cr
#'  Item names come from itemname or paste("Q",1:nitems,sep=""). \cr
#'  ncat must be given. \cr
#'
#'  When Uc is given: \cr
#'   Item name comes from colnames(Uc) or itemname. \cr
#'   ncat can be calculated from Uc.\cr
#'
#'  \cr
#'  Note on the iSQUAREM: \cr
#'  Try always=1 with maxalpha=-1 or less first. \cr
#'  Changing to reset1=1 and reset2=2 or increasing nSQUAREM may help. \cr
#'  If it seems not working, use always=0 with maxalpha=1 or less. \cr
#'  Changing to reset1=1 and reset2=2 may help. \cr
#'  If all of the above fail, be patient and use SQUAREM=0.
#'  \cr
#'
#'
#' @examples
#' #
#' #### In the following examples, maxiter is set to 20 which is
#' #### not large enough to obtain convergence.
#' ####
#' #
#' #
#' #
#' set.seed(1701)
#'
#' param=paramB1[c(1:3,7:9,13:15),]
#' thmin=-2; thmax=2; npoint=5
#' N=1000
#'
#' # discrete theta
#' # theta0=seq(thmin,thmax,length=npoint)
#' theta0=c(-2, -1, 0, 2, 3)
#' theta=unlist(lapply( theta0, rep, round(N/npoint) ))
#' res2 <- gendataIRT( 1, paramB1, theta=theta, compress=1 )
#' Uc=as.data.frame(res2$U)
#' ncat=res2$ncat
#' type=res2$type
#'
#' nclass=5
#' res1 <- uLRT( Uc, nclass=nclass, estrho=1, monotone=1, alpha=20
#' , maxiter=20, plot=1, print=1 )
#'
#' # normal theta
#' Uc2 <- gendataIRT( 1, paramB1, npoints=N, thdist="rnorm", compress=1 )$U
#' Uc2 <- as.data.frame(Uc2)
#'
#' nclass=5
#' res1 <- uLRT( Uc2, nclass=nclass, estrho=1, monotone=1
#' , maxiter=20, plot=1, print=1 )
#'
#' @export
#'

uLRT <- function( Uc, U=NULL, groupvar=NULL, idvar=NULL
                  , ncat=NULL, type=NULL, itemname=NULL, nclass=5
                  , V=NULL, fixeditems=NULL, baseform=1, monotone=0
                  , rho=NULL, estrho=0
                  , alpha=rep(1,nclass)
                  , vmin=0.0001, vmax=1-vmin
                  , maxiter=200, eps=1e-7, epsd=1e-6
                  , maxiter2=1, eps2=0, nstrict=9, maxiter22=5
                  , minp1=1e-1, maxabsparam=20
                  , SQUAREM=3, nSQUAREM=1, minalpha=-999, maxalpha=-1
                  , always=1, reset1=0, reset2=1
                  , print=2, plot=0, smallP=0, debug=0 ){
 # Nonparametric IRT (NTT or LRT)
 # Shin-ichi Mayekawa
 # uIRT modifed: 20170619,20
 # LRT renamed as uLRT: 20170620
 # monotone: 20170621
 # prior alpha: 20170622
 # get_monotone -> mbinreg: 20170803
 # nIRT renamed as uLRT: 20170803
 # mbinreg by pava: 20170805
 #
 #
 #
 # Args:
 #
 #  Uc       n x nitems+1 or +2    compressed item resopnse data frame
 #           in BILOG-MG's expanded format
 #  groupvar name of the grouping varible contained as a column of Uc
 #  idvar    name of the id varible contained as a column of Uc
 #
 #  ncat     nitems x 1  # of categories for each item
 #  type     nitems x 1  item type as   "Bn" | "G" | "PN" | "P"
 #           If NULL, P will be used.
 #
 #  param    nitems x max(ncat)  initial value data frame for item param
 #           Only the subset of the items can be given.
 #           param$name and the colname(Uc) will be used to match the items.
 #  msn      nG x 3  initial value matrix of (mean, std, n) for each group
 #
 #  fixeditems list of items whose item parameters are to be fixed
 #           to the values given in the param data frame.
 #           item numbers or item names
 #
 #  baseform base form number whose theta distribution is fixed at the values
 #           given in msn[baseform,]
 #
 #           If there are fixed items, the values of the item paramters
 #           given in the param data frame must be on the baseform scale.
 #
 #  estrho    =1  to estimate multigroup theta means
 #
 #  theta    nclass x 1   discrete theta points
 #  nclass, # of discrete theta points between (thmin, thmax)
 #
 #  maxiter  max # of iterations
 #  eps      eps for the relative improvement of lmlh
 #  epsd     eps for the max. abs. diff. of msn
 #           Seems this is important.
 #  maxiter2 max # of iterations for ordinal_reg and smn when llll <= nstrict
 #  maxiter22  max # of iterations for ordinal_reg and smn when llll > nstrict
 #  eps2     eps for the relative improvement of llh in ordinal_reg and smn
 #  nstrict  maxiter2 will be reduced to maxiter22 after nstrict iterations
 #
 #  SQUAREM  = 1 to emply SQUAREM to accelerate
 #  minalpha minimum value of alpha parameter:  default is -6.
 #           If minalpha > 0,  alpha will be fixed to minalpha.
 #  nSQUAREM iteration number from which SQUAREM update begins
 #
 #
 #
 #  Values:
 #   list of
 #     V       estimated parameter data frame
 #     msn     mu, sigma, n  of theta distriution for each group
 #     (theta, thd) theta distribution
 #     H       n x nclass   posterior of theta given data
 #     (EAP, poststd)   n x 2    estimated theta and its std
 #
 #
 #  Needs:
 #   dummy_expand, irf, dirf_p, ordinal_reg, smn
 #
 #
 #
 #  Comments:
 #   When theta distribution is estimated, mu and sigma for the baseform group
 #   is fixed to (0,1).
 #
 #   When some of the item parameters are fixed,
 #   theta distribution of the baseform group will also be estimated
 #   if estrho == 1.
 #
 #
 #   Method SqS3 is implemented for SQUAREM.
 #     alpha > minalpha is enforced.
 #
 #
 #
 #
 #
 #
 #   If SQUAREM == 3 and minalpha == -6,  it seems maxiter2=1 seems to be OK.
 #   However, for the first few iterations, use maxiter2=10 or so.
 #   maxiter2=1 may cause problem when n or nitems are large.
 #
 #




 getlmlh <- function( U, group, V, rho, alpha1 ){
  # calculation of log marginal likelihood and the posterior expectation
  # Shin-ichi Mayekawa
  # 20121024
  # output H: 121112
  # limits of exp or log: 121117,18,19
  # return R: 121217
  # comments: 20161127
  # rowmin and 720, not 745: 20161214
  # 20170619
  #
  # Args:
  #   U      n x sum(ncat)  expanded data
  #          U[i,fromP[j]:toP[j]]=0 if Uc[i,j] is missing
  #   group  n x 1
  #   V      npoint x sum(ncat)   zeroth category included
  #   rho    npoint x nG  class prob for each group (sums to unity)
  #
  # Values:
  #   lmlh   log marginal likelihood
  #   logP   nclass x sum(ncat) log of IRF
  #   H      n x nclass  posterior expectation of class indicator given obs
  #   R      nclass x sum(ncat) x nG array of post exp of U at class
  #
  #  Note
  #     min arg to 1/1e-x is 308
  #     min arg to log is 1e-308:
  #     max arg to exp is 709.782:
  #

  nG=length(unique(group))
  logP=log(V)
  H=matrix(0,nrow(U),nrow(V))
  lmlh=0
  R=array(NA,c(ncol(H),ncol(U),nG))
  for( g in 1:nG ){
   locg=which(group==g)
   # eULP=exp(U[locg,,drop=0]%*%t(logP)) * matrix(1,length(locg))%*%t(rho[,g])
   ULP=U[locg,,drop=0]%*%t(logP) + matrix(1,length(locg))%*%t(log(rho[,g]))

   rowmin=apply(ULP,1,min)
   ULP=ULP-rowmin  - 600

   ULP[ULP > 709]=709
   eULP=exp(ULP)
   rseULP=rowSums(eULP)
   rseULP[rseULP < 1e-307]=1e-307
   rseULP[rseULP > 10e307]=10e307
   lmlh=lmlh+sum( log( rseULP ) + rowmin + 600 ) + sum(alpha1*log(rho[,g]))
   H[locg,]=eULP/rseULP
   R[,,g]=t(H[locg,])%*%U[locg,,drop=0]
  }
  # Print(rseULP,H,fmt="8.5"); # stop()
  rm(ULP,eULP,rseULP)
  return( list(lmlh=lmlh, logP=logP, H=H, R=R) )

 } # end of getlmlh



 fit <- function( v ){
  # returns minus log marginal likelihood etc for iSQUAREM
  # Shin-ichi Mayekawa
  # 20160102,03
  # 20170620
  #

  # Everything except for "newp" comes from the environment
  # in which this function is defined. i.e., uLRT.

  # reshape parameters
  V=matrix(v,nclass)

  # evaluate lmlh
  temp=getlmlh( U, group, V, rho, alpha1 )
  names(temp)=c("critval",names(temp)[-1])
  temp[[1]]=-temp[[1]]
  return( temp )
 } # end of fit


 const <- function( v ){
  # enforces constraints on the parameters for iSQUAREM
  # Shin-ichi Mayekawa
  # 20160102,03
  # 20170620

  # Everything except for "newp" comes from the environment
  # in which this function is defined. i.e., uLRT.

  # reshape parameters
  V=matrix(v,nclass)
  V[V<vmin]=vmin
  V[V>vmax]=vmax
  v=c(V)

  return( v )
 } # end of const




mbinreg <- function( f, n, x=1:length(f), method=4, print=0, plot=0 ){
 # find monotone pi which maximize the binominal likelihood
 # Shin-ichi Mayekawa
 # 20170621
 # initial by isoreg: 20170803
 # Fisher Scoring Method: 20170803,04ShinbashiSta
 # gpava: 20170804
 # bugfix: 20170808
 #

 # method
 mm=substr(toupper(method),1,1)
 if( mm %in% c("P","C") ) method=3
 else if( mm  %in% c("S","D") ) method=4
 else if( method <= 3 ) method=3
 else method=4

 # vectorize
 x=c(x); f=c(f); n=c(n)
 nobs=length(x)
 if( length(n) == 1 ) n=rep(n,nobs)
 if( nobs != length(f) || nobs != length(n) ){
  cat("\n\n error1:(mbinreg) x and  ( f, n ) must be the same size. \n\n")
  return()
 }


 # unrestricted ml = observed prob
 p0=f/n


 # sort (x,fn,id) by x
 fn=cbind(f,n)
 id=1:nobs
 od1=order(x)
 x=x[od1]; fn=fn[od1,]; id=id[od1]; p0=p0[od1]

 # tie info
 tie0=matrix(table(x))
 tie=matrix(unlist(mapply( rep,tie0,tie0 )),,1)

 if( print >= 1 ){
  cat("\n\nMonotone Binary Regression of (f, n) on x. \n\n")
  cat("  # of observations =", nobs,"\n")
 }


 # Here (x,fn,id,tie) is sorted by x

 phat=fn[,1]/fn[,2]

 # discrete or continuous
 if( method == 4 ){

  # discrete processing
  # replace f/n by sum(f)/sum(n) corresponding to tied x

  phat=unlist(
   by( fn, x, function(z) rep(sum(z[,1])/sum(z[,2]),length(z[,1]) ) ) )

   # if method = 2, above phat is the final result.

 } # end of 4
 else if( method == 3 ){

  # continuous ordinal or nominal
  # sort fn corresponding to tied x by fn

#  comments('
#  i=1
#  while( i <= nobs ){
#   to=i+tie[i]-1
#   if( to > i ){
#    od=order(phat[i:to])
#    fn[i:to,]=fn[i:to,][od,]
#    phat[i:to]=fn[i:to,1]/fn[i:to,2]
#    id[i:to]=id[i:to][od]
#   }
#   i=to+1
#  }
#  ')

  # seems do.call is faster
  # all2=Reduce( rbind, by( cbind(fn,id,phat,x), x
  #                   , function(z){ od=order(z[,4]); z=z[od,]; return(z) } ) )
  all2=do.call( rbind, by( cbind(fn,id,phat,x), x
              , function(z){ od=order(z[,4]); z=z[od,]; return(z) } ) )
  f=all2[,1]; n=all2[,2]; id=all2[,3]; phat=all2[,4]; x=all2[,5]
  rm(all2)
  fn=cbind(f,n)

 } # end of 3


 # enforce monotonicity: The Core Part

 # replace non-monotonic part by the average
 for( i in 2:nobs ){
  for( j in 1:(i-1) ){
   if( phat[i-j] > phat[i] ){
    phat[(i-j):i]=sum(fn[(i-j):i,1])/sum(fn[(i-j):i,2])
   }
   else break      #### new 111130
  }
 }


 # sort everything back to the original order
 xyhs=cbind(id, x,fn,phat,tie);
 colnames(xyhs)=c("ids","xs","fs","ns", "yhats","ties")
 od=order(id)
 x=x[od]; fn=fn[od,]; phat=phat[od]; tie=tie[od]; id=id[od]

 llh=sum(fn[,1]*log(phat) + (fn[,2]-fn[,1])*log(1-phat))


 if( print >= 1 ){
  cat("\n x, (f, n) and phat in original and sorted order\n")
  Print( method, llh )
  xyh=cbind(id,x,fn,phat,tie);
  colnames(xyh)=c("id","x","f","n","phat","tie")
  cat("Result in the Original Order\n")
  Print( xyh, fmt="9.5" )
  cat("Result in the Sorted Order by x\n")
  Print( xyhs, fmt="9.5" )
 }

 # In order to plot the method=3 result nicely,
 # we must sort (x,fn,phat) by x and phat.
 if( plot > 0 ){
  title=paste("Monotone Binary Regression")
  title2=paste("method = ", as.character(method), sep=" ")
  title3=paste(":  LLH = ",llh)
  title=paste( title, title2, title3, sep=" ")
  xyh=cbind(x,fn,phat)
  od=order(x,phat)
  xyh=xyh[od,,drop=F]
  p0=xyh[,2]/xyh[,3]
  Print(x,p0)
  ylim=c(0,1)
  plot( xyh[,1],p0, xlab="", ylab="", ylim=ylim )
  par(new=T)
  plot( xyh[,1],xyh[,4],type="s", xlab="", ylab="", ylim=ylim )
  par(new=T)
  plot( xyh[,1],xyh[,4],type="p", pch="M", xlab="X", ylab="p-obs and p-hat"
        , ylim=ylim, main=title )
 }

 return( phat )

} # end of mbinreg






 #########################################################################
 # Main Starts Here.
 #########################################################################


 #
 # Process either Uc or U data frame.
 #

 if( is.null(U) ){

  # U is NOT given.
  Uin=0

  # const
  n=nrow(Uc)

  # chech if dataset exists
  dfname=deparse(substitute(Uc))
  dfname=unlist( strsplit(dfname,"[",fixed=1) )[1]
  if( !exists(dfname) ){
   cat("\n\nerror1(uLRT) **** data frame ", dfname
       , " does not exist. ****\n\n")
   return(NULL)
  }
  else if( !is.data.frame(Uc) ){
   cat("\n\nerror1(uLRT) **** ", dfname, " is not a data frame. ****\n\n")
   return(NULL)
  }

  # grouping variable
  if( is.null(groupvar) ){
   groupvar="none"
   group=rep(1,nrow(Uc))
   baseform=1
  }
  else{
   locgvar=which(colnames(Uc)==groupvar)
   if( length(locgvar) == 0 ){
    cat("\n\n error1(uLRT): group variable does not exist.\n\n\n")
    Print(groupvar)
    return( NULL )
   }
   # make group an integer variable containing 1,2,...,nG
   group=Uc[,locgvar]
   Uc=Uc[,-locgvar,drop=0]
   ugroup=unique(group)
   nG=length(unique(group))
   names(ugroup)=paste("group",1:nG,sep="")
   group0=group
   for( g in 1:nG ) group0[group0 == ugroup[g]]=g
   group=group0
   rm(group0)
  }
  nG=length(unique(group))
  nobs=table(group)
  if( nG == 1 ) baseform=1

  # id variable
  if( is.null(idvar) ){
   idvar="none"
   id=paste("s",1:nrow(Uc),sep="")
  }
  else{
   locidvar=which(colnames(Uc)==idvar)
   if( length(locidvar) == 0 ){
    cat("\n\n error1(uLRT): ID variable does not exist.\n\n\n")
    Print(idvar)
    return( NULL )
   }
   #
   id=Uc[,locidvar]
   Uc=Uc[,-locidvar,drop=0]
  }

  # base form
  if( baseform <= 0  |  nG < baseform ){
   cat("\n\n error1(uLRT): baseform must be in [1,",nG, "]\n\n\n")
   return( NULL )
  }

  # const
  nitems=ncol(Uc)
  if( is.null(itemname) ) itemname=colnames(Uc)

  # missing data count
  pmiss=sum(is.na(Uc))/(n*nitems)

  # dummy expand the compressed data
  #   U[i,fromP[j]:toP[j]]=0 if Uc[i,j] is missing
  temp=dummy_expand( Uc )
  rm(Uc)
  U=temp$U
  fromP=temp$fromP; toP=temp$toP; ncat=temp$ncat
  rm(temp)

 } # end of U not given
 else{

  # U is given.
  Uin=1

  # const
  n=nrow(U)
  nitems=length(ncat)

  # check if ncat is given
  if( is.null(ncat) ){
   cat("\n\nerror1(uLRT): *** When U is given, ncat must be given. *** \n")
   return(NULL)
  }
  if( is.null(itemname) ){
   itemname=paste("Q",1:length(ncat),sep="")
  }

  # chech if dataset exists
  dfname=deparse(substitute(U))
  dfname=unlist( strsplit(dfname,"[",fixed=1) )[1]
  if( !exists(dfname) ){
   cat("\n\nerror1(uLRT) **** data frame ", dfname
       , " does not exist. ****\n\n")
   return(NULL)
  }
  else if( !is.data.frame(U) ){
   cat("\n\nerror1(uLRT) **** ", dfname, " is not a data frame. ****\n\n")
   return(NULL)
  }

  # grouping variable
  if( is.null(groupvar) ){
   groupvar="none"
   group=rep(1,nrow(U))
   baseform=1
  }
  else{
   locgvar=which(colnames(c)==groupvar)
   if( length(locgvar) == 0 ){
    cat("\n\n error1(uLRT): group variable does not exist.\n\n\n")
    Print(groupvar)
    return( NULL )
   }
   # make group an integer variable containing 1,2,...,nG
   group=c[,locgvar]
   c=c[,-locgvar,drop=0]
   ugroup=unique(group)
   nG=length(unique(group))
   names(ugroup)=paste("group",1:nG,sep="")
   group0=group
   for( g in 1:nG ) group0[group0 == ugroup[g]]=g
   group=group0
   rm(group0)
  }
  nG=length(unique(group))
  nobs=table(group)
  if( nG == 1 ) baseform=1

  # id variable
  if( is.null(idvar) ){
   idvar="none"
   id=paste("s",1:nrow(U),sep="")
  }
  else{
   locidvar=which(colnames(U)==idvar)
   if( length(locidvar) == 0 ){
    cat("\n\n error1(uLRT): ID variable does not exist.\n\n\n")
    Print(idvar)
    return( NULL )
   }
   #
   id=U[,locidvar]
   U=U[,-locidvar,drop=0]
  }

  # base form
  if( baseform <= 0  |  nG < baseform ){
   cat("\n\n error1(uLRT): baseform must be in [1,",nG, "]\n\n\n")
   return( NULL )
  }

  # replace missing by 0
  U[is.na(U)]=0

  # location index
  fromto=get_range( ncat )
  fromP=fromto[,1]; toP=fromto[,2]
  U=as.matrix(U)

  # missing data count
  Ucc=mapply( function(f,t){rowSums(U[,f:t])}, fromP, toP )

  pmiss=sum(Ucc == 0)/(n*sum(ncat))
  rm(Ucc)


 } # end of U given


 ########################################################################
 ####################### Here, U is ready. ##############################
 ########################################################################


 # item type
 if( is.null(type) ){
  type=rep("B",nitems)
  type[which(ncat>2)]="P"
 }
 else{
  type=type[1:nitems]
  if( !( (length(grep("^B",type[ncat == 2])) == sum(ncat == 2))
         &&  (length(grep("^B",type[ncat == 3])) == 0)   ) ){
   cat("\n\n error1(uLRT): Item type does not match the # of categories.\n\n")
   Print(param)
   return(NULL)
  }
 }

 # of parameters of the model:
 np=ncat
 np[which(type == "B3")]=3
 np[which(type == "N")]=2*(ncat[which(type == "N")]-1)
 maxnp=max(np)
 if( length(grep("^B[2]*$",type)) == nitems ) maxnp=maxnp+1

 if( debug > 0 )
  Print(  nitems, nG, nclass, "/", type, ncat, np )

 # initial parameter data frame
 paramORG=NULL
 # create an empty data frame
 paramhead=as.data.frame( cbind( name=itemname , type=type )
                          , stringsAsFactors=0 )
 paramnum=matrix(NA,nitems,maxnp+1
             , dimnames=list(itemname,c("ncat", paste("p",1:maxnp,sep=""))))
 param=cbind(paramhead, paramnum)
 rownames(param)=paste("item",1:nitems,sep="")
 if( debug > 0 ) Print(" initial param with NA", param)


 # insert the defined part of paramORG into param
 initemlist=NULL
 if( !is.null(paramORG) ){
  for( j in 1:nrow(paramORG) ){
   paramj=paramORG[j,]
   namej=paramj$name
   typej=paramj$type
   npj=paramj$ncat
   if( length(grep("^B[[:digit:]]*$",typej)) > 0 ) npj=3
   else if( typej == "N" ) npj=2*(npj-1)
   if( !any(is.na(paramj[3:(3+npj)])) ){
    locj=which( itemname == namej )
    if( length( locj ) > 0 ){
     param[locj,1:(3+npj)]=paramj[1:(3+npj)]
     initemlist=c(initemlist,locj)
    }
   }
  }
 }
 rm(paramORG)
 if( debug > 0 )
  Print(  "after inserting paramORG", type, ncat, np, param )


 # default initial parameter values
 for( j in 1:nitems ){
  if( is.na(param$ncat[j]) ){
   # initial parameter values or fixed values are not given for item j.
   param$ncat[j]=ncat[j]; param$p1[j]=1
   if( type[j] %in% c("G", "PN", "P") ){
    # equally spaced b parameters
    param[j,5:(np[j]+3)]=(1:(np[j]-1))-np[j]*(np[j]-1)/2/(np[j]-1)
   }
   else if( type[j]  %in% c("B", "B2") ){
    param$p2[j]=0; param$p3[j]=0
   }
   else if( type[j] == "B3" ){
    param$p2[j]=0; param$p3[j]=0.3   # do not use 0 for 3PLM.
   }
   else if( type[j] == "N" ){
    #Print("****",j)
    param[j,4:(4+ncat[j]-1-1)]=1          # a paramters
    param[j,(4+ncat[j]-1):(4+np[j]-1)]=0  # b paramters
   }
  } # end of initial not given for item j
 }

 # item paramters fixed or not
 if( length(fixeditems) > 0 && length(initemlist) > 0 ){
  if( is.character(fixeditems) )
   fixeditems=which( itemname %in% fixeditems )
  fixeditemlist=intersect(fixeditems,initemlist)
  if( length(fixeditemlist) > 0 ){
   FixedItems=1
  }
 }
 else{
  FixedItems=0; fixeditemlist=NULL
 }
 nonfixeditemlist=setdiff(1:nitems,fixeditemlist)

 # location of non-fixed params excluding the c-parameters of B2 items.
 pp=param[nonfixeditemlist,4:(4+maxnp-1)]
 ppt=param$type[nonfixeditemlist]
 pp$p3[which(ppt == "B" | ppt == "B2")]=NA
 locSQEM=which(!is.na(pp))
 # if( debug ) Print(pp,locSQEM)
 rm(pp,ppt)
 # param[nonfixeditemlist,4:(4+maxnp-1)][locSQEM]


 # c-parameters of B2 items : all items, not restrocted to nonfixeditems
 locB2=which(param$type == "B" | param$type == "B2")
 cB2=param$p3[locB2]
 locB3=which(param$type == "B3")

 # # of paramters to be estimated
 nparam1=nclass*(sum(ncat)-nitems)-sum(np[fixeditemlist])
 nparam2=0
 if( estrho == 1 ) nparam2=nparam2+nG*(nclass-1)
 if( FixedItems ){
  if( estrho == 1 ) nparam2=nparam2+1  # for the baseform
 }
 nparam=nparam1+nparam2



 # class distribution
 if( length(alpha) == 1 ){
  if( alpha >= 0 ){
   alpha=rep(alpha,nclass)
  }
  else{
   dd=dnorm(seq(-3,3,length=nclass))
   alpha=-alpha*dd/sum(dd)
  }
 }
 alpha[alpha<1]=1

 alpha1=alpha-1
 if( is.null(rho) ){
  rho=matrix(0,nclass,nG)
  for( g in 1:nG ){
   temp=alpha1+1
   rho[,g]=temp/sum(temp)
  }
 }
 else if( is.vector(rho) ){
  rho=matrix(rho,,1)
 }
 classname=paste("class",1:nclass,sep="")
 rownames(rho)=classname


 # what to estimate
 estms=matrix(0,nG,1)
 if( estrho == 1  )
  for( g in 1:nG )
   if( FixedItems == 1  ||  g != baseform )
    estms[g,]=estrho
 locms=which(estms == 1)


 if( print > 0 ){
  cat("\n\n\nuLRT: parameter estimation of nonparametric IRT models\n")
  cat("\n Input data frame name = ", dfname, "  (compressed=",1-Uin,")\n"
      , sep="")
  cat("\n id varaible =", idvar, ",  group varaible =", groupvar,"\n")
  cat(" # of subjects =", n, ",  # of items =", nitems,"\n")
  cat(" # of groups =", nG, "\n")
  if( nG > 1 ){
   cat("  ");  print(ugroup)
  }
  cat(" percentage of missing data =", pmiss,"\n")
  cat(" # of latent classes =", nclass,"\n")
  cat(" estimation of class probabilities\n")
  cat("   estrho =", estrho, "\n")
  cat(" prior distribution of rho: Dirichle with ", alpha, "\n")
  cat("   baseform =", baseform, "\n")
  if( FixedItems ){
   cat(" The item paramters of the following items are fixed. \n")
   print(fixeditemlist)
  }
  cat(" total # of parameters to be estimated =", nparam
      ," (",nparam1," + ",nparam2 ,") \n")
  cat(" min and max values of estimated p1 parameter =(", vmin,",", vmax,")\n")
  cat(" max abs value of estimated parameter =", maxabsparam,"\n")
  cat(" max # of iterations =", maxiter," with eps =",eps
      , " and epsd =", epsd,"\n")
  cat(" max # of iterations in m-steps =", maxiter2,"/",maxiter22
      ," with eps =",eps2, "and nstrict =", nstrict,"\n")
  if( SQUAREM){
   cat("\n SQUAREM =", SQUAREM," will be used after", nSQUAREM
       , "iterations with: \n")
   cat( " minalpha = ", minalpha, ",  maxalpha = ",  maxalpha
        , ",  always = ", always, ",  reset1 = ", reset1
        , ",  reset2 =", reset2, "\n\n")
  }
  cat(" print level =", print,"\n")
  if( length(initemlist) > 0 ){
   cat(" initial item parameter values of the following items are given.\n")
   print(initemlist)
  }
  if( print >= 2 ){
   cat(" inital parameter value \n")
   Print(param)
   if( print >= 3 ){
    Print(rho,digits=1)
   }
  }
 }


 #initial
 V=irf( param, npoints=nclass, print=0 )$ICRF
 rownames(V)=classname
 cnv=colnames(V)
 temp=getlmlh( U, group, V, rho, alpha1 )
 lmlhp=temp$lmlh
 lmlh=lmlhp
 logP=temp$logP
 H=temp$H
 R=temp$R
 rm(temp)
 paramnum=as.matrix(param[nonfixeditemlist,4:(4+maxnp-1)])
 paramnump=paramnum

 # design matrix for monotone
 A=matrix(0,nclass,nclass)
 for( i in 1:nclass ){
  A[i,1:i]=1
 }


 if( print > 0 ) Print("initial", lmlhp)
 converged=0


 ###################################################################(1)
 # store param
 # initialize inline SQUAREM
 rmv(iSQUAREM)
 pp=c( V )
 iSQUAREM=generate_iSQUAREM( pp, print=0 )
 param0=V
 rho0=rho
 ###################################################################(1)


 iter_hist=matrix(c(0,lmlhp,0,0),1,4)
 colnames(iter_hist)=c("llll","lmlh","maxadp","maxadrho")
 Vp=V
 rhop=rho
 conv=0

 # main iteration
 for( llll in 1:maxiter ){


  # e-step
  # expected value of U at class
  # Rall is nclass x sum(ncat)
  Rall=apply(R,c(1,2),sum)

  if( llll <= nstrict ) maxiter21=maxiter2
  else maxiter21=maxiter22


  if( debug > 0 ) Print("*top of m-step for items")
  # m-step for item parameters
  for( j in 1:nitems ){
   if( !( j %in% fixeditemlist) ){
    Rj=Rall[,fromP[j]:toP[j],drop=0]
    Nj=rowSums(Rj)
    # Print(j,Rj,Nj)
    V[,fromP[j]:toP[j]]=Rj/Nj
    if( monotone ){
     # monotone
     if( any( diff(Rj[,2]/Nj) < 0 ) ){
      pij=mbinreg( Rj[,2], Nj, 1:nclass )
      if( any( diff(pij) < 0 ) ){
       Print(j,pij)
       stop()
      }
      V[,fromP[j]:toP[j]]=cbind(1-pij,pij)
      # Print(llll,temp$code)
     }
    }
    # lmlh1=getlmlh( U, group, V, rho, alpha1 )$lmlh
    # Print(llll,j,lmlh,lmlh1,lmlh1-lmlh)
   }
  } # end of j loop

  if( monotone ){
   for( j in 1:nitems ){
    if( any( diff(V[,2*j]) < 0 ) ){
     Print(j, V[,2*j], Rall[,fromP[j]:toP[j],drop=0])
    }
   }
  }
  #  if( any(V<vmin) ) Print(V, fmt="7.3")
  V[V<vmin]=vmin
  V[V>vmax]=vmax

  # Print(V)

  if( debug > 0 ) Print("*top of m-step for class dist")
  # m-step for the class distribution
  locg=which(group==baseform)
  if( estrho == 1 && llll >= 2 ){
   for( g in 1:nG ){
    if( FixedItems != 1 ){
     locg=which(group==g)
     Rg=colSums(H[locg,])+alpha1
     rho[,g]=Rg/sum(Rg)
    }
   }
  } # end of class dist



  if( debug > 0 ) Print("*top of convergence check")
  # convergence:  moved to this place from above 121110
  temp=getlmlh( U, group, V, rho, alpha1 )
  lmlh=temp$lmlh; H=temp$H; R=temp$R; logP=temp$logP; rm(temp)
  lmlhimpr=(lmlh-lmlhp)/abs(lmlhp)
  maxadp=max(abs(paramnump-param[nonfixeditemlist,4:(4+maxnp-1)]), na.rm=1)
  maxadrho=max(abs(rhop-rho))
  if( print >= 2 ) Print( llll, lmlh, lmlhp, lmlhimpr,maxadp,maxadrho
                          , fmt=c("i4",".6") )
  if( lmlhimpr < eps  && maxadp < epsd  && maxadrho < epsd ){
   converged=1
   break
  }

  iter_hist=rbind(iter_hist,c(llll,lmlh,maxadp,maxadrho))


  ###################################################################(2)
  if( SQUAREM  &  llll >= nSQUAREM ){
   # inline SQUAREM
   pp=c( V )
   # Print(locSQEM,locms, nparam1, nparam, param)
   res=iSQUAREM( pp
                 , enforce_constraints=const, badness_of_fit=fit
                 , SQUAREM=SQUAREM, minalpha=minalpha, maxalpha=maxalpha
                 , bof_value=NULL, always=always, reset1=reset1, reset2=reset2
                 , print=0, debug=0 )
   #Print(pp,res$param)
   # reshape
   V=matrix(res$param,nclass)
   colnames(V)=cnv;  rownames(V)=classname
   temp=getlmlh( U, group, V, rho, alpha1 )
   lmlh=temp$lmlh; H=temp$H; R=temp$R; logP=temp$logP; rm(temp)
  }
  ###################################################################(2)


  # next
  lmlhp=lmlh
  rhop=rho
  Vp=V
  paramnump=param[nonfixeditemlist,4:(4+maxnp-1)]


 } # end of llll loop




 iter=llll
 if( print >= 1 && converged )
  cat("\nIteration converged with ", iter, " iterations.  eps= ",eps,"\n")

 if( print > 0 ){
  cat("\nuLRT: parameter estimation of nonparametric IRT models\n")
  cat("\n Input data frame name = ", dfname, "  (compressed=",1-Uin,")\n"
      , sep="")
  cat("\n id varaible =", idvar, ",  group varaible =", groupvar,"\n")
  cat(" # of subjects =", n, ",  # of items =", nitems,"\n")
  cat(" # of groups =", nG, "\n")
  if( nG > 1 ){
   cat("  ");  print(ugroup)
  }
  cat(" percentage of missing data =", pmiss,"\n")
  cat("\n # of classes =", nclass,"\n")
  cat(" estimation of class distribution\n")
  cat("   estrho =", estrho, "\n")
  cat(" prior distribution of rho: Dirichle with ", alpha, "\n")
  cat("   baseform =", baseform, "\n")
  if( FixedItems ){
   cat(" initial item parameter values of the following items are given.\n")
   print(fixeditemlist)
  }
  cat(" total # of parameters to be estimated =", nparam
      ," (",nparam1," + ",nparam2 ,") \n")
  cat(" min and max values of estimated p1 parameter =(", vmin,",", vmax,")\n")
  cat(" max abs value of estimated parameter =", maxabsparam,"\n")
  cat("\n iteration terminated with", iter," iterations\n")
  cat(" max # of iterations in m-steps =", maxiter2,"/",maxiter22
      ," with eps =",eps2, "and nstrict =", nstrict,"\n")
  if( SQUAREM){
   cat("\n SQUAREM =", SQUAREM," will be used after", nSQUAREM
       , "iterations with: \n")
   cat( " minalpha = ", minalpha, ",  maxalpha = ",  maxalpha
        , ",  always = ", always, ",  reset1 = ", reset1
        , ",  reset2 =", reset2, "\n\n")
  }
  cat("\n log marginal likelihood maximized =", lmlh,"\n\n")
  cat(" max # of iterations =", maxiter," with eps =",eps
      , " and epsd =", epsd,"\n")
  cat(" max # of iterations in m-steps =", maxiter2," with eps =",eps2
      , "and nstrict =", nstrict,"\n")
  cat(" print level =", print,"\n")
  cat("\n estimated item parameters\n")
  Print( V, fmt="7.3" )
  cat("\n estimated class distribution \n")
  Print(rho)
  if( print >= 3 ){
   # cat(" estimated # of persons at each of theta points\n")
   # Print(NN)
  }
 }

 if( plot ){
  lmlhf=formatC(lmlh, wid=12, digits=4, format="f")
  title=paste("nonparametric IRT:  # of items =", nitems
              , "   log marginal likelihood =", lmlhf )
  sub=paste("prior Dirichle alpha = ", paste(alpha, collapse=" "), sep="")
  matplot( 1:nclass, V[,toP], type="l"
           , main=title, sub=sub, xlab="class", ylab="prob" )
 }

 if( plot >= 1 ){
  for( g in 1:(nG) ){
   title="latent class distribution"
   if( nG >= 2 ) title=paste(title, " for group ", g,sep="")
   title=paste(title, ":   log marginal likelihood =", lmlhf, sep="")
   sub=paste("prior Dirichle alpha = ", paste(alpha, collapse=" "), sep="")
   barplot(height=rho[,g], names.arg=classname
           , main=title, sub=sub, space=0, las=2)
   }
 }


 return(
  list( V=V, rho=rho, ncat=ncat, type=type, fromP=fromP, toP=toP
        , nitems=nitems, nclass=nclass, n=n
        , converged=converged, lmlh=lmlh,  iter_hist=iter_hist
        , H=H, id=id )
        )


} # end of uLRT












comments(
'

nclass=5

res1 <- uLRT( Uc, maxiter=1000, print=1, nclass=nclass, estrho=1, monotone=0 )

resm1 <- uLRT( Uc, maxiter=1000, print=1, nclass=nclass, estrho=1, monotone=1 )

res1$V[,seq(2,ncol(Uc),2)]
resm1$V[,seq(2,ncol(Uc),2)]


res1 <- uLRT( Uc2, maxiter=1000, SQUAREM=3, print=1, nclass=nclass, estrho=1 )

res2 <- uLRT( Uc2, maxiter=1000, SQUAREM=3, print=1, nclass=nclass, estrho=1
, rho=rep(1/nclass, nclass) )





')

