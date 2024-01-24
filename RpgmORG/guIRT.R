#' Item Parameter Estimation of Unidimensional IRT
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
#' @param   ncat     nitems x 1 #  of categories for each item or NULL.
#' @param   type     nitems x 1 vector of item types consisting of:
#'          "Bn" | "G" | "PN" | "P" \cr
#'          If NULL, P will be used.
#' @param itemname nitems x 1 vector of item names or NULL.
#' @param DinP = 0 to exclude 1.7 from logistic function.
#' @param param  nitems x max(ncat)  initial value data frame for
#'  item param\cr
#'         Only the subset of the items can be given.\cr
#'         param$name and the colname(Uc) will be used to match the items.
#' @param msn  nG x 3  initial value matrix of (mean, std, n) for each group
#' @param   fixeditems list of items whose item parameters are to be fixed\cr
#'   to the values given in the param data frame.\cr
#'  Item numbers or item names\cr
#' @param  baseform base form number whose theta distribution is
#' fixed at the values given in msn[baseform,] \cr
#'  If there are fixed items, the values of the item paramters
#'  given in the param data frame must be on the baseform scale.\cr
#' @param   theta    npoints x 1   discrete theta points
#' @param thd NOT used.
#' @param   estmu    = 1  to estimate multigroup theta means
#' @param   estsigma = 1  to estimate mul tigroup theta std
#' @param   npoints # of discrete theta points between (thmin, thmax)
#' @param thmin Minimum value of theta points to be generated.
#' @param thmax Maximum value of theta points to be generated
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
#' Scandinavian Journal of Statistics, Vol. 35: 335-353.
#'
#' @return A list of:
#' param Estimated Item Parameters in a data frame \cr
#' msn Estimated means, standard deviations, and n for each group. \cr
#' theta Discrete theta points used. \cr
#' thd Theta distribution \cr
#' converged = 1 if converged, = 0 otherwise. \cr
#' lmlh Log Marginal LIkelihood maximized \cr
#' iter_hist A list of iteration history \cr
#' H Posterior distribution of theta given U.  \cr
#' NN Estimated # of persons at each of theta points \cr
#' EAP Estimated ability \cr
#' pstdtd Estimated posterior std for theta \cr
#' id Id variable \cr
#'
#' @details
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
#' # Binary Items: simgle group analysis from compressed data: Uc
#' # generate 3PLM data
#' set.seed(1701)
#' Uc <- gendataIRT( 1, paramB2, npoints=500, thdist="rnorm", compress=1 )$U
#' Uc <- as.data.frame(Uc)
#'
#' # 2PLM analysis
#' itemtype <- rep("B2",nrow(paramB2))
#' res1 <- guIRT( Uc, type=itemtype, maxiter=200 )
#'
#' # 3PLM analysis
#' itemtype <- rep("B3",nrow(paramB2))
#' res1 <- guIRT( Uc, type=itemtype, maxiter=200 )
#'
#'
#' # Mixed Type Items: simgle group analysis from compressed data: Uc
#'
#' # generate data
#' set.seed(1701)
#' Uc <- gendataIRT( 1, paramS1, npoints=500, thdist="rnorm"
#'                     , compress=1 )$U
#' Uc <- as.data.frame(Uc)
#' itemtype <- paramS1$type
#' res1 <- guIRT( Uc, type=itemtype, maxiter=20 )
#'
#' # convert compressed data to uncompressed data: U
#' temp=dummy_expand( Uc )
#' U=data.frame( temp$U )
#' ncat=temp$ncat
#' rm(temp)
#' res2 <- guIRT( U=U, ncat=ncat, type=itemtype, maxiter=20 )
#'
#' # should be identical
#' Print(res1$param, res2$param)
#'
#'
#' # fixed parameter values
#' paramF <- paramS1[2:3,]; fixeditems=c("Q2","Q3")
#' res1 <- guIRT( Uc, type=itemtype, maxiter=20, param=paramF
#'          , fixeditems=fixeditems )
#'
#'
#' # multi-group analysis with different theta distributions
#' set.seed(1701)
#' indata1 <- gendataIRT( 1, paramS1, npoints=500, thdist="rnorm"
#'                      , compress=1 )$U
#' indata1 <- as.data.frame(indata1,row.names=NULL)
#' indata1 <- data.frame(group="G1",indata1, stringsAsFactors=0
#' ,row.names=NULL)
#' indata2 <- gendataIRT( 1, paramS1, npoints=500, thdist="rnorm", compress=1
#'                      , thmean=1, thstd=1 )$U
#' indata2 <- as.data.frame(indata2,row.names=NULL)
#' indata2 <- data.frame(group="G2",indata2, stringsAsFactors=0
#' ,row.names=NULL)
#' indata12 <- rbind(indata1,indata2)
#' itemtype <- paramS1$type
#' # This will not converge: increase maxiter.
#' res1 <- guIRT( indata12, groupvar="group", type=itemtype, maxiter=10
#' , baseform=1, estmu=1, estsigma=1, minalpha=-2, SQUAREM=3, plot=1 )
#'
#'
#' # example of U matrix input with unequal # of trials per item.
#' set.seed(1701)
#' npoints=2000
#' param=paramS2
#' nitems=nrow(param)
#' ncat=param$ncat
#' Nmat=matrix(sample(1:9,npoints*nitems,replace=1),npoints,nitems)
#' pmiss=0.2
#' for( j in 1:nitems ){
#'   Nmat[sample(1:npoints,npoints*pmiss),j]=0
#' }
#' temp <- gendataIRT( 1, param, Nmat=Nmat, npoints=npoints
#'                     , thdist="rnorm", compress=0 )
#' ncat=temp$ncat
#' UN <- as.data.frame(temp$U)
#' res3 <- guIRT( U=UN, ncat=ncat, maxiter=20 )
#'
#'
#' @export
#'

guIRT <- function( Uc, U=NULL, groupvar=NULL, idvar=NULL
                  , ncat=NULL, type=NULL, itemname=NULL, DinP=1
                  , param=NULL, msn=NULL, fixeditems=NULL, baseform=1
                  , theta=NULL, thd=NULL, psi=NULL, psd=NULL
                  , thetapsi=NULL, thpsd=NULL
                  , estmu=0, estsigma=0
                  , npointth=21, thmin=-4, thmax=4
                  , npointps=21, psmin=0, psmax=9
                  , maxiter=200, eps=1e-7, epsd=1e-6
                  , maxiter2=20, eps2=1e-4, nstrict=9, maxiter22=5
                  , minp1=1e-1, maxabsparam=20
                  , SQUAREM=3, nSQUAREM=1, minalpha=-999, maxalpha=-1
                  , always=1, reset1=0, reset2=1
                  , print=2, plot=0, smallP=0, debug=0 ){
 # unidimensional IRT parameter estimation
 # Shin-ichi Mayekawa
 # 121024,25
 # multi-group: 121026
 # group and id as a column of Uc: 121028,29
 # baseform: 121031
 # print: 121031,1101
 # output lmlh: 121103
 # DinP: 121110(London)
 # epsd, location of conv. check, plot: 121110(London)
 # initial msn: 121112
 # calculate H only in getlmlh: 121112
 # smallP to ordinal_reg: 121115
 # limits to args of exp or log: 121117
 # initial p3=c=0.3 for type="B3" items: 121118
 # group names etc: 121119
 # nominal response model: 121123,24
 # fixed parameters: 121124,25,26
 # inline SQUAREM: 121126,27
 # epsd for maxadp: 121128
 # positive minalpha: 121129
 # nstrict, nSQUAREM: 121130
 # check if the input dataframe exists: 121215
 # estimate baseform theta dist if fixeditems: 121216
 # avoid 3d array and exclude fixed items from SQUAREM: 121216
 # include ms in SQUAREM: 121217
 # calculation of R in getlmlh: 121217
 # new print: 121230
 # maxabsparam for ordinal_reg: 130204
 # baseform cannot be 0: 20150615
 # check if param0 is valid: 20150615
 # iSQUAREM: 20160101,02,20160104
 # U as arg: 20161103,05
 # minp1: 20161115
 # comments of getlmlh: 20161127
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
 #  estmu    =1  to estimate multigroup theta means
 #  estsigma =1  to estimate multigroup theta std
 #
 #  theta    npointth x 1   discrete theta points
 #  npointth, # of discrete theta points between (thmin, thmax)
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
 #     param   estimated parameter data frame
 #     msn     mu, sigma, n  of theta distriution for each group
 #     (theta, thd) theta distribution
 #     H       n x npoints   posterior of theta given data
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
 #   if estmu == 1.
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


 dgbeta <- function( x, paramab, rangex, paramms=NULL ){
  # density of generalized beta
  # Shin-ichi Mayekawa
  # 20160812cot,13cot
  #
  # Args:
  #   x        range of the score
  #   paramab  (a,b) shape parameter vector
  #   rangex   (minx-margin, maxx+margin)
  #            When x consists of integers, margin may be 0.5 .
  #            Margin can be 0 if x is continuous.
  #   paramms  (mean,std) parameter vector.
  #
  # x is in (minx-margin,maxx+margin) and y is usual beta in (0,1)
  #  x=(maxx-minx+2*margin)*y + minx-margin
  #  y=(x-minx+margin)/(maxx-minx+2*margin)
  #
  minscore=rangex[1]; maxscore=rangex[2]
  if( !is.null(paramms) ){
   paramab=ms2ab( paramms, rangex )
  }
  a=paramab[1]; b=paramab[2]
  pdf=dbeta( (x-minscore)/(maxscore-minscore), a, b )
  pd=pdf / (maxscore-minscore)^(a+b-1)

  # for univarage fitting, this is OK.
  # pdf=(x-minscore)^(a-1)*(maxscore-x)^(b-1)

  # Print(x,(x-minscore)/(maxscore-minscore),rangex)

  return( pdf )
 } # end of dgbeta


 ab2ms <- function( ab, rangex ){
  # get (mean,std) from (a,b) of beta distribution
  # Shin-ichi Mayekawa
  # 20160813cot,18
  #

  minscore=rangex[1]; maxscore=rangex[2]

  if( length(ab) == 2 ){
   ab=matrix(ab,2,1)
  }
  nclass=ncol(ab)
  ms=matrix(0,2,nclass)
  colnames(ms)=paste("c",1:nclass,sep="")
  rownames(ms)=c("mean","std")

  for( k in 1:nclass ){
   a=ab[1,k]; b=ab[2,k]
   meangb=(maxscore-minscore)*a/(a+b) + minscore
   vargb=(maxscore-minscore)^2*(a*b)/(a+b)^2/(a+b+1)
   if( vargb <= 0 ) vargb=0
   stdgb=sqrt(vargb)
   ms[,k]=c( meangb, stdgb )
  }

  # if( nclass == 1 ) ms=c(ms)
  return( ms )

 } # end of ab2ms


 ms2ab <- function( ms, rangex ){
  # get (alpha, beta) from (mean, std) of beta distribution
  # Shin-ichi Mayekawa
  # 20160812cot,13cot
  # (m,s) -> ms: 20160813cot
  # return NULL if quad > 0: 20160817
  #

  minscore=rangex[1]; maxscore=rangex[2]

  if( length(ms) == 2 ){
   ms=matrix(ms,2,1)
  }
  nclass=ncol(ms)
  res=matrix(0,2,nclass)
  colnames(res)=paste("c",1:nclass,sep="")
  rownames(res)=c("alpha","beta")

  for( k in 1:nclass ){
   m=ms[1,k]; v=ms[2,k]^2

   # is valid?
   # quad=(m^2-(maxscore+minscore)*m+maxscore*minscore+v)/(maxscore-minscore)
   # d=(maxscore+minscore)^2-4*(maxscore*minscore+v)
   # Print(k,m,v,quad)
   # if( quad >= 0 ){
   # cat("\nInput parameter (m,v) are invalid.\n")
   # Print( m,v,quad )
   # return(NULL)
   # }
   #Print(k,m,v,d)
   # translate
   quad=(m^2-(maxscore+minscore)*m+maxscore*minscore+v)/(maxscore-minscore)
   alpha=-(m-minscore)*quad/v
   beta=(m-maxscore)*quad/v
   if( alpha <= 0 | beta <= 0 ){
    cat("\nInput parameter (m,v) are invalid.\n")
    Print( m,v,quad,alpha,beta )
    # return(NULL)
   }
   res[,k]=c(alpha,beta)

   # check
   # mm=(maxscore-minscore)*alpha/(alpha+beta)+minscore
   # vv=(maxscore-minscore)^2*(alpha*beta)/(alpha+beta)^2/(alpha+beta+1)
   # Print(m,v,mm,vv )
  }
  # Print(rangex,d,ms,res)
  # if( nclass == 1 ) res=c(res)
  return( res )

 } # end of ms2ab








 getlmlh <- function( U, group, param, theta, thd ){
  # calculation of log marginal likelihood and the posterior expectation
  # Shin-ichi Mayekawa
  # 20121024
  # output H: 121112
  # limits of exp or log: 121117,18,19
  # return R: 121217
  # comments: 20161127
  #
  # Args:
  #   U      n x sum(ncat)  expanded data
  #          U[i,fromP[j]:toP[j]]=0 if Uc[i,j] is missing
  #   group  n x 1
  #   param  parameter data frame: name, type, ncat, p1, p2, ...
  #   theta  npoint x 1   discrete theta points
  #   thd    npoint x nG  prob at theta points for each group (sums to unity)
  #
  # Values:
  #   lmlh   log marginal likelihood
  #   logP   npoints x sum(ncat) log of IRF
  #   H      n x npoints  posterior expectation of theta given obs
  #   R      npoints x sum(ncat) x nG array of post exp of U at theta
  #
  #  Note
  #     min arg to 1/1e-x is 308
  #     min arg to log is 1e-308:
  #     max arg to exp is 709.782:
  #

  nG=length(unique(group))
  logP=log( irf( param, theta, DinP=DinP, print=0 )$ICRF )
  H=matrix(0,nrow(U),length(theta))
  lmlh=0
  R=array(NA,c(ncol(H),ncol(U),nG))
  for( g in 1:nG ){
   locg=which(group==g)
   # eULP=exp(U[locg,,drop=0]%*%t(logP)) * matrix(1,length(locg))%*%t(thd[,g])
   ULP=U[locg,,drop=0]%*%t(logP)
   ULP[ULP > 709]=709
   eULP=exp(ULP)*matrix(1,length(locg))%*%t(thd[,g])
   rseULP=rowSums(eULP)
   rseULP[rseULP < 1e-307]=1e-307
   lmlh=lmlh+sum( log( rseULP ) )
   H[locg,]=eULP/rseULP
   R[,,g]=t(H[locg,])%*%U[locg,,drop=0]
  }
  rm(ULP,eULP,rseULP)

  return( list(lmlh=lmlh, logP=logP, H=H, R=R) )

 } # end of getlmlh



 fit <- function( newp ){
  # returns minus log marginal likelihood etc for iSQUAREM
  # Shin-ichi Mayekawa
  # 20160102,03
  #

  # Everything except for "newp" comes from the environment
  # in which this function is defined. i.e., guIRT.

  # reshape parameters
  paramnum[locSQEM]=newp[1:nparam1]
  param0[nonfixeditemlist,4:(4+maxnp-1)]=paramnum
  if( length(locms) > 0 ){
   msn0[locms]=newp[(nparam1+1):nparam]
   thd0=matrix(0,npointth,nG)
   for( g in 1:nG ){
    temp=exp(-0.5*( (theta-msn[g,1])/msn[g,2] )^2);
    thd0[,g]=temp/sum(temp)
   }
  }
  # evaluate lmlh
  temp=getlmlh( U, group, param0, theta, thd0 )
  names(temp)=c("critval",names(temp)[-1])
  temp[[1]]=-temp[[1]]
  return( temp )
 } # end of fit


 const <- function( newp ){
  # enforces constraints on the parameters for iSQUAREM
  # Shin-ichi Mayekawa
  # 20160102,03

  # Everything except for "newp" comes from the environment
  # in which this function is defined. i.e., guIRT.

  # reshape parameters
  paramnum[locSQEM]=newp[1:nparam1]
  param0[nonfixeditemlist,4:(4+maxnp-1)]=paramnum
  # check if param0 is valid
  loc_a=which( param0$p1 <= minp1 )
  if( length(loc_a) > 0 ) param0[loc_a,]$p1=minp1
  loc_c=which( type == "B3"  &  param0$p3 <= 0 )
  if( length(loc_c) > 0 ) param0[loc_c,]$p3=0
  # reshape
  newp[1:nparam1]=c(as.matrix(param0[nonfixeditemlist,4:(4+maxnp-1)])[locSQEM])
  return( newp )
 } # end of const


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
   cat("\n\nerror1(guIRT) **** data frame ", dfname
       , " does not exist. ****\n\n")
   return(NULL)
  }
  else if( !is.data.frame(Uc) ){
   cat("\n\nerror1(guIRT) **** ", dfname, " is not a data frame. ****\n\n")
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
    cat("\n\n error1(guIRT): group variable does not exist.\n\n\n")
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
    cat("\n\n error1(guIRT): ID variable does not exist.\n\n\n")
    Print(idvar)
    return( NULL )
   }
   #
   id=Uc[,locidvar]
   Uc=Uc[,-locidvar,drop=0]
  }

  # base form
  if( baseform <= 0  |  nG < baseform ){
   cat("\n\n error1(guIRT): baseform must be in [1,",nG, "]\n\n\n")
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
   cat("\n\nerror1(guIRT): *** When U is given, ncat must be given. *** \n")
   return(NULL)
  }
  if( is.null(itemname) ){
   itemname=paste("Q",1:length(ncat),sep="")
  }

  # chech if dataset exists
  dfname=deparse(substitute(U))
  dfname=unlist( strsplit(dfname,"[",fixed=1) )[1]
  if( !exists(dfname) ){
   cat("\n\nerror1(guIRT) **** data frame ", dfname
       , " does not exist. ****\n\n")
   return(NULL)
  }
  else if( !is.data.frame(U) ){
   cat("\n\nerror1(guIRT) **** ", dfname, " is not a data frame. ****\n\n")
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
    cat("\n\n error1(guIRT): group variable does not exist.\n\n\n")
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
    cat("\n\n error1(guIRT): ID variable does not exist.\n\n\n")
    Print(idvar)
    return( NULL )
   }
   #
   id=U[,locidvar]
   U=U[,-locidvar,drop=0]
  }

  # base form
  if( baseform <= 0  |  nG < baseform ){
   cat("\n\n error1(guIRT): baseform must be in [1,",nG, "]\n\n\n")
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
 ####################### Hire, U is ready. ##############################
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
   cat("\n\n error1(guIRT): Item type does not match the # of categories.\n\n")
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
  Print(  nitems, nG, npointth, npointps, "/", type, ncat, np )

 # initial parameter data frame
 paramORG=param
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
 nparam1=sum(np)-sum(np[fixeditemlist])
 nparam2=0
 if( estmu == 1 ) nparam2=nparam2+nG-1
 if( estsigma == 1 ) nparam2=nparam2+nG-1
 if( FixedItems ){
  if( estmu == 1 ) nparam2=nparam2+1  # for the baseform
  if( estsigma == 1 ) nparam2=nparam2+1  # for the baseform
 }
 nparam=nparam1+nparam2


 # generate theta
 if( is.null(theta) ){
  theta=matrix(seq(thmin,thmax,length.out=npointth),,1)
 }
 else{
  npointth=length(as.vector(theta))
 }
 # theta distribution
 thd=matrix(0,npointth,nG)
 rownames(thd)=format(theta,2)
 if( is.null(msn) ){
  # make all N(0,1)
  msn=matrix(c(0,1,1),nG,3,byrow=1)
 }
 for( g in 1:nG ){
  temp=exp(-0.5*( (theta-msn[g,1])/msn[g,2] )^2);
  thd[,g]=temp/sum(temp)
 }
 msn[,3]=nobs
 colnames(msn)=c("mu","sigma","n")
 rownames(msn)=paste("group",1:nG,sep="")

 # generate psi
 if( is.null(psi) ){
  psi=matrix(seq(psmin,psmax,length.out=npointps),,1)
 }
 else{
  npointps=length(as.vector(psi))
 }
 # psi distribution
 psd=matrix(0,npointps,nG)
 rownames(psd)=format(psi,2)
 if( is.null(msnps) ){
  # make all N(0,1)
  msnps=matrix(c(1,8,1),nG,3,byrow=1)
 }
 rangeps=c(0,5)
 for( g in 1:nG ){
  temp=dgbeta( psi, paramab=msnps[g,], rangex=rangeps )
  psd[,g]=temp/sum(temp)
 }
 msn[,3]=nobs
 colnames(msnps)=c("a","b","n")
 rownames(msnps)=paste("group",1:nG,sep="")






 # what to estimate
 estms=matrix(0,nG,2)
 if( (estmu == 1  ||  estsigma == 1) )
  for( g in 1:nG )
   if( FixedItems == 1  ||  g != baseform )
    estms[g,]=c(estmu,estsigma)
 locms=which(estms == 1)


 if( print > 0 ){
  cat("\n\n\nuIRT: parameter estimation of unidimensional IRT models\n")
  cat("\n Input data frame name = ", dfname, "  (compressed=",1-Uin,")\n"
      , sep="")
  cat("\n id varaible =", idvar, ",  group varaible =", groupvar,"\n")
  cat(" # of subjects =", n, ",  # of items =", nitems,"\n")
  cat(" # of groups =", nG, "\n")
  if( nG > 1 ){
   cat("  ");  print(ugroup)
  }
  cat(" percentage of missing data =", pmiss,"\n")
  cat(" # of theta points =", npointth, " in [", thmin, ",", thmax,"]","\n")
  cat(" estimation of theta distribution\n")
  cat("   estmu =", estmu, ",  estsigma =", estsigma, "\n")
  cat("   baseform =", baseform, "\n")
  if( FixedItems ){
   cat(" The item paramters of the following items are fixed. \n")
   print(fixeditemlist)
  }
  cat(" total # of parameters to be estimated =", nparam
      ," (",nparam1," + ",nparam2 ,") \n")
  cat(" minimum value of estimated p1 parameter =", minp1,"\n")
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
  cat(" D in P =", DinP,"\n")
  if( length(initemlist) > 0 ){
   cat(" initial item parameter values of the following items are given.\n")
   print(initemlist)
  }
  if( print >= 2 ){
   cat(" inital parameter value \n")
   Print(param)
   Print(msn)
   if( print >= 3 ){
    Print(thd,digits=1)
   }
  }
 }


 #initial
 temp=getlmlh( U, group, param, theta, thd )
 lmlhp=temp$lmlh
 logP=temp$logP
 H=temp$H
 R=temp$R
 rm(temp)
 paramnum=as.matrix(param[nonfixeditemlist,4:(4+maxnp-1)])
 paramnump=paramnum
 NN=matrix(0,npoints+1,nG); rownames(NN)=c(format(theta,2),"total")


 if( print > 0 ) Print("initial", lmlhp)
 converged=0
 msnp=msn
 alpha=-1


 ###################################################################(1)
 # store param
 # initialize inline SQUAREM
 rmv(iSQUAREM)
 pp=c( as.matrix(param[nonfixeditemlist,4:(4+maxnp-1)])[locSQEM], msn[locms] )
 iSQUAREM=generate_iSQUAREM( pp )
 param0=param
 msn0=msn
 thd0=thd
 ###################################################################(1)


 iter_hist=matrix(c(0,lmlhp,0,0,-1),1,5)
 colnames(iter_hist)=c("llll","lmlh","maxadp","maxadms","alpha")



 # main iteration
 for( llll in 1:maxiter ){


  # e-step
  # expected value of U at theta point
  # Rall is npoints x sum(ncat)
  Rall=apply(R,c(1,2),sum)

  if( llll <= nstrict ) maxiter21=maxiter2
  else maxiter21=maxiter22


  if( debug > 0 ) Print("*top of m-step for items")
  # m-step for item parameters
  for( j in 1:nitems ){
   if( !( j %in% fixeditemlist) ){
    Rj=Rall[,fromP[j]:toP[j]]
    paramj=param[j,]
    eps21=eps2
    printx=print-10
    temp=ordinal_reg( Rj, theta, type=type[j], param=paramj, DinP=DinP
                      , maxiter=maxiter21, eps=eps21
                      , maxabsparam=maxabsparam, minp1=minp1
                      , print=printx, smallP=smallP )
    param[j,1:ncol(temp$param)]=temp$param
   }
  } # end of j loop


  if( debug > 0 ) Print("*top of m-step for theta dist")
  # m-step for the theta distribution
  locg=which(group==baseform)
  # NN[1:npoints,baseform]=rowSums( t(H[locg,,drop=0])%*%U[locg,,drop=0] )
  NN[1:npoints,baseform]=rowSums( R[,,baseform] )
  if( (estmu == 1  ||  estsigma == 1) && llll >= 2 ){
   for( g in 1:nG ){
    if( FixedItems == 1 ||  g != baseform ){
     locg=which(group==g)
     Rg=R[,,g]
     Nj=rowSums(Rg)
     NN[1:npoints,g]=Nj
     temp=smn( theta, Nj, maxiter=maxiter2, eps=eps2, print=print-10
               , estmu=estmu, estsigma=estsigma
               , mu=msn[g,1], sigma=msn[g,2] )
     thd[,g]=temp$P
     msn[g,1]=temp$mu
     msn[g,2]=temp$sigma
    }
   }
   if( nG > 1) NN[npoints+1,]=colSums(NN[1:npoints,])
  } # end of theta dist

  if( debug > 0 ) Print("*top of convergence check")
  # convergence:  moved to this place from above 121110
  temp=getlmlh( U, group, param, theta, thd )
  lmlh=temp$lmlh; H=temp$H; R=temp$R; logP=temp$logP; rm(temp)
  lmlhimpr=(lmlh-lmlhp)/abs(lmlhp)
  maxadp=max(abs(paramnump-param[nonfixeditemlist,4:(4+maxnp-1)]), na.rm=1)
  maxadms=max(abs(msnp-msn))
  if( print >= 2 ) Print( llll, lmlh, lmlhp, lmlhimpr,maxadp,maxadms
                          , fmt=c("i4",".6") )
  if( lmlhimpr < eps  && maxadp < epsd  && maxadms < epsd ){
   converged=1
   break
  }

  iter_hist=rbind(iter_hist,c(llll,lmlh,maxadp,maxadms,alpha))


  ###################################################################(2)
  if( SQUAREM  &  llll >= nSQUAREM ){
   # inline SQUAREM
   pp=c(
    as.matrix(param[nonfixeditemlist,4:(4+maxnp-1)])[locSQEM], msn[locms] )
   # Print(locSQEM,locms, nparam1, nparam, param)
   res=iSQUAREM( pp
                 , enforce_constraints=const, badness_of_fit=fit
                 , SQUAREM=SQUAREM, minalpha=minalpha, maxalpha=maxalpha
                 , bof_value=NULL, always=always, reset1=reset1, reset2=reset2
                 , print=print, debug=0 )
   #Print(pp,res$param)
   # reshape
   paramnum[locSQEM]=res$param[1:nparam1]
   param[nonfixeditemlist,4:(4+maxnp-1)]=paramnum
   if( length(locms) > 0 ){
    msn[locms]=res$param[(nparam1+1):nparam]
    thd=matrix(0,npointth,nG)
    for( g in 1:nG ){
     temp=exp(-0.5*( (theta-msn[g,1])/msn[g,2] )^2);
     thd[,g]=temp/sum(temp)
    }
   }
   temp=getlmlh( U, group, param, theta, thd )
   lmlh=temp$lmlh; H=temp$H; R=temp$R; logP=temp$logP; rm(temp)
  }
  ###################################################################(2)


  # next
  lmlhp=lmlh
  paramnump=param[nonfixeditemlist,4:(4+maxnp-1)]
  msnp=msn



 } # end of llll loop




 iter=llll
 if( print >= 1 && converged )
  cat("\nIteration converged with ", iter, " iterations.  eps= ",eps,"\n")

 if( print > 0 ){
  cat("\nuIRT: parameter estimation of unidimensional IRT models\n")
  cat("\n Input data frame name = ", dfname, "  (compressed=",1-Uin,")\n"
      , sep="")
  cat("\n id varaible =", idvar, ",  group varaible =", groupvar,"\n")
  cat(" # of subjects =", n, ",  # of items =", nitems,"\n")
  cat(" # of groups =", nG, "\n")
  if( nG > 1 ){
   cat("  ");  print(ugroup)
  }
  cat(" percentage of missing data =", pmiss,"\n")
  cat("\n # of theta points =", npointth, " in [", thmin, ",", thmax,"]","\n")
  cat(" estimation of theta distribution\n")
  cat("   estmu =", estmu, ",  estsigma =", estsigma, "\n")
  cat("   baseform =", baseform, "\n")
  if( FixedItems ){
   cat(" initial item parameter values of the following items are given.\n")
   print(fixeditemlist)
  }
  cat(" total # of parameters to be estimated =", nparam
      ," (",nparam1," + ",nparam2 ,") \n")
  cat(" minimum value of estimated p1 parameter =", minp1,"\n")
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
  cat("\n estimated item parameters (DinP=", DinP,") \n")
  Print( param, fmt="7.3" )
  cat("\n estimated mean and std of theta distribution \n")
  Print(msn)
  if( print >= 3 ){
   cat(" estimated # of persons at each of theta points\n")
   Print(NN)
  }
 }

 if( plot >= 1  &  nG > 1 ){
  maxthd=max(thd)
  title="theta distribution for each group"
  sub=NULL
  for( g in 1:(nG) ){
   plot(theta,thd[,g]
        , xlim=c(min(theta),max(theta)), ylim=c(0,maxthd), type="l", ylab="")
   par(new=1)
  }
  plot(theta,thd[,nG], main=title, sub=sub
       , xlim=c(min(theta),max(theta)), ylim=c(0,maxthd), type="l", ylab="prb")
  par(new=0)
 }

 EAP=H%*%theta
 poststd=sqrt(H%*%theta^2-EAP^2)
 return(
  list( param=param, msn=msn, theta=theta, thd=thd
        , converged=converged, lmlh=lmlh,  iter_hist=iter_hist
        , H=H, NN=NN, EAP=EAP, pstdtd=poststd, id=id )
 )


} # end of guIRT












comments(
'


# 2PLM Items: simgle group analysis from compressed data: Uc
set.seed(1701)
Uc <- gendataIRT( 1, paramB1, npoints=500, thdist="rnorm", compress=1 )$U
Uc <- as.data.frame(Uc)
res1 <- guIRT( Uc, maxiter=200 )


# Binary Items: simgle group analysis from compressed data: Uc
# generate 3PLM data
set.seed(1701)
Uc <- gendataIRT( 1, paramB2, npoints=500, thdist="rnorm", compress=1 )$U
Uc <- as.data.frame(Uc)

# 2PLM analysis
itemtype <- rep("B2",nrow(paramB2))
res1 <- guIRT( Uc, param=paramC, type=itemtype, maxiter=200 )

# 3PLM analysis
itemtype <- rep("B3",nrow(paramB2))
res1 <- guIRT( Uc, type=itemtype, maxiter=200 )












# simgle group analysis
set.seed(1701)
indata <- gendataIRT( 1, paramS1, npoints=500, thdist="rnorm", compress=1 )$U
indata <- as.data.frame(indata)
itemtype <- paramS1$type
res1 <- guIRT( indata, type=itemtype, maxiter=200, SQUAREM=3, always=1
              , reset1=0, reset2=1 )



# simgle group analysis
set.seed(1701)
indata <- gendataIRT( 1, paramS1, npoints=500, thdist="rnorm", compress=1 )$U
indata <- as.data.frame(indata)
itemtype <- paramS1$type
res1 <- guIRT( indata, type=itemtype, maxiter=200 )


# fixed parameter values
paramF <- paramS1[2:3,]; fixeditems=c("Q2","Q3")
res1 <- guIRT( indata, type=itemtype, maxiter=20, param=paramF
         , fixeditems=fixeditems )




# multi group analysis
set.seed(1701)
indata1 <- gendataIRT( 1, paramS1, npoints=500, thdist="rnorm", compress=1 )$U
indata1 <- as.data.frame(indata1,row.names=NULL)
indata1 <- data.frame(group="G1",indata1, stringsAsFactors=0,row.names=NULL)
indata2 <- gendataIRT( 1, paramS1, npoints=500, thdist="rnorm", compress=1
                     , thmean=1, thstd=1 )$U
indata2 <- as.data.frame(indata2,row.names=NULL)
indata2 <- data.frame(group="G2",indata2, stringsAsFactors=0,row.names=NULL)
indata12 <- rbind(indata1,indata2)
itemtype <- paramS1$type
res1 <- guIRT( indata12, "group", type=itemtype, maxiter=200, baseform=1
            , estmu=1, estsigma=1, minalpha=-2, SQUAREM=3, plot=1 )


')

