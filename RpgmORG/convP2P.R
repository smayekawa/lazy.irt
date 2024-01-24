#' Convert Partial Credit Item Parameters in Standard Format (step parameters)
#' to Partial Credit Item Parameters in Nominal Format
#'
#'
#' @param param Item Parameter Data Frame
#' @param print = 1 to print result
#' @param DinP = 1 to include D=1.7 in logistic function
#' @param ctype = 1 to convert to type PN with c-parameters
#'
#'
#' @details
#'     ICRF of category k of item j at theta, namely, \eqn{p_{jk}(\theta)}\cr
#'     is defined as\cr
#' \eqn{p_{jk}(\theta) = Ez_{jk}(\theta)
#' / \sum_{k=0}^{ncat[j]-1} Ez_{jk}(theta)}\cr
#'       where\cr
#' \cr
#'    if type = "P" , \cr
#'     \eqn{Ez_{jk}(\theta) =}
#'    \eqn{exp( 1.7^{DinP} \sum_{h=0}^k a*_j \; (\theta - b*_{jh}) )} \cr
#'                                       \eqn{, k=0,1, ..., ncat[j]-1} \cr
#'    where \eqn{b*_{jk}} is the step parameter with \eqn{b*_{0j}=0}. \cr
#' \cr
#'    if type = "PN" and ctype=0, \cr
#'     \eqn{Ez_{jk}(\theta) = exp( 1.7^{DinP} a_j \; k (\theta - b_{jk}) )} \cr
#'                                       \eqn{, k=0,1, ..., ncat[j]-1} \cr
#'    with \eqn{b_{j0}=0}. \cr
#' \cr
#'    if type = "PN" and ctype=1, \cr
#'     \eqn{Ez_{jk}(\theta) = exp(  a_j k \theta + c_{jk}) )} \cr
#'                                     \eqn{, k=0,1, ..., ncat[j]-1} \cr
#'    with \eqn{c_{j0}=0}. \cr
#' \cr
#'     Note that if \code{DinP=1}, 1.7 will be used. \cr
#' \cr
#'     The step parameter, \eqn{b*_{jk}}, is  the value of theta where
#'      \eqn{P_{jk-1}(\theta)} and \eqn{P_{jk}(\theta)} intersect.\cr
#'
#'
#' @return
#' Partial Credit Item Parameter Data Frame in Nominal Format
#'
#' @examples
#' temp <- irf( paramS1[4,], plot=1 )
#' temp <- irf( convP2PN(paramS1[4,]), plot=1 )
#'
#' temp1 <- convP2PN(paramS1[4,])
#' temp2 <- convPN2P(temp1)
#' Print(paramS1[4,],temp1,temp2)
#'
#' temp1 <- convP2PN(paramS1[4,], ctype=1)
#' temp2 <- convPN2P(temp1, ctype=1)
#' Print(paramS1[4,],temp1,temp2)
#'
#'
#' @export
#'
convP2PN <- function( param, print=0, DinP=1, ctype=0 ){
 # convert the original parametrization ( no 1.7 and step parameters)
 # to x-intercept parameters with 1.7
 # Shin-ichi Mayekawa
 # 20121015
 # DinP option:121109(London)
 # 121113
 # ctype option: 20150614
 # roxygen2: 20231103cot
 #
 #    ICRF of category k of item j at theta, namely, p_{jk}(theta)
 #    is defined as
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0,ncat[j]-1} Ez_{jk}(theta)
 #    where
 #
 #    if type = "PN" and ctype=0,
 #     Ez_{jk}(theta) = exp( 1.7^DinP a_j k (theta - b_{jk}) )
 #                                              , k=0,1, ..., ncat[j]-1
 #    with b_{j0}=0.
 #
 #    if type = "PN" and ctype=1,
 #     Ez_{jk}(theta) = exp(  a_j k theta + c_{jk}) )
 #                                              , k=0,1, ..., ncat[j]-1
 #    with c_{j0}=0.
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

 if( ctype == 0 ){
  b=bs
  b=t( apply( bs,1,cumsum)*(1/(1:ncol(bs))) )
  if( DinP == 0 ) param0[,locparam]=param0[,locparam]/1.7
  param0[,(locparam+1):(locparam+max(ncat)-1)]=b
 }
 else{
  c=-t( apply( bs,1,cumsum) )*param0[,locparam]
  param0[,(locparam+1):(locparam+max(ncat)-1)]=c
  if( DinP == 1 )
   param0[,locparam:(locparam+max(ncat)-1)]=
   param0[,locparam:(locparam+max(ncat)-1)]*1.7
 }
 param0$type="PN"
 param[locP,]=param0

 if( print > 0 ){
  cat( "\n item parameters with D in P =", DinP
       , "  and ctype =", ctype, "\n" )
  Print(param0)
 }

 return( param )

} # end of convP2PN




#' Convert Partial Credit Item Parameters in Nominal Format to
#' Standard Format (step parameters)
#'
#'
#' @param param Item Parameter Data Frame
#' @param print = 1 to print result
#' @param DinP = 1 to include D=1.7 in logistic function
#' @param ctype = 1 to convert to type PN with c-parameters
#'
#'
#' @details
#'     ICRF of category k of item j at theta, namely, \eqn{p_{jk}(\theta)}\cr
#'     is defined as\cr
#' \eqn{p_{jk}(\theta) = Ez_{jk}(\theta)
#' / \sum_{k=0}^{ncat[j]-1} Ez_{jk}(theta)}\cr
#'       where\cr
#' \cr
#'    if type = "P" , \cr
#'     \eqn{Ez_{jk}(\theta) =}
#'    \eqn{exp( 1.7^{DinP} \sum_{h=0}^k a*_j \; (\theta - b*_{jh}) )} \cr
#'                                       \eqn{, k=0,1, ..., ncat[j]-1} \cr
#'    where \eqn{b*_{jk}} is the step parameter with \eqn{b*_{0j}=0}. \cr
#' \cr
#'    if type = "PN" and ctype=0, \cr
#'     \eqn{Ez_{jk}(\theta) = exp( 1.7^{DinP} a_j \; k (\theta - b_{jk}) )} \cr
#'                                       \eqn{, k=0,1, ..., ncat[j]-1} \cr
#'    with \eqn{b_{j0}=0}. \cr
#' \cr
#'    if type = "PN" and ctype=1, \cr
#'     \eqn{Ez_{jk}(\theta) = exp(  a_j k \theta + c_{jk}) )} \cr
#'                                     \eqn{, k=0,1, ..., ncat[j]-1} \cr
#'    with \eqn{c_{j0}=0}. \cr
#' \cr
#'     Note that if \code{DinP=1}, 1.7 will be used. \cr
#' \cr
#'     The step parameter, \eqn{b*_{jk}}, is  the value of theta where
#'      \eqn{P_{jk-1}(\theta)} and \eqn{P_{jk}(\theta)} intersect.\cr
#'
#'
#'
#' @return
#' Partial Credit Item Parameter Data Frame in Standard Format
#'
#' @examples
#' temp1 <- convP2PN(paramS1[4,])
#' temp2 <- convPN2P(temp1)
#' Print(paramS1[4,],temp1,temp2)
#'
#' temp1 <- convP2PN(paramS1[4,], ctype=1)
#' temp2 <- convPN2P(temp1, ctype=1)
#' Print(paramS1[4,],temp1,temp2)
#'
#' @export
#'
convPN2P <- function( param, print=0, DinP=1, ctype=0 ){
 # convert the x-intercept parametrization ( with 1.7 and b-type parameters)
 # to the original step parameters w/o 1.7
 # Shin-ichi Mayekawa
 # 20121015
 # DinP option:121109(London)
 # ctype: 20150614
 # bugfix for ncat=2: 20221205dnc
 # roxygen2: 20231103cot
 #
 #    ICRF of category k of item j at theta, namely, p_{jk}(theta)
 #    is defined as
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0,ncat[j]-1} Ez_{jk}(theta)
 #      where
 #
 #    if type = "PN" and ctype=0,
 #     Ez_{jk}(theta) = exp( 1.7^DinP a_j k (theta - b_{jk}) )
 #                                              , k=0,1, ..., ncat[j]-1
 #    with b_{j0}=0.
 #
 #    if type = "PN" and ctype=1,
 #     Ez_{jk}(theta) = exp(  a_j k theta + c_{jk}) )
 #                                              , k=0,1, ..., ncat[j]-1
 #    with c_{j0}=0.
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

 if( ctype == 1 ){
  # Here, above b is actually the c parameter.
  b=-b/param0[,locparam]*(1/(1:ncol(b)))
  param0[,locparam]=param0[,locparam]/1.7
 }

 bs=b
 for( j in 1:nrow(param0) ){
  if( ncat[j] >= 3 ){
   for( k in 2:(ncat[j]-1) ){
    bs[j,k]=k*b[j,k]-(k-1)*b[j,k-1]
   }
  }
 }

 if( DinP == 0 ) param0[,locparam]=1.7*param0[,locparam]
 param0[,(locparam+1):(locparam+max(ncat)-1)]=bs
 param0$type="P"
 param[locPN,]=param0

 if( print > 0 ){
  cat( "\n item parameters with D in P =", DinP
       , "  and ctype =", ctype, "\n" )
  Print(param0)
 }

 return( param )

} # end of convPN2P





#' Convert Partial Credit Item Parameters in Standard Format  (step parameters)
#'  to Nominal Format
#'
#'
#' @param param Item Parameter Data Frame for Partial Credit Items
#' @param print = 1 to print result
#' @param DinP = 1 to include D=1.7 in logistic function
#'
#' @details
#'     ICRF of category k of item j at theta, namely, \eqn{p_{jk}(\theta)}\cr
#'     is defined as\cr
#' \eqn{p_{jk}(\theta) = Ez_{jk}(\theta)
#' / \sum_{k=0}^{ncat[j]-1} Ez_{jk}(theta)}\cr
#'       where\cr
#' \cr
#'    if type = "P" , \cr
#'     \eqn{Ez_{jk}(\theta) =}
#'    \eqn{exp( 1.7^{DinP} \sum_{h=0}^k a*_j \; (\theta - b*_{jh}) )} \cr
#'                                       \eqn{, k=0,1, ..., ncat[j]-1} \cr
#'    where \eqn{b*_{jk}} is the step parameter with \eqn{b*_{0j}=0}. \cr
#' \cr
#'    if type = "PN" and ctype=0, \cr
#'     \eqn{Ez_{jk}(\theta) = exp( 1.7^{DinP} a_j \; k (\theta - b_{jk}) )} \cr
#'                                       \eqn{, k=0,1, ..., ncat[j]-1} \cr
#'    with \eqn{b_{j0}=0}. \cr
#' \cr
#'    if type = "PN" and ctype=1, \cr
#'     \eqn{Ez_{jk}(\theta) = exp(  a_j k \theta + c_{jk}) )} \cr
#'                                     \eqn{, k=0,1, ..., ncat[j]-1} \cr
#'    with \eqn{c_{j0}=0}. \cr
#' \cr
#'     Note that if \code{DinP=1}, 1.7 will be used. \cr
#' \cr
#'     The step parameter, \eqn{b*_{jk}}, is  the value of theta where
#'      \eqn{P_{jk-1}(\theta)} and \eqn{P_{jk}(\theta)} intersect.\cr
#'
#'
#'
#' @export
#'

convP2N <- function( param, print=0, DinP=1 ){
 # convert the original parametrization ( no 1.7 and step parameters)
 # to a_jk theta + c_jk type nominal response model w/o 1.7
 # where a_jk = a_k * k
 # Shin-ichi Mayekawa
 # 20121122,1202
 # roxygen2: 20231103cot
 #
 #    ICRF of category k of item j at theta, namely, p_{jk}(theta)
 #    is defined as
 #     p_{jk}(theta) = Ez_{jk}(theta) / sum_{k=0}^{ncat[j]-1} Ez_{jk}(theta)
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



