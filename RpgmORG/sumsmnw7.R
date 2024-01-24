#' Distribution of the Weighted Sum of Several Independent
#' Scored Multinomial Distributions
#'
#'
#' @param P matrix of probabilities   (max # of categories x # of r.v.)
#' @param V matrix of domain values   (max # of categories x # of r.v.) or NULL
#' @param w vector of weights         (1 x # of r.v.)
#' @param ncat max # of categories       (1 x # of r.v)    or NULL
#' @param compress = 1 to remove the zero probability categories
#' @param print = 1 to print result
#' @param plot = 1 to plot result
#' @param debug = 1 to print intemediate result
#'
#' @details
#'  ( V[,i], P[,i], w[i] ), i=1,2,...,n is the set of\cr
#'  ( domain or category weight ,probability, and weight ) for the i-th r.v.\cr
#' \cr
#' Non integer V and w will be first converted to integer by linear
#' transformation and converted back at the very end. \cr
#' \cr
#'  ncat[i] = max # of categories for the i-th r.v and\cr
#'  P[(ncat[i]+1):nrow(P),i] == NA\cr
#'\cr
#'   This program calculats the distrobution of\cr
#'     \eqn{X = sum_{i=1}^n w[i] X_i} \cr
#'    where\cr
#'     \eqn{X_i} is distributed as Scored Multinomial with (V[,i],P[,i])
#'  ,i=1,2,...,n\cr
#'     \eqn{V[1,i] <=  X_i  <= V[ncat[i],i]  or  0  <=  X_i  <=  ncat[i]}\cr
#'   That is, this program calculates the probability that\cr
#'     \eqn{Pr( X = a )} ,\cr
#'    where\cr
#'     \eqn{sum_{i=1}^n V[1,i]*w[i] <= a <=  sum_{i=1}^n V[ncat[i],i]*w[i]}\cr
#'
#' @return A matrix of (score, prob)
#'
#'
#' @examples
#' # category x variable matrix of probability:  colSums(P)=c(1,1,1...)
#' P <- matrix(c(1,2,3,2,   1,2,1,0,  1,2,0,0), 4,3)
#' P <- t(t(P)/colSums(P))
#' ncat <- c(4,3,2)
#' # category x variable matrix of category weight
#' V <- NULL
#' # variable weight
#' w <- c(.5,1,1)
#' res <- sumsmnw( P, V, w, compress=0, print=1, plot=1, ncat=ncat )
#'
#' # category x variable matrix of category weight
#' V <- matrix(c(0,1,2,3,   1,2,3,0,  1,1,0,0), 4,3)
#' # variable weight
#' w <- c(.5,1,1)
#' res <- sumsmnw( P, V, w, compress=0, print=1, plot=1, ncat=ncat )
#'
#' @references Mayekawa, S., & Arai, S. (2008).
#' Distribution of the Sum of Scored Multinomial Random Variables
#' and Its Application to the Item Response Theory.
#' In K. Shigemasu, A. Okada, T.Imaizumi, & T. Hoshino (Eds.)
#' New Trends in Psychometrics. Tokyo: University Academic Press.
#'
#' @export
#'

sumsmnw <- function( P, V=NULL, w=rep(1,ncol(P)), ncat=NULL, compress=0
                   , print=0, plot=0, debug=0 ){
 # sum of independent scored multinomial random variables
 # Shin-ichi Mayekawa
 # iml version: 000823 -- 080106
 # 120204,05
 # use V%*%diag(w) : 120205
 # bugfix: 120206
 # plot: 120207
 # print ncat: 120210
 # comments added: 130804
 # bug fix: 20150305
 # ncat: 20150305,07
 # roxgen: 20171205
 #
 #
 # Args:
 #
 #  P     matrix of probabilities   (max # of categories x # of r.v.)
 #  V     matrix of domain values   (max # of categories x # of r.v.) or NULL
 #  w     vector of weights         (1 x # of r.v.)
 #  ncat  max # of categories       (1 x # of r.v)    or NULL
 #
 #  ( V[,i], P[,i], w[i] ), i=1,2,...,n is the set of
 #  ( domain ,probability, and weight ) for the i-th r.v.
 #
 #  ncat[i] = max # of categories for the i-th r.v and
 #  P[(ncat[i]+1):nrow(P),i] == NA
 #
 #   This program calculats the distrobution of
 #     X = sum_{i=1}^n w[i] X_i
 #    where
 #     X_i is distributed as Scored Multinomial with (V[,i],P[,i]),i=1,2,...,n
 #     V[1,i] <=  X_i  <= V[ncat[i],i]   or  0  <=  X_i  <=  ncat[i]
 #   That is, this program calculates the probability that
 #     Pr( X = a ) ,
 #    where
 #     sum_{i=1}^n V[1,i]*w[i] <= a <=  sum_{i=1}^n V[ncat[i],i]*w[i]
 #
 #

 # const
 n=ncol(P)   # # of random variables to be added
 maxncat=nrow(P)
 if( is.null( ncat ) ) ncat=colSums(!is.na(P)) # # of categories for each r.v.
 else
  for( i in 1:n )
   if( ncat[i] < maxncat ) P[(ncat[i]+1):maxncat,i]=NA

 # normalize P so that colsums=1
 locna=which( is.na(P) )
 if( length(locna) > 0 ){
  P[locna]=0
 }
 P=t(t(P)/colSums(P))
 if( length(locna) > 0 ){
  P[locna]=NA
 }
 P0=P
 if( debug ) Print(ncat,P)

 # weight for each r.v.
 if( is.matrix(w) ) w=as.vector(w)

 # category score: v-score: 0,1,...,ncat
 if( is.null(V) ) V=matrix(1:nrow(P)-1,nrow(P),ncol(P))
 V0=V

 if( debug ) Print(P,V,n,ncat)

 # combine the cells with duplicate v-scores
 # There ara nuch room for improvement.
 P1=matrix(NA,nrow(P),ncol(P))
 V1=matrix(0,nrow(P),ncol(P)) # cannot have NA in V
 chksum=numeric(n)
 for( k in 1:n ){
  # sort (p[,k],v[,k]) according to v[,k] and find tie info
  pk=P[,k][1:ncat[k]]; vk=V[,k][1:ncat[k]]
  od=order(vk); pk=pk[od]; vk=vk[od]
  tie=as.matrix(table(vk)); nc=nrow(tie)
  to=cumsum(tie); from=to-tie+1
  if( debug ) Print(k,pk,vk,tie,from,to,nc)
  pk1=numeric(nc)
  for( kk in 1:nc )
   pk1[kk]=sum(pk[from[kk]:to[kk]])
  P1[1:nc,k]=pk1; V1[1:nc,k]=unique(vk)
  chksum[k]=sum(pk1)
 }
 ncat1=colSums(!is.na(P1))  # new # of categories for each r.v.
 if( debug ) Print(ncat1,P1,V1,chksum)

 P=P1; V=V1;

 # make V nonnegative integer by subtracting min and mulplying dd
 # maxv=apply(V,2,max,na.rm=1)
 # minv=apply(V,2,min,na.rm=1)
 minv=V[1,]  # min is the first value, not 0 nor NA
 maxv=V[matrix(cbind(ncat1,1:n),,2)]  # max is the last NA value

 # combine V and w into V2
 V2=V-matrix(minv,nrow(V),n,byrow=1)
 V2=V2%*%diag(w)   # if V has NA, this fails.
 if( debug ) Print(minv,maxv,V,V2)

 found=0; dd=0;  maxdec=0;
 for( d in c(1,2,3,4,5,10,15,20,25,40,50) ){
  #  if( all(d*V2-floor(d*V2) == 0) ){
  if( isTRUE(all.equal(d*V2,floor(d*V2))) ){
   V2=d*V2; found=1; dd=d
   break
  }
 }
 if( !found ){
  maxdec=-1
  for( k in 1:n ){
   if( debug ) Print(k,decp(V[,k]))
   maxdec=max(maxdec,decp(V1[,k])[,3])
  }
  if( any( maxdec > 2 ) ){
   cat("\n\nerror:(sumsmnw) category values too small.  V=", V,"\n")
   return( NULL )
  }
  dd=10^maxdec
  V2=dd*V2
 }

 V2[is.na(P)]=NA
 minv2=V2[1,]
 maxv2=V2[matrix(cbind(ncat1,1:n),,2)]
 ncat2=maxv2+1
 if( debug ) Print(ncat1,V2,maxv2)

 # Now, minv2[k]=0 <= V2[,k] <= ncat1[k]*w[k]=maxv2[k]
 # and there are ncat2[k] distinct values of V2[,k].

 # score of the sum
 score=( 0:sum(maxv2) )/dd + sum(minv*w)

 # expand P1 according to V2
 P2=matrix(0,max(maxv2)+1,n)
 rownames(P2)=0:max(maxv2)
 if( debug ) Print(P1,V2,ncat1)
 for( k in 1:n ){
  P2[V2[1:ncat1[k],k]+1,k]=P1[1:ncat1[k],k]
 }

 if( debug ){
  Print(n,ncat,ncat1, maxdec,dd,minv,maxv,minv2,maxv2,w,ncat2)
  Print(P0,P1,P2, "/", V0,V1,V2)
 }


 # main body
 p3=P2[1:ncat2[1],1]
 for( k in 2:n ){
  if( debug ) Print(k,p3,P2[1:ncat2[k],k])
  p3=sumsmnw12( p3,P2[1:ncat2[k],k],1 )[,2] # always use unweighted sum.
 }
 chksum=sum(p3)

 if( compress == 1 ){
  loc=which(p3 > 0)
  p3=p3[loc]; score=score[loc]
 }

 if( print > 0 ){
  offset=minv*w
  cat("\n\nSum of Scored Multinomial Distributions")
  cat("\n # of Scored Multinomials to be added =", n, "\n")
  cat("\nInput category weight and probability matrices and variable weight\n")
  Print(ncat, V0,P0,w)
  if( !isTRUE(all.equal( P0[!is.na(P0)],P1[!is.na(P0)] )) ){
   cat("\nAfter collecting duplicate v-values\n")
   Print(V1,P1,w)
  }
  cat("\nV <- V%*%diag(w) and After expantion of P. \n")
  Print(V2,P2,"/",offset,dd,maxdec,compress)
  cat("\nThe result\n")
  Print(score, p3, chksum, fmt="6.2, 5.3")
 }

 if( plot ){
  barplot(p3, names.arg=score, xlab="score", ylab="probability"
          , main="Plot of the Score Distribution")
 }

 res=cbind(score,p3)
 colnames(res)=c("score","p")
 rownames(res)=NULL
 return( res )


} # end of sumsmnw



#' Distribution of the Weighted Sum of Two Independent
#' Scored Multinomial Distributions
#'
#'
#' @param p1 Probability of the fist random variable
#' @param p2 Probability of the second random variable
#' @param w2 Integer weight to the second random variable
#' @param eps small number
#' @param print = 1 to print result
#'
#' @details
#'    This program calculats the distrobution of\cr
#'      X12 = X1 + w2 X2\cr
#'     where\cr
#'      Xi  is distributed as Scored Multinomial with  pi,  i=1,2.\cr
#'    That is, this program calculates the probability that\cr
#'      p12[a] = Pr( X12 = a ) , a=0,1,..., n12=(n1-1)+w2*(n2-1)\cr
#' \cr
#'  Note\cr
#'    length(pi) = ni = #'  of categories of the i-th r.v.\cr
#'               = mi + 1   where mi is the max score of the i-th r.v.\cr
#' \cr
#'
#' @return A matrix of (score, prob)
#'
#'
#' @examples
#' sumsmnw12( 1:3, 1:2, 2, print=1 )
#'
#' @references Mayekawa, S., & Arai, S. (2008).
#' Distribution of the Sum of Scored Multinomial Random Variables
#' and Its Application to the Item Response Theory.
#' In K. Shigemasu, A. Okada, T.Imaizumi, & T. Hoshino (Eds.)
#' New Trends in Psychometrics. Tokyo: University Academic Press.
#'
#' @export
#'

sumsmnw12 <- function( p1, p2, w2=1, eps=1e-8, print=0 ){
 # sum of two independent scored multinomial random variables
 # Shin-ichi Mayekawa
 # iml version: 000823 -- 080106
 # 120204
 # Z matrix: 120207
 #
 # Args:
 #   p1 and p2   scored multinomial probabilities
 #    the domains are 0:n1=length(p1) and 0:n2=length(p2), resp
 #
 #   w2          integer weight for the 2nd r.v.
 #
 #   This program calculats the distrobution of
 #     X12 = X1 + w2 X2
 #    where
 #     Xi  is distributed as Scored Multinomial with  pi,  i=1,2.
 #   That is, this program calculates the probability that
 #     p12[a] = Pr( X12 = a ) , a=0,1,..., n12=(n1-1)+w2*(n2-1)
 #
 # Note
 #   length(pi) = ni = # of categories of the i-th r.v.
 #              = mi + 1   where mi is the max score of the i-th r.v.
 #

 # const
 if( is.matrix(p1) ) p1=as.vector(p1)
 if( is.matrix(p2) ) p2=as.vector(p2)
 n1=length(p1); n2=length(p2)
 n12=(n1-1)+w2*(n2-1)+1
 score=1:n12-1

 p1=p1/sum(p1)
 p2=p2/sum(p2)

 if( w2 != floor(w2) ){
  cat("\n\nerror:(sumsmnw12) w2 must be an integer.  w2=", w2,"\n")
  return( NULL )
 }

 # main body
 Z=matrix(0,n12,n2)
 rownames(Z)=score
 for( j in 1:n2 ){
  Z[,j]=c( rep(0,w2*(j-1)), p2[j]*p1, rep(0,w2*(n2-j)) )
 }
 p12=rowSums(Z)
 chksum=sum(p12)

 if( abs(chksum-1) > eps ){
  cat("\n\nwarning(sumsmnw12): Big Trouble!! \n")
  cat(" Sum of probabilities calculated is not equa to 1.", chksum, "\n\n")
 }
 res=cbind(score,p12)
 colnames(res)=c("score","p")
 rownames(res)=rep("",n12)

 if( print > 0 ){
  cat("\n\nSum of Two Scored Multinomial Random Variables\n")
  Print(n1,p1,n2,p2,w2, digits=3)
  if( print >= 2 ){
    cat("\n Z matrix whose row sum is the desired probability.\n")
    print(Z,digits=3)
   }
  Print(score,p12,chksum)
  }

 return( res )

} # end of sumsmnw12




# 20180703
#

min=0; max=5
x=min:max
mean=(max+min)/2
std=(max-min)/2/2
p=graded_prob( x, mean, std, method=1 )

P=cbind(x1=p, x2=p)
V=cbind(x1=x, x2=x)
w=c(1,1)
# ncat=c(max,max)+1

res=sumsmnw( P, V, w, compress=0, print=1, plot=1 )
ms=mands( res[,1], res[,2])
Print(ms)



















comments(
'
         V=matrix(c(
         1, 1,
         1, 0,
         1, 1
         ),3,byrow=1)
         P=matrix(c(
         0.353, 0.042,
         0.248, 0.506,
         0.399, 0.452
         ),3,byrow=1)
         w=c(3,2)

         res=sumsmnw( P, V, w=w, compress=0, print=2, debug=1 )
         Print(res, digits=2)
')




comments(
'
         p1=c(1,2,3,4,NA,NA); p1=p1/sum(p1[!is.na(p1)])
         v1=c(0,1,2,3,NA,NA)
         p2=c(2,3,4,3,2,1); p2=p2/sum(p2[!is.na(p2)])
         v2=c(0,1,2,3,4,5)
         p3=c(2,3,NA,NA,NA,NA); p3=p3/sum(p3[!is.na(p3)])
         v3=c(0,1,2,3,4,5)-2

         P=cbind(p1,p2,p3)
         V=matrix(1:nrow(P)-1,nrow(P),ncol(P))
         V[is.na(P)]=NA


         #V[1,1]=1
         #V[1,2]=-1
         #V[3,1]=-2

         #V=V/2
         #V[,3]=V[,3]/3

         nc=6; nrv=10
         P=matrix(1:nc,nc,nrv)
         V=P-1
         P=P%*%diag(1/colSums(P))

         Print( P,V )


         res=sumsmnw( P, V, compress=1, print=1, plot=1, debug=0 )
         Print(res, digits=2)

 ')


comments(
'
# category x variable matrix of probability:  colSums(P)=c(1,1,1...)
P=matrix(c(1,2,3,2,1,   1,2,1,0,0,  1,2,0,0,0), 5,3)
P=t(t(P)/colSums(P))
# category x variable matrix of category weight
V=matrix(c(0,1,2,3,4,   0,1,2,0,0,  0,1,0,0,0), 5,3)
# variable weight
w=c(.5,1,1)

res=sumsmnw( P, V, w, compress=1, print=1, plot=1, debug=0 )
Print(res)

# negative category x variable matrix of category weight
V=matrix(c(0,1,2,3,4,   0,1,2,0,0,  0,1,0,0,0), 5,3)
V[,1]=V[,1]-1
# decimal variable wieght
w=c(1,0.5,1)
res=sumsmnw( P, V, w, compress=1, print=1, plot=1, debug=0 )
Print(res)
')


comments(
 '
# category x variable matrix of probability:  colSums(P)=c(1,1,1...)
P=matrix(c(1,2,3,2,   1,2,1,0,  1,2,0,0), 4,3)
P=t(t(P)/colSums(P))
ncat=c(4,3,2)
# category x variable matrix of category weight
V=NULL
# variable weight
w=c(.5,1,1)
res=sumsmnw( P, V, w, compress=0, print=1, plot=1, ncat=ncat )
Print(res)

# category x variable matrix of category weight
V=matrix(c(0,1,2,3,   1,2,3,0,  1,1,0,0), 4,3)
# variable weight
w=c(.5,1,1)
res=sumsmnw( P, V, w, compress=0, print=1, plot=1, ncat=ncat )
Print(res)
')

