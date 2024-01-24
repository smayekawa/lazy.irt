
get_monotone <- function( r, n, A, bb=NULL, print=0 ){
 # find monotone pi which maximize the binominal likelihood
 # Shin-ichi Mayekawa
 # 20170621
 #
 # pi = logistic( A %*% expb )
 # where expb=c( b[1],exp(b[-1]) )
 #

 get_ll <- function( bb ){
  # binominal log lihelihood for nlm
  bb1=c(bb[1],exp(bb[-1]))
  pi=logistic1( A%*%bb1 )
  ll=sum(r*log(pi) + (n-r)*log(1-pi))
  return( -ll )
 } # end of get_ll

 # # of elements
 np=length(c(r))

 # unrestricted mle of pi
 pi0=r/n
 ll0=sum(r*log(pi0) + (n-r)*log(1-pi0))

 if( print ){
  pi00=isoreg(pi0)$yf
  ll00=sum(r*log(pi00) + (n-r)*log(1-pi00))
 }

 # initial value of bb
 if( is.null(bb) ){
  expbb=solve(A)%*%logit1( pi0 )
  expbb1=expbb[1]
  expbb[expbb<0]=1e-9
  bb=c(expbb1,log(expbb[-1]))
 }

 # bb1=c(bb[1],exp(bb[-1]))
 # pi=logistic1( A%*%exp(bb1) )
 # Print(expbb,bb,pi)

 # get pi as A%*%bb by native nlm
 res=nlm( get_ll, bb )
 bb=res$estimate
 bb1=bb[1]
 pi=logistic1( A%*%c(bb1,exp(bb[-1])) )
 ll=-res$minimum
 code=res$code

 if( print ){
  Print(r,n, pi0, pi00, pi, bb)
  Print(ll0,ll00,ll, code)
 }

 return( named_list(pi,bb,ll,ll0,code) )

} # end of get_monotone




comments('
np=5
A=matrix(0,np,np)
for( i in 1:np ){
 A[i,1:i]=1
}


r=c(1, 2, 4.05, 4, 5)
n=c(10, 10, 10, 10, 10)

resm <- get_monotone( r, n, A, print=1 )

')

np=5
A=matrix(0,np,np)
for( i in 1:np ){
 A[i,1:i]=1
}

RR="
  88.703782  17.727079
 161.885134 100.874207
  56.971263 162.736071
   3.118323 208.933854
   5.321499 193.728789
"; RR=cards(RR)

r=RR[,2]
n=rowSums(RR)


resm <- get_monotone( r, n, A, print=1 )



