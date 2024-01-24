
get_monotone <- function( r, n, A, bb=NULL ){
 # find monotone pi which maximize the binominal likelihood
 # Shin-ichi Mayekawa
 # 20170621
 #

 get_ll <- function( bb ){
  # binominal log lihelihood for nlm
  pi=A%*%exp(bb)
  ll=sum(r*log(pi) + (n-r)*log(1-pi))
  return( -ll )
 } # end of get_ll

 # # of elements
 np=length(c(r))

 # unrestricted mle of pi
 pi0=r/n
 ll0=sum(r*log(pi0) + (n-r)*log(1-pi0))

 pi00=isoreg(pi0)$yf
 ll00=sum(r*log(pi00) + (n-r)*log(1-pi00))

 # initial value of bb
 if( is.null(bb) ){
  expbb=solve(A)%*%pi0
  expbb[expbb<0]=1e-2
  bb=log(expbb)
 }

 # get pi as A%*%bb by native nlm
 res=nlm( get_ll, bb )
 bb=res$estimate
 pi=A%*%exp(bb)
 ll=-res$minimum
 code=res$code

 Print(r,n, pi0, pi00, pi, bb)
 Print(ll0,ll00,ll, code)

 return( pi )

} # end of get_monotone



np=5
A=matrix(0,np,np)
for( i in 1:np ){
 A[i,1:i]=1
}


r=c(1,2,4.6,4,5)
n=c(10,10,9,10,10)

resm <- get_monotone( r, n, A )





