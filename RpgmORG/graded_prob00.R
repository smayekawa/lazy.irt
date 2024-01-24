
graded_prob <- function( x, mean, std, min=mean-5*std, max=mean+5*std
                         , truncate=1, method=1){
 # calculate prob for discrete x
 # Shin-ichi Mayekawa
 # 20180629
 #

 #
 # method = 1 to use the difference of pnorms
 #        = 2 to use nornalized dnorm at x
 #

 if( method == 1 ){
  brk=midp2brk( x, min=min, max=max )
  cp=pnorm( brk, mean, std )
  pp=cumsuminv(cp)
  p=pp[-1]
  p[p<0]=0; p[p>1]=1
  p=p/sum(p)
 }
 else if( method == 2 ){
  p=dnorm( x, mean, std )
  p=p/sum(p)
 }

 # Print( x, brk, cp, pp, p )
 # Print(sum(p))

 # Print( x, p, p2, p-p2 )

 return( p )

} # end of graded_prob




graded_prob( 0:10,  5,5, 0,10 )


graded_prob( 0:3,  1.5,1, 0,3 )

graded_prob( 0:3,  2.9,1, 0,3 )
graded_prob( 0:3,  3,1, 0,3 )
