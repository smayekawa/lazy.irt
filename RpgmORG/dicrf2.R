
dirf2 <- function( param, theta, print=1, plot=0 ){

 dicrf <- function( theta ){
  res=vec( dirf( param, theta, print=0 )$dICRF )
  return( res )
 }

 d2=JacobianMat( theta, dicrf )

 Print(d2, sum(d2))

} # end of dirf2


theta=seq(-4,4,length=21)

dirf2( paramA1[8,], -1 )
dirf2( paramA1[8,], 0 )
dirf2( paramA1[8,], 1 )





