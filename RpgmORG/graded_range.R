
graded_range <- function( conf_level=c(68,90,95), wid=1:2, rho=c(0.7,0.8,0.9)
                          , alpha=NULL ){
 # calculate the range of graded score
 # Shin-ichi Mayekawa
 # 20170409
 #

 nwid=length(wid)
 nrho=length(rho)

 if( !is.null(alpha) ) conf_level=100*(1-alpha)
 else{
  if( any(conf_level < 1) ) conf_level=100*conf_level
  alpha=(100-conf_level)/100
 }
 ncl=length(conf_level)

 res=NULL
 for( k in 1:ncl ){
  za2=qnorm(1-alpha[k]/2)

  for( i in wid ){
   for( j in rho ){

    #Print(k,i,j)

    r1=round( 2*i/(za2*sqrt(1-j)) )
    r2=round( 2*i/(za2*sqrt((1-j)*j)) )

    cat("CL =", conf_level[k]
        , ", wid =",i, ", rho = ",j, ",  r1 =", r1, ", r2 =", r2, "\n")

    res=rbind(res,cbind(CL=conf_level[k],wid=i,rho=j,r1=r1,r2=r2))
   }
  }
 }

 res=as.data.frame(res)

 return(res)


} # end of graded_range



res=graded_range()

comments('
graded_range( conf_level=68, rho=c(0.6,0.7,0.8,0.9) )

graded_range( conf_level=90, rho=c(0.6,0.7,0.8,0.9) )
')
