
graded_range <- function( conf_level=c(68,90,95), wid=1:2
                          , rho=c(0.7,0.8,0.9,0.95), rangeX=2, alpha=NULL ){
 # calculate the range of graded score
 # Shin-ichi Mayekawa
 # 20170409
 # rangeX: 20170530dnc,0609dnc
 #
 # conf_level  confidence level of the interval of T: 100*(1-alpha)
 # wid         width of the interval:  make this equal to 1 or 2.
 # rho         reliability coefficient
 # alpha       alpha to calculate conf_level
 # rangeX      range of X measured by sigma_X
 #                where the distribution concentrated: 2 or larger.
 #             This is denoted as n in RN1701.
 #    Prob( mu_X - n*sigma_X < X < mu_X + n*sigma_X ) === 1-alpha
 #    The length of above interval is equal to 2*n*sigma_X.
 #    Solve 2*n = z_{1-alpha/2}
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

    r1=round( rangeX*i/(za2*sqrt(1-j)) )
    r2=round( rangeX*i/(za2*sqrt((1-j)*j)) )

#    cat("rangeX=", rangeX, ",  CL =", conf_level[k]
#        , ", wid =",i, ", rho = ",j, ",  r1 =", r1, ", r2 =", r2, "\n")

    res=rbind(res,cbind(rangeX=rangeX,CL=conf_level[k],wid=i,rho=j,r1=r1,r2=r2))
   }
  }
 }

 res=as.data.frame(res)

# print(res)

 return(res)


} # end of graded_range



res=graded_range( rangeX=2)


res=graded_range( rangeX=3)


# graded_range( conf_level=68, rho=c(0.6,0.7,0.8,0.9) )

# graded_range( conf_level=90, rho=c(0.6,0.7,0.8,0.9) )

