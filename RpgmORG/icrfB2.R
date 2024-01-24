icrfB <- function( param, theta, maxZ=700, zero=0, print=0, plot=0 ){
 # calculation of Binary Logistic ICRF
 # Shin-ichi Mayekawa
 # 120201
 # maxZ: 120202
 # zero: 120209
 # renamed as irfB -> icrfB: 120210
 # dataframe input: 120213,14
 # checkparam: 120224,29
 # when nrow(param) == 1: 120308
 # bugfix: 120908
 # when c-parameter is NA: 120921
 #
 # Args:
 #  param    nitems x 3 item parameters matrix or data frame
 #  theta    npoints x 1 theta values
 #  zero     = 1 to include 1-P as the 0-th category probability
 #
 # Values:
 #   P       length(theta) x nitems   ICRF matrix when zero=0
 #           length(theta) x 2*nitems ICRF matrix when zero=1
 #
 #
 # Needs:
 #   checkparam
 #
 
  # argument name
  paramname=as.character(substitute(param))
  isdf=0
  if( is.data.frame(param) ){
   # param and weight given as data frame
   isdf=1
   param=checkparam( param, "B", "icrfB" )
   if( is.null(param) ){
    cat("\n\nInput Item Parameter ERROR.\n\n")
    return()
   }
  }

  # const
  nitems=nrow(param)
  iname=rownames(param)
  if( length(which(colnames(param) == "ncat")) > 0 )
   param=param[,-which(colnames(param) == "ncat"),drop=F]
  a=param[,1]; b=param[,2]; c=param[,3]

  # if the c-parameter is NA   120921
  c[is.na(c)]=0

  if( is.matrix(theta) ) theta=as.vector(theta)
  npoints=length(theta)
  
  if( length(a) > 1 ) da=diag(a)
  else da=matrix(a,1,1)
  Z=-1.7*outer(theta,b,"-")%*%da
  Z[which(Z > maxZ)]=maxZ
  Z[which(Z < -maxZ)]=-maxZ
  cc=matrix(c,npoints,nitems,byrow=1)
  P=cc+(1-cc)/(1+exp(Z))
  colnames(P)=rownames(param)
  if( zero == 1 ){
   P0=matrix(0,npoints,2*nitems)
   P0[,2*(1:nitems)]=P
   P0[,2*(1:nitems)-1]=1-P
   catname=outer(paste(iname,"_",sep=""),0:1,paste,sep="")
   colnames(P0)=matrix(t(catname),1)
   P=P0
  }
  rownames(P)=format(theta,digits=3)
  if( print > 0 ){
   cat("\nIRF of Binary Logistic Items \n ")
   if( isdf )
    cat(" item parameter data frame name =", paramname,"\n")
   else
    cat(" item parameter matrix name =", paramname,"\n")
   cat("  # of item parameters =",nitems,
       ",  # of theta points =", npoints,
       ",  zero category =",zero, "\n")
   Print(param)
   Print(P, digits=3)
  }
  if( plot > 0 ){
   matplot(theta,P, type = "l", ylab="IRF"
         , xlim=c(min(theta),max(theta)), ylim=c(0,1)
         , main="Plot of ICRF of Binary Items")
  }
  return( P )
  
} # end of icrfB




theta=seq(-3,3,length.out=11)
icrfB(pp,theta, plot=1, print=2, zero=1)


comments('
theta=seq(-3,3,length.out=11)
P=icrfB(param,theta,print=1, plot=1, zero=0)
')


comments('
paramB="
    a   b    c
Q1  1    -1    0
Q2  1     0    0
Q3  1     1    0
Q4  1.5  -1    0
Q5  1.5   0    0
Q6  1.5   1    0
Q7  1    -1    0.2
Q8  1     0    0.2
Q9  1     1    0.2
"; paramB=lazy.matrix(paramB,header=1)

theta=seq(-3,3,length.out=11)
P=icrfB(paramB,theta,print=1,plot=0, zero=1)
')



