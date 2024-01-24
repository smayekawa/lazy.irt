
icrf <- function( paramj ){
 # to be used in JacobianMat
 # everything else is from outside
 paramjdf[1,3+(1:ncatj)]=paramj
 vpjhat=c( irf( paramjdf, theta, print=0 )$ICRF )
 return(vpjhat)
} # end of icrf


paramG=paramS1[3,]
paramP=paramS1[4,]

theta=seq(-3,3,length=21)

param=paramG

paramjdf=param
ncatj=paramjdf$ncat
Jac=dirf_p( param, theta, print=0, zero=1 )$Jack
JacN=JacobianMat( param[1,-(1:3)], icrf, ..eps.. = 1e-06 )
dif=Jac-JacN
Print(max(abs(dif)))

