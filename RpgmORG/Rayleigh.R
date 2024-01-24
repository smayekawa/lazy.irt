# Rayleigh quotient
#  q = t(x)%*%A%*%x /  t(x)%*%B%*%x
#  where B is psd.

np=5
a=1:np
A=a%*%t(a)

B=diag(np)

eigA=eigen(A)
evec=eigA$vector[,1]

res=optim( evec+1, function(x){ t(x)%*%A%*%x / t(x)%*%B%*%x }, method="BFGS" )


res2=optim( -evec, function(x){ -t(x)%*%A%*%x / t(x)%*%B%*%x }, method="BFGS" )


