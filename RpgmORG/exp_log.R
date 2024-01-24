
# max arg to exp is 709.782:
#  x  =   709.782 ,   exp(x)  =   1.796412e+308

for( x in seq(708,710,.1) ){
 Print( x, exp(x) )
}

# min arg to exp is -745
#  x  =   -745 ,   exp(x)  =   4.940656e-324

for( x in seq(720,750,.1) ){
 Print( -x, exp(-x) )
}



# min arg to log is 1e-308:
#  x  =   308 ,   log(1/10^x)  =   -709.1962

for( x in seq(305,310,0.1) ){
 Print( x, log(1/10^x) )
}


# max arg to log is 10e308:
#  x  =   308 ,   log(10^x)  =   709.1962

for( x in seq(305,310,0.1) ){
 Print( x, log(10^x) )
}



# max arg to logit( logistic(x) ) is 36:
#   x  =   36 ,   logistic  =   1 ,   logit  =   36.04365

# min arg to logit( logistic(x) ) is -709:
#   x  =   -709 ,   logistic  =   1.216781e-308 ,   logit  =   -709

for( x in 30:40 ){
 logistic=1/(1+exp(-x))
 logit=log(logistic/(1-logistic))
 Print(x,logistic,logit)
}

for( x in -710:-700 ){
 logistic=1/(1+exp(-x))
 logit=log(logistic/(1-logistic))
 Print(x,logistic,logit)
}



# min arg to 1/1e-x is 308
# x  =   -308 ,   p  =   1e-308 ,   rp  =   1e+308

for( x in -300:-330 ){
 p=10^x
 rp=1/p
 logit=log(p/(1-p))
 Print(x,p,rp,logit)
}



# min arg to logistic( logit(1e-x) ) is 323:
#  x  = -323 ,   p = 9.881313e-324 , logistic = 0, logit = -743.7469

# max arg to logistic( logit(1-1e-x) ) is 16.25:
#   x  =   -16.25 ,   p  =   1 ,   logistic  =   1 ,   logit  =   36.7368


for( x in seq(-16,-17,-.01) ){
 p=1-10^x
 logit=log(p/(1-p))
 logistic=1/(1+exp(-logit))
 Print(x,p,logistic,logit)
}





