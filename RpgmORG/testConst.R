# classical test theory
# 20170604,05
#
# T is the true score
# E is the error score
# s = E(E|T) = the conditional expectation of E given T.
#
# It is assumed that E(s)=0, and E(T s)=0
#
# The expectation is calculated assuming discrete uniform distribution:
#   E(T) = sum(T_i) / n
#   E(s) = sum(s_i) / n
#   E(T s) = sum( T_i* s_i) / n
#
# Therefore,  E(s)=0, and E(T s)=0 means:
#   sum(s_i)=0  and sum( T_i* s_i)=0
# or
#   B s = 0
# where
#   B = 1,   1,   ..., 1
#       T_1, T_2, ..., T_n
#
# Convert B to A form and calculate s=A gamma
# where gamma is free.
#


n=11
one=rep(1,n)/n
T=(1:n)-1
B=rbind(one,T)
res=BtoA( B, rep(0,nrow(B)) )
A=res$A
Print(B,res$A, fmt="7.4")

gamma=c(rep(0,ncol(A)-1),1)
s=A%*%gamma
Print(s, B%*%s)

plot( T,s )
Print(mean(s),var(T,s))


ss=(T-mean(T))^2
gamma=solve(t(A)%*%A)%*%t(A)%*%ss
s=A%*%gamma

X=T+s/n

plot( T,s )
Print(mean(s),var(T,s))
Print(X,T,s, B%*%s)
plot( T,X )
















