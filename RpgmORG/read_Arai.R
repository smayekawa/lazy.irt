loc=sample(1:100000,10000)

theta=read.csv("RpgmORG/testD/pp.wge", header=0, sep="")[,2]
theta=theta[loc]
# hist(theta)

U=read.csv("RpgmORG/testD/pp1.csv", header=0)[,-1]
U=U[loc,]
X=rowSums(U)
hist(X, nclass=55)
