a=1:4
A=a %*% t(a)
ev=eigen(A)
eva=diag(ev$values)
eve=ev$vectors

Afromev=eve%*%eva%*%t(eve)

d=A-Afromev
