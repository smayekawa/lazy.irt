func0 <- function(x, y){ sum(sin(x))+y  }
res=grad(func0 , (0:10)*2*pi/10, y=10)



func2 <- function(x) c(sin(x), cos(x))
x <- (0:1)*2*pi
jacobian(func2, x)
jacobian(func2, x, "complex")


func3 <- function( x ){
 y1=x
 y2=x^2
}








f2 <- function( x ){
 f1 <- a[1]*x[1] + a[2]*x[2]^2
 f2 <- b[1]*x[1]^2 + b[2]*x[2]
 f3 <- a[1]*x[1]^2 + b[1]*x[2]^2
 return( c(f1,f2,f3) )
}
# analytic Jacobian
df2dx <- function( x ){
 df1dx1 <- a[1] ; df1dx2 <- 2*a[2]*x[2]
 df2dx1 <- 2*b[1]*x[1]; df2dx2 <- b[2]
 df3dx1 <- 2*a[1]*x[1]; df3dx2 <- 2*b[1]*x[2]
 Jac <- matrix(c(df1dx1, df2dx1, df3dx1, df1dx2, df2dx2, df3dx2), 3)
 return( Jac )
}
a=1:2; b=2:3; x=1:2
Jac=JacobianMat( x, f2 )
Print(df2dx(x), Jac)
