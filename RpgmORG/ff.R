
library(codetools)

funclist=ls(package:lazy.irt)
res=sapply( funclist, find_functions, packages="lazy.irt" )


funclist=ls("package:lazy.irt")
res <- sort( unique( unlist( lapply( funclist, find_functions ) ) ) )


res <- unique( unlist( lapply(
 funclist, find_functions, packages=c("lazy.mat","lazy.tools","lazy.irt") )) )

res1 <- unique( unlist( lapply(
 funclist, find_functions, packages=c("lazy.irt","lazy.tools") ) ) )

res2 <- unique( unlist( lapply(
 funclist, find_functions, packages=c("lazy.irt","lazy.mat") ) ) )


res11 <- unique( unlist( lapply( res1, find_functions, packages=NULL ) ) )
res22 <- unique( unlist( lapply( res2, find_functions, packages=NULL ) ) )

res111 <- unique( unlist( lapply( res11, find_functions, packages=NULL ) ) )
res222 <- unique( unlist( lapply( res22, find_functions, packages=NULL ) ) )

find_functions(trim)
