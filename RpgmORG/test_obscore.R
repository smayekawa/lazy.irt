cards <- '
name type ncat p1 p2 p3
Q1 B 2 0.5 0 0
Q2 B 2 1.0 0 0
Q3 B 2 1.5 0 0
'
param <- read.table( text=cards, header=1, stringsAsFactors=0)

cards2 <- '
name type ncat w v0 v1
Q1 B 2   1   0 1
Q2 B 2   2   0 1
Q3 B 2   4   0 1
'
weight <- read.table( text=cards2, header=1, stringsAsFactors=0)


# temp=irf( param, plot=1, zero=0 )


cards <- '
name type ncat p1 p2 p3
Q1 B 2 0.5 0 0
Q2 B 2 0.5 0 0
Q3 B 2 1.0 0 0
Q4 B 2 1.0 0 0
'
param <- read.table( text=cards, header=1, stringsAsFactors=0)

cards2 <- '
name type ncat w v0 v1
Q1 B 2   1   0 1
Q2 B 2   1   0 1
Q3 B 2   2   0 1
Q3 B 2   2   0 1
'
weight2 <- read.table( text=cards2, header=1, stringsAsFactors=0)


cards3 <- '
name type ncat w v0 v1
Q1 B 2   1   0 1
Q2 B 2   2   0 1
Q3 B 2   4   0 1
Q3 B 2   8   0 1
'
weight3 <- read.table( text=cards3, header=1, stringsAsFactors=0)
weight3$w <- rev(weight3$w)

# temp=irf( param, plot=1, zero=0 )


res2=obscore( param, weight=weight2, plot=1 )
res3=obscore( param, weight=weight3, plot=1 )

Print(as.matrix(res2$obs_stat), fmt="8.3")
Print(as.matrix(res3$obs_stat), fmt="8.3")


