



# tiny set of binary items
param=paramB1
maxscore=sum(param$ncat-1)
# param$p1=1
out_obscore <- obscore( param, npoints=161, thmin=-4, thmax=4 )

res=flatten_SEM( out_obscore, sigma=1, plot=1, print=1 )
res2c=graded_info( res$out_obscore, brk=res$brk_x2uc[,2]
                  , scorey=res$brk_x2uc[-1,1], plot=1 )
res2=graded_info( res$out_obscore, brk=res$brk_x2u[,2], plot=1 )
# plot(res$brk_x2u)

res31=graded_info( res$out_obscore, ncat=res$ncat, plot=1 )
res32=graded_info( res$out_obscore, ncat=res$ncat, plot=1, method=2 )







# out_obscore <- obscore( pp, npoints=291, thmin=-5, thmax=5 )

pp=rbind(paramB1,paramB1,paramB1)
pp$name=1:nrow(pp)


paramB11=paramB1
paramB11$p2=paramB11$p2/2
# paramB11$p1=paramB11$p1/1.4
pp=rbind(paramB11,paramB11,paramB11)
pp$name=1:nrow(pp)

comments('

 pp=paramB2
 pp$p3=0

')


out_obscore2 <- obscore( pp, npoints=21, thmin=-4, thmax=4 )
res3=graded_info( out_obscore2, ncat=3, plot=1 )
res3=graded_info( out_obscore2, ncat=5, plot=1 )

res=flatten_SEM( out_obscore2, sigma=1, plot=1, print=1 )
res2c=graded_info( res$out_obscore, brk=res$brk_x2uc[,2]
                  , scorey=res$brk_x2uc[-1,1], plot=1 )
res2=graded_info( res$out_obscore, brk=res$brk_x2u[,2], plot=1 )

res31=graded_info( res$out_obscore, ncat=res$ncat, plot=1 )
res32=graded_info( res$out_obscore, ncat=res$ncat, plot=1, method=2 )


