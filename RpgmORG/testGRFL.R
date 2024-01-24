

pp=rbind(paramB1,paramB1,paramB1)
pp$name=1:nrow(pp)

out_obscore=obscore( pp, weight=NULL, npoints=121, print=0, plot=3 )

res=graded_info( out_obscore, ncat=3, method=1, plot=1 )
res=graded_info( out_obscore, ncat=5, method=1, plot=1 )
res=graded_info( out_obscore, ncat=9, method=1, plot=1 )



res2=flatten_SEM( param=pp, sigma=1, t_by=0.5, plot=1, print=1, npoints=151 )
res21=graded_info( res2$out_obscore, brk=res2$brk_x2u[,2], plot=1 )

res3=flatten_SEM( param=pp, sigma=0.5, t_by=0.5, plot=1, print=1, npoints=151 )
res31=graded_info( res2$out_obscore, brk=res3$brk_x2u[,2], plot=1 )

