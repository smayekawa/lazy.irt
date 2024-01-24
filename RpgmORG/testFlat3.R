

pp=rbind(paramB1,paramB1,paramB1)
pp$name=1:nrow(pp)


out_obscore2 <- obscore( pp, npoints=161, thmin=-4, thmax=4 )
res3=graded_info( out_obscore2, ncat=5, plot=1 )




res=flatten_SEM( out_obscore2, sigma=1, plot=1, print=1 )
res2c=graded_info( res$out_obscore, brk=res$brk_x2uc[,2]
                   , scorey=res$brk_x2uc[-1,1], plot=1 )


res2=graded_info( res$out_obscore, brk=res$brk_x2u[,2], plot=1 )


res31=graded_info( res$out_obscore, ncat=res$ncat, plot=1 )








plot(res3$tcc,res3$stdx_t)
plot(res3$TRFy_t,res3$stdy_t)




# creating Table{tab-Py-t}
out_obscore2 <- obscore( pp, npoints=21, thmin=-4, thmax=4 )
res3=graded_info( out_obscore2, ncat=5, plot=1 )

Print(res3$tcc)
Print(res3$tcc[5:17])
round(res3$TRFy_t[5:17],1)










