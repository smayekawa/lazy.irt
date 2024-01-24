ts=cbind(res$t,res$s)
np=nrow(ts)
ts[np,2]=ceiling(ts[np,2])
Print(ts)
plot(ts)
maxs=ts[np,2]
Print(maxs)

midps=gen_midp( 0, maxs, maxs+1 ) # consecutive integers
midps=gen_midp( 0, maxs, 12)
#midps=round(midps)


midpt=interpol( ts[,2], ts[,1], midps )[,2]
brk=midp2brk( midpt )

Print( midps, round(midps), midpt, brk)

res2=graded_info( res$out_obscore, brk=brk, plot=1, scorey=round(midps) )

# res2=graded_info( res$out_obscore, brk=brk, plot=1, scorey=midpu )

# resx=graded_info( res$out_obscore, method=1, ncat=5, plot=1 )

