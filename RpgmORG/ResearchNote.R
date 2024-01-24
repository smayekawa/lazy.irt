# DNC Research Note 17-01

thmin=-4
thmax=4

# thmin=-3; thmax=3 # results in 0-24 graded score, not 0-25. 20170530dnc

# parameter data frame
pp=rbind(paramB1,paramB1,paramB1)
pp$name=1:nrow(pp)

out_obscore2 <- obscore( pp, npoints=161, thmin=thmin, thmax=thmax )
res3=graded_info( out_obscore2, ncat=5, plot=1 )

# flatten
res=flatten_SEM( out_obscore2, sigma=1, plot=1, print=1 )
res2c=graded_info( res$out_obscore, brk=res$brk_x2uc[,2]
                   , scorey=res$brk_x2uc[-1,1], plot=1 )

res2=graded_info( res$out_obscore, brk=res$brk_x2u[,2], plot=1 )

# simple grading with 26 categories
res31=graded_info( res$out_obscore, ncat=res$ncat, plot=1 )


# Grading of S

# t: true score of X,  s: true score of Y
ts=cbind(res$t,res$s)
np=nrow(ts)
ts[np,2]=ceiling(ts[np,2])
Print(ts)
plot(ts)
maxs=ts[np,2]
Print(maxs)

# new midpoints of S
midps=gen_midp( 0, maxs, 19)
# correspoinding midpoints of t and associated breaks
midpt=interpol( ts[,2], ts[,1], midps )[,2]
brk=midp2brk( midpt )

# Print( midps, round(midps), midpt, brk)
res2=graded_info( res$out_obscore, brk=brk, plot=1, scorey=round(midps) )


# new midpoints of S
midps=gen_midp( 0, maxs, 16)
# correspoinding midpoints of t and associated breaks
midpt=interpol( ts[,2], ts[,1], midps )[,2]
brk=midp2brk( midpt )

# Print( midps, round(midps), midpt, brk)
res2=graded_info( res$out_obscore, brk=brk, plot=1, scorey=round(midps) )



# new midpoints of S
midps=gen_midp( 0, maxs, 12)
# correspoinding midpoints of t and associated breaks
midpt=interpol( ts[,2], ts[,1], midps )[,2]
brk=midp2brk( midpt )

# Print( midps, round(midps), midpt, brk)
res2=graded_info( res$out_obscore, brk=brk, plot=1, scorey=round(midps) )



