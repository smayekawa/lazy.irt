
resff=find_functions( func=ls("package:lazy.irt"), unlist=0, nonative=0 )

res_P=find_calling_functions( c("icrfP","icrfPN","icrfPN0","dicrfP","dicrfPN","dicrfPN0"), resff )

res_B=find_calling_functions( c("icrfB","dicrfB"), resff )$res

res_G=find_calling_functions( c("icrfG","dicrfG"), resff )$res


res_plot=find_calling_functions( c("plot","matplot"), resff )$res

