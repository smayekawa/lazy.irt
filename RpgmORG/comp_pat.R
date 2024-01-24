nitem=4
pat <- gen01pat( nitem, sort=1 )
pat%*%weight3$w

pat_s=cbind(pat, rowSums(pat), pat%*%weight2$w, pat%*%weight3$w)
colnames(pat_s)=c( paste("Q",1:nitem,sep=""), "sum", "wsum2", "wsum3" )
Print(pat_s)

