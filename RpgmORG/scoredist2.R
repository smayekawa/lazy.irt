ScoreDist <- function( U, p, w, omitmiss=0, print=0, plot=0 ){
 # calculate score dist from U, p and w
 # Shin-ichi Mayekawa
 # 20220119
 #

 # Args:
 #  U    binary response pattern matrix
 #  p    prob of each pattern (row of U)
 #  w    weight vector
 #  omitmiss = 1 to omit the zero freq cells
 #
 # Return
 #  ff  Freq dist of the weighted score as a matrix (score, freq)
 #

 score=U%*%w
 f=table2freqdist(xtabs(p~score))
 if( omitmiss == 0 ){
  ff=data.frame(score=(0:sum(w)))
  ff=merge( ff, f, by.x="score", by.y="score", all=TRUE)
  ff[is.na(ff)]=0
 }
 else ff=f

 ff=as.matrix(ff)

 tt=paste("prob=", paste(prob,collapse=", "), ",  w=", paste(w,collapse=", "))
 if( print ){
  cat(tt,"\n\n")
  Print(U,p,score, ff)
 }
 if( plot ){
  plot_freqdist(ff, title=tt)
 }

 return( ff )

} # end of ScoreDist




nitem=4
prob=rep(0.5,nitem)
pmat=rbind(prob,1-prob)
w1=rep(1,nitem)
w2=c(1,1,2,2)
w3=c(2,2,3,3)
w4=c(1,1,4,4)
w5=c(1,2,3,4)

# simulate sumsmnw

U=gen01pat( nitem, sort=1 )
p=exp( U%*%log(prob) + (1-U)%*%log(1-prob) )

res1=ScoreDist( U, p, w1 )
res2=ScoreDist( U, p, w2 )
res3=ScoreDist( U, p, w3 )
res4=ScoreDist( U, p, w4 )
res5=ScoreDist( U, p, w5 )
