ScoreDist <- function( U, prob, w, omitmiss=0, print=0, plot=0 ){
 # calculate score dist from U, p and w
 # Shin-ichi Mayekawa
 # 20220119
 #

 # Args:
 #  U        binary response pattern matrix (2^nitems x nitems)
 #  prob     prob of correct answer for each item
 #  w        item weight vector
 #  omitmiss = 1 to omit the zero freq cells
 #
 # Return
 #  ff  Freq dist of the weighted score as a matrix (score, freq)
 #

 p=exp( U%*%log(prob) + (1-U)%*%log(1-prob) )
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
  cat("\n\n", tt,"\n\n")
  Print(U,p,score, ff)
 }
 if( plot ){
  plot_freqdist(ff, title=tt)
 }

 return( ff )

} # end of ScoreDist




# irf
P=matrix(c(
   0.5, 0.7, 0.9
 , 0.3, 0.5, 0.7
 , 0.3, 0.5, 0.7
 , 0.1, 0.3, 0.5
), 4,3, byrow=1 )
dimnames(P)=list(paste("Q",1:nrow(P),sep=""),c("L","M","H"))
H=c(0.2, 0.6, 0.2)

nitem=nrow(P)
nclass=ncol(P)
cname=colnames(P)


# item weight
w1=rep(1,nitem)
w2=c(1,1,2,2)
w3=c(2,2,3,3)
w4=c(1,1,4,4)
w5=c(1,2,3,4)

# All possible item response patterns and the associated probs
U=gen01pat( nitem, sort=1 )

w=rev(w1)

score=0:sum(w)
cprob=matrix(0,length(score),nclass)
rownames(cprob)=score
colnames(cprob)=cname
for( q in 1:nclass){
 pq=P[,q]

  res=ScoreDist( U, pq, w, print=1, plot=1 )
  cprob[,q]=res[,2]

}


mprob=cprob%*%H

Print(cprob, mprob)

plot_freqdist( score, mprob )




comments(
'



'
)



