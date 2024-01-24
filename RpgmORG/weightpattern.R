
from=c(1,2,3,4)
nitem=4
total=10
res=NULL
for( i in 1:10000 ){
 x=sample( from, nitem, replace=1 )
 if( sum(x) == total ){
  res=rbind(res,sort(x))
 }
}
res=unique(res)
Print(res)





from=c(2,4,5)
nitem=25
total=100
res=NULL
for( i in 1:100000 ){
 x=sample( from, nitem, replace=1 )
 if( sum(x) == total ){
  res=rbind(res,sort(x))
 }
}
res=unique(res)
Print(res)







