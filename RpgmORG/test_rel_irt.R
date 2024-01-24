param100=data.frame(name=paste("Q",1:100,sep=""), type="B"
                    ,ncat=2,p1=1,p2=0,p3=0,stringsAsFactors=0)
res <- rel_irt( param100, thmin=-4, thmax=4, npoint=121 )
