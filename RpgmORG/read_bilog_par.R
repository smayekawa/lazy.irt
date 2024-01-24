
library(lazy.irt)


read_blg_par <- function( infile ){
 res=read.table( infile, skip=4, fill=TRUE, stringsAsFactors=0
                , col.names=c("name","subtest","intercept","intercept_se"
                              , "a", "a_se","b","b_se","disp","disp_se"
                              , "c", "c_se", "DRIFT", "DRIFT_se"
                              , "location"))
 return( res )
} # end of read_blg_par




# read bilog .par file
res=read_blg_par( "RpgmORG/Hashimoto/C12302PL.PAR" )

# divide it into two.
locE=which(substr(res$name,1,1) == "E")
locP=which(substr(res$name,1,1) == "P")

# Ethics
paramE=data.frame( name=res[locE,"name"], type="B", ncat=2
                   , res[locE,c("a","b","c")], stringsAsFactors=0 )
colnames(paramE)[4:6]=c("p1","p2","p3")
weightE=create_weight_df( paramE )$weight
w=c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,3,4,3,3,3,3,4)
weightE$w=w

# Polieco
paramP=data.frame( name=res[locP,"name"], type="B", ncat=2
                   , res[locP,c("a","b","c")], stringsAsFactors=0 )
colnames(paramP)[4:6]=c("p1","p2","p3")
weightP=create_weight_df( paramP )$weight
w=c(3,4,3,4,3,3,4,3,3,4,3,3,4,3,3,3,3,3,4,3,3,3,4,3,3,3,4,4,3,2,2)
weightP$w=w


# true score eauating
res1=tseq( paramP, paramE, weight1=weightP, weight2=weightE, print=2, plot=1
      , round=1 )

# observed score
res2=oseq( paramP, paramE, weight1=weightP, weight2=weightE, print=2, plot=1
      , round=1 )

# comparison
matplot( res1$x1, cbind(res1$x2_1, res2$x2_1), type="l", col=c("black","red")
         , main="comparison of tseq(black) vs oseq(red)")


