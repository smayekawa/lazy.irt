# multi-group analysis with different theta distributions
set.seed(1701)
indata1 <- gendataIRT( 1, paramS3, npoints=5000, thdist="rnorm"
                       , compress=1 )$U
indata1 <- as.data.frame(indata1,row.names=NULL)
indata1 <- data.frame(group="G1",indata1, stringsAsFactors=0,row.names=NULL)
indata2 <- gendataIRT( 1, paramS3, npoints=5000, thdist="rnorm", compress=1
                       , thmean=1, thstd=1 )$U
indata2 <- as.data.frame(indata2,row.names=NULL)
indata2 <- data.frame(group="G2",indata2, stringsAsFactors=0,row.names=NULL)
indata3 <- gendataIRT( 1, paramS3, npoints=5000, thdist="rnorm", compress=1
                       , thmean=-0.5, thstd=0.7 )$U
indata3 <- as.data.frame(indata3,row.names=NULL)
indata3 <- data.frame(group="G3",indata3, stringsAsFactors=0,row.names=NULL)
indata123 <- rbind(indata1,indata2,indata3)
itemtype <- paramS3$type
# This will not converge: increase maxiter.
res1 <- uIRT( indata123, "group", type=itemtype, maxiter=1000, baseform=1
              , estmu=1, estsigma=1, maxalpha=-1, SQUAREM=3, plot=1, always=1 )




paramL3=rbind(paramS3,paramS3,paramS3,paramS3,paramS3,paramS2,paramS2,paramS2)
set.seed(1701)
indata1 <- gendataIRT( 1, paramL3, npoints=5000, thdist="rnorm"
                       , compress=1 )$U
indata1 <- as.data.frame(indata1,row.names=NULL)
indata1 <- data.frame(group="G1",indata1, stringsAsFactors=0,row.names=NULL)
indata2 <- gendataIRT( 1, paramL3, npoints=5000, thdist="rnorm", compress=1
                       , thmean=1, thstd=1 )$U
indata2 <- as.data.frame(indata2,row.names=NULL)
indata2 <- data.frame(group="G2",indata2, stringsAsFactors=0,row.names=NULL)
indata3 <- gendataIRT( 1, paramL3, npoints=5000, thdist="rnorm", compress=1
                       , thmean=-0.5, thstd=0.7 )$U
indata3 <- as.data.frame(indata3,row.names=NULL)
indata3 <- data.frame(group="G3",indata3, stringsAsFactors=0,row.names=NULL)
indata123 <- rbind(indata1,indata2,indata3)

itemtype <- paramL3$type
# This will not converge: increase maxiter.
res1 <- uIRT( indata123, "group", type=itemtype, maxiter=1000, baseform=1
              , estmu=1, estsigma=1, maxalpha=-1, SQUAREM=0, plot=1, always=1 )













