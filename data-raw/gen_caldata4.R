

set.seed(234)

eps=.3

skip=0;
infile="RpgmOLD/param_BGP5.dat"

infile="RpgmOLD/param_BGP52.dat"

infile="RpgmOLD/param_BGP53.dat"

infile="RpgmOLD/param_BGP54.dat"

infile="RpgmOLD/param_BGPN1.dat"

infile="RpgmOLD/param_BGP55.dat"

infile="RpgmOLD/param_B1.dat" ; skip=1

paramC=read.param( infile, print=1, skip=skip )
paramC=cbind(gfid=rep("form1",nrow(paramC)),paramC)
paramC$gfid=as.character(paramC$gfid)
paramC$name=as.character(paramC$name)



cn=colnames(paramC)
nitems=nrow(paramC)
nc=ncol(paramC)
locB=grep("^B[[:digit:]]*$", paramC$type)
locB2=which(paramC$type == "B2")
locG=which(paramC$type=="G")
locPN=which(paramC$type=="PN")
loca=which(cn == "p1")
Print(locB,locB2,locG,locPN,loca)
paramp3=paramC$p3

even=seq(2,nitems,2)
odd=seq(1,nitems,2)


paramC1=paramC[-c(3,11,20),]


paramC2=paramC
paramC2$gfid="form2"
paramC2[,loca]=paramC[,loca]*.75
paramC2[,(loca+1):nc]=paramC[,(loca+1):nc]/.75  -1/.75
paramC2[,loca:nc]=paramC2[,loca:nc] + eps*runif(nrow(paramC2)*(nc-loca+1))
paramC2[locB,loca+2]=paramC[locB,loca+2] + eps*runif(length(locB))
paramC2[locB2,loca+2]=paramp3[locB2]
paramC2=paramC2[c(odd,2,6),]

paramC3=paramC
paramC3$gfid="form3"
paramC3[,loca]=paramC[,loca]/.75
paramC3[,(loca+1):nc]=paramC[,(loca+1):nc]*.75  +.75
paramC3[,loca:nc]=paramC3[,loca:nc] + eps*runif(nrow(paramC3)*(nc-loca+1))
paramC3[locB,loca+2]=paramC[locB,loca+2]+ eps*runif(length(locB))
paramC3[locB2,loca+2]=paramp3[locB2]
paramC3=paramC3[c(even,1),]


Print(paramC1,paramC2,paramC3)

res=rbind(paramC1,paramC2,paramC3)

paramCal1=res
rownames(paramCal1)=NULL

paramCal2=convPN2P(paramCal1)

Print(paramCal1, paramCal2, digits=2)


# read.cal.Arai.R
# res=cala.core( indata=paramCal1, baseform=1, maxiter=1000, eps=1e-8, print=3 )
