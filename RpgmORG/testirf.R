
temp=c("
 name type ncat p1 p2 p3 p4 p5 p6
       Q1 G 6  1, -3,-1,0,1,4
       Q2 P 6  1, -3,-1,0,1,4
       ")
param=cards(temp,header=1)

temp=c("
 name type ncat p1 p2 p3
 Q1 G 3  1, 0,1
 Q2 P 3  1, 0,1
   ")
param=cards(temp,header=1)

irf2=irf( param, print=0, plot=1, npoints=21 )

pPfromG=fitP2G(param[1,], plot=1)$paramP
pGfromP=fitG2P(pPfromG, plot=1)

pGfromP=fitG2P(param[2,], plot=1)
pPfromG=fitP2G(pGfromP, plot=1)$paramP



temp=c("
 name type ncat p1 p2 p3 p4
 Q1 G 4  1,  0,1,4
 Q2 P 4  1,  0,1,4
   ")
param=cards(temp,header=1)

irf2=irf( param, print=0, plot=1, npoints=21 )

pPfromG=fitP2G(param[1,], plot=1)$paramP
pGfromP=fitG2P(pPfromG, plot=1)

pGfromP=fitG2P(param[2,], plot=1)
pPfromG=fitP2G(pGfromP, plot=1)$paramP
















