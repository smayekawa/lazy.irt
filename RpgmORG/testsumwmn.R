npoints=11; thmin=-4; thmax=4
theta=seq(thmin,thmax,length.out=npoints)
temp=irf( paramS2, theta, weightS21, print=0, debug=0, plot=0 )
icrf=temp$ICRF



Pk="
   Q1         Q2         Q3         Q4         Q5         Q6         Q7          Q8
c0 0.50000000 0.82200631 0.40000000 0.67642779 0.15446527 0.03229546 0.16960095  0.03090173
c1 0.50000000 0.17799369 0.60000000 0.32357221 0.69106947 0.46770454 0.66079810  0.46909827
c2         NA         NA         NA         NA 0.15446527 0.46770454 0.16960095  0.46909827
c3         NA         NA         NA         NA         NA 0.03229546         NA  0.03090173
"; Pk=cards(Pk, header=1)

v="
   r1 r2 r3 r4 r5 r6 r7 r8
v0  0  0  0  0  0  0  0  0
v1  1  1  1  1  1  1  1  1
v2 .  .  .  .   2  2  2  2
v3 .  .  .  .   .  3  .  3
"; v=t( cards(v, header=1) )
w=c(1,2,1,2,1,2,1,2)

Print(Pk,v,w)

temp=sumsmnw( Pk, t(v), w, compress=0, ncat=c(2,2,2,2,3,4,3,4)
         , print=2, plot=0, debug=0 )
