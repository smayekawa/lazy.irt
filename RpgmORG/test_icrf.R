
thmin=-9
thmax=9
theta=seq(thmin,thmax,length=51)

paramVV=data.frame( name=11, type="B", ncat=2
     , p1=0.365273841812623, p2=0.561401148919932, p3=0, stringsAsFactors=0 )
# set_number=3
temp=irf(paramVV, theta=theta, plot=1, zero=1 )

temp=irf(paramVV, theta=theta, plot=1, zero=0 )

