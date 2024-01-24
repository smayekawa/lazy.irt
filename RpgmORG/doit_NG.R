# must follow doit.R
# 
doit_NG <- function( reprange, interpol_method="constant", HKmethod=1
                  ,  print=1, plot=1  ){
 
 for( rep in reprange ){
  
  from_named_list( SampleFreqs[[rep]] )
  
  res=rpm( freq4, interpol_method=interpol_method, HKmethod=HKmethod
           , title=rep, print=print, plot=plot )
  
 } # end of repoop
 

} # end of doit_NG




graph2pdf( "RpgmORG/doit_NG.pdf" )



reprange=check3_NG$rep


interpol_method="linear"; HKmethod=1
doit_NG( reprange, interpol_method=interpol_method, HKmethod=HKmethod )


graph2pdf( close=1 )