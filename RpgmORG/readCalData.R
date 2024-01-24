

infile="RpgmOLD/test11.dat"
paramCal2=read.table(file=infile, header=1, stringsAsFactors=0, skip=2 )
paramCal2=convPN2P(paramCal2)
use_data(paramCal1, overwrite=1)


infile="RpgmOLD/test2.dat"
paramCal1=read.table(file=infile, header=1, stringsAsFactors=0 )
use_data(paramCal1, overwrite=1)


Print(paramCal1,paramCal2)

