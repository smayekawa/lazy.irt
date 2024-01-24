library(devtools)
paramS3=read.param( "d:/rpgm/packages/lazy.irt/data-raw/param_BGP5.dat")
weightS3=read.table( "d:/rpgm/packages/lazy.irt/data-raw/weight_BGP5.dat"
                     , header=1)
use_data(paramS3, overwrite=1)
use_data(weightS3, overwrite=1)
