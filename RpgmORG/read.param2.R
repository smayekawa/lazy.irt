read.param <- function( infile, skip=0, nrows=-1, na.strings="NA"
                      , sep="", print=0 ){
 # reading IRT parameter files
 # Shin-ichi Mayekawa
 # 120208,09,12
 # clean input data frame with NA: 120224
 # use of col.names: 120224
 # a dropped: 120225
 # sep: 120228,29
 # type Bxx items: 120919,21
 # simplified version: 121130
 #
 #
 # It is assumed that the data file has the following quantities
 #   in the following order.
 #
 #   item_name, item_type, ncat, p1, p2, ..., p[ncat+1]
 #
 #
 #   When item_type == "B", ( ncat == 2 )
 #     p1 is the discrimination parameter
 #     p2 is the difficulty parameter
 #     p3 is the gussing parameter
 #
 #   When item_type == "G",
 #     p1 is the discrimination parameter
 #     p2 is the category threshold parameter b_{j1}
 #     p3 is the category threshold parameter b_{j2}
 #
 #     p_ncat is the category threshold parameter b_{j,[ncat-1]}
 #
 #   When item_type == "P",
 #     p1 is the discrimination parameter
 #     p2 is the step parameter b_{j1}
 #     p3 is the step parameter b_{j2}
 #
 #     p_ncat is the step parameter b_{j,[ncat-1]}
 #
 #   When item_type == "PN,
 #     p1 is the slope parameter
 #     p2 is the intercept parameter c_{j1}
 #     p3 is the intercept parameter c_{j2}
 #
 #     p_ncat is the intercept parameter c_{j,[ncat-1]}
 #
 #   When item_type == "N,
 #     p1 - p_(ncat[j]-1) are the slope parameters
 #     p_ncat[j] - p_2*(ncat[j]-1) are the intercept parameters
 #
 #   Regardless of the item_type, there will be ncat[j]
 #   item paramters for item j, except for the Binary Items
 #   which have the gusseing param as the ncat[j]+1st,
 #   and the Nominal Items which have 2*(ncat[j]-1) item parameters.
 #
 #
 # Value: as data.frame
 #
 #   name type ncat  p1  p2  p3   .....
 #

 param=read.table( infile, header=1, fill=1
                 , na.strings=na.strings, skip=skip, flush=1, sep=sep
                 , stringsAsFactors=0
                 )

 locnotna=which( !apply( param, 2, function(x) all(is.na(x)) ) )
 param=param[,locnotna]
 cname=colnames(param)

 # in case the header has a variable 'name'.
 if( cname[1] == "name" ){
  rownames(param)=param[,1]
   param=param[,-1]
 }

 # param names from the row names
 param=cbind(name=rownames(param),param)
 cname=colnames(param)

 nitems=nrow(param)
 ncol=ncol(param)


 # structural error
 if( !(cname[1] %in% c("name","name")) ){
  cat("\n\n** error1(read.param) Must have a variable named 'name'"
     , " as the first variable.\n\n")
  return()
 }
 if( !(cname[2] == "type") ){
  cat("\n\n** error1(read.param) Must have a variable named 'type'"
    , " as the second variable.\n\n")
  return()
 }
 if( !(cname[3] == "ncat") ){
  cat("\n\n** error1(read.param) Must have a variable named 'ncat'"
    , " as the third variable.\n\n")
  return()
 }

 if( length(grep("^p[[:digit:]]$",cname)) == 0 ){
  cat("\n\n** error1(read.param) Must have a variable named 'p_num_'"
    , " at colums 4, 5, ....\n\n")
  return()
 }
 ncat=param$ncat
 locNO=which( param$type  %in% c("N","NO","O") )
 if( !any(cname %in% paste("p",max(ncat),sep="")) ){
  cat("\n\n** error1(read.param) Too few item paramters given."
      , " Must have 'p1' through ", paste("p",max(ncat),sep=""), ".\n\n")
  return()
 }
 if( length(locNO) > 0 )
  if( !any(cname %in% paste("p",max(ncat[locNO]-1)*2,sep="")) ){
   cat("\n\n** error1(read.param) Too few item paramters given."
       , " Must have 'p1' through "
       , paste("p",max(ncat[locNO]-1)*2,sep=""), ".\n\n")
   return()
  }



 # item names and types
 type=param$type
 locB=grep("^B[[:digit:]]*$", type)

 # use 0 for the c-parameter
 if( length(locB) > 0 ) param[locB,6][is.na(param[locB,6])]=0

 if( print ){
  cat("\n\n The following item parameter data frame was created from ",
      infile,".\n\n")
  Print(param)
 }

 return( param )

} # end of read.param




comments('
infile="d:/Rpgm/calr/param_BGPN0.dat"
param=read.param( infile, skip=1, print=1 )
str(param)


infile="d:/Rpgm/calr/param_BGP5.dat"
param=read.param( infile, print=1 )
str(param)
')


comments('
infile="d:/Rpgm/calr/param_BGPN1.dat"
param=read.param( infile, print=1 )
str(param)

#irf(param,plot=1)
')

infile="RpgmOLD/param_BGP5.dat"
infile="RpgmOLD/param_B1.dat"
param=read.param( infile, print=1, skip=0 )
str(param)




