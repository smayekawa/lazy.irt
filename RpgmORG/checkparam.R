checkparam <- function( param, type=NULL, modulename="checkparam" ){
 # check parameter data frame
 # Shin-ichi Mayekawa
 # 120224,27,29
 # name of giid: 120306
 # remove form/gfid: 120308
 # new types: 120911
 # type Bxx items: 120921
 # bugfix: 20150611
 #
 # Args:
 #
 #   param       parameter data frame name
 #   type        item type to be selected or NULL to select all
 #               B, G, P, PN or N  or ALL
 #   modulename  name of the module from which checkparam is called
 #
 #
 # Value:
 #
 #  parameter matrix of the form
 #     ncat p1 p2 ....
 #   or
 #  parameter data frame if type == "ALL"
 #
 #  If param is a matrix, param wll be returned as it is.
 #  If anything is wrong, NULL will be returned.
 #

 err=paste("\nerror(",modulename,")",sep="")

 if( is.null(param) ){
  cat(err, " Empty data frame or matrix:", substitute(param),"\n")
  return(NULL)
 }

 # param given as data frame or not
 if( is.data.frame(param) ){
  #
  #  Here, param is a data frame
  #
  # check if it is a genuine data frame
  if( is.null(param$name) & is.null(param$giid) ){
   cat(err,
       " item param/weight data frame must have a variable named 'name'"
     , " or 'giid'.\n")
   return(NULL)
  }
  else if( is.null(param$name) ){
   param=data.frame(name=param$giid,param)
   param=param[,-which(colnames(param) == "giid")]
  }

  if( !is.null(param$gfid) ) param=param[,-which(colnames(param) == "gfid")]
  if( !is.null(param$form) ) param=param[,-which(colnames(param) == "form")]

  if( is.null(param$type) ){
   cat(err,
       " item param/weight data frame must have a variable named 'type'.\n")
   return(NULL)
  }
  if( is.null(param$ncat) ){
   cat(err,
       " item param/weight data frame must have a variable named 'ncat'.\n")
  return(NULL)
  }
  if( is.null(param$p1) ){
   cat(err,
       " item param/weight data frame must have a variable named 'p1'.\n")
   return(NULL)
  }

  if( !is.null(type) ){
   if(  type == "ALL" ){
    # return the data frame as it is:   assuming checkparam is called from irf.
    return(param)
   }
   if( type == "B" ) loc=grep("^B[[:digit:]]*$", param$type)
   else              loc=which( param$type == type )
   if( length(loc) <= 0 ){
    cat(err,
        " There are no type '",type,"' items in the data frame.\n")
    return(NULL)
   }
   param=param[loc,,drop=0]
  }

  # convert data frame to matrix
  locparamnum=grep("^p[[:digit:]]",colnames(param))
  lenp=ncol(param)
  parammat=as.matrix( cbind(param$ncat ,param[,locparamnum]) )
  colnames(parammat)=c("ncat",colnames(param)[locparamnum])

  if( !is.null(param$name) ) rownames(parammat)=param$name
  else rownames(parammat)=param$giid

  return(parammat)

 } # end of data frame
 else{
  # param is a matrix
  return(param)
 }



} # end of checkparam




res=checkparam(param,  , "checkparam")
Print(res,class(res),str(res))

checkparam(res,  , "checkparam")
