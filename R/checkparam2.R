#' Checking and Subsetting the Parameter Data Frame (New Version)
#'
#' @param param Item Parameter Data Frame
#' @param type Item type
#' @param modulename name of the module which calls this function
#' @param tomat = 0 not to convert the result to a matrix.
#' @param printerror = 0 to suppress error messages.
#'
#' @details
#' This function checks if the input param is a valid parameter data frame
#' with name, type, ncat, p1, p2, p3, ...,  variables.\cr
#' When type = (B | B3 | Bn | Bn3 | G | Gn | P |  PN | N) is given,
#' the items with the given types will be selected
#' and the parameter data frame or matrix (w/o name and type) will be returned.
#' In the resulting matrix, the item names will be stored as the row names.
#' \cr
#' Note that \code{type="B"} means both "B" and "B3".
#' To select 2PLM items, use \code{type="B2"}. \cr
#' Also, \code{type="Bn"} means both "Bn" and "Bn3".
#' To select 2PNM items, use \code{type="Bn2"}. \cr
#' \cr
#' If \code{type="ALL"}, regardless of \code{tomat}, the input data frame
#' will be returned after validity check.
#'
#' Note that, when \code{type} contains more than 2 types, the result will be
#' sorted according to the order of the appearance of the types in \code{type}.
#'
#'
#' @examples
#' checkparam2(paramA1, type="B" )
#' checkparam2(paramA1, type=c("B","P"), tomat=0 )
#' checkparam2(paramA1, type=c("Bn2","G"), tomat=0 )
#'
#' @return
#' A parameter matrix or data frame.
#'
#' @export
#'

checkparam2 <- function( param, type=NULL, modulename="checkparam2"
                      , tomat=1, printerror=1 ){
 # check parameter data frame
 # Shin-ichi Mayekawa
 # 120224,27,29
 # name of giid: 120306
 # remove form/gfid: 120308
 # new types: 120911
 # type Bxx items: 120921
 # bugfix: 20150611
 # tomat etc: 20161123
 # Bn and Gn: 20180211,15
 # error -> error1: 20230326
 # vectorize type arg: 20230327
 # check if param exists: 20230329
 #
 # Args:
 #
 #   param       parameter data frame name
 #   type        item type to be selected or NULL to select all
 #               B, B2, B3, Bn, Bn3, G, Gn, P, PN or N  or ALL
 #   modulename  name of the module from which checkparam2 is called
 #   tomat       = 1 or 0
 #
 #
 # Value:
 #
 #  parameter matrix of the form
 #     ncat p1 p2 ....
 #   or
 #  parameter data frame if type == "ALL" or tomat == 0
 #
 #  If param is a matrix, param wll be returned as it is.
 #  If anything is wrong, NULL will be returned.
 #

 dfname=deparse(substitute(param))
 exist=exists( dfname, envir=parent.frame() )

 # return null if param is empty
 if( !exist || is.null(param) ){
    cat("\nerror1:(checkparam) Empty data frame or matrix:"
        , substitute(param),"\n")
   return(NULL)
 }

 # param given as a data frame or not
 if( is.data.frame(param) ){
  #
  #  Here, param is a data frame
  #

  # error message root
  err=paste("\nerror1(",modulename,")",sep="")

  # check if it is a genuine data frame
  # Must have name or giid(cal) column.
  if( is.null(param$name) & is.null(param$giid) ){
   cat(err,
       " item param/weight data frame must have a variable named 'name'"
     , " or 'giid'.\n")
   return(NULL)
  }
  else if( is.null(param$name) ){
   # rename giid(cal) as name and remove giid(cal)
   param=data.frame(name=param$giid,param)
   param=param[,-which(colnames(param) == "giid")]
  }

  # remove gfid(cal) and form(cal) columns
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
       " item param/weight data frame must have variable named ",
       "'p1', 'p2', and 'p3'.\n")
   return(NULL)
  }

  if( !is.null(type) ){

   if( any( toupper(type) == "ALL" ) ){
    # return the data frame as it is:
    # assuming checkparam is called from
    #  irf, dirf, iif, obscore, info_func, rel_irt.
    # set unspecified c-parameter to zero
    if( is.null(param$p3) ) param$p3=0
    loc=grep("(^B$)|(^Bn$)", param$type)
    param$p3[loc]=0
    return(param)
   }

   lent=length(type)
   loc=NULL
   for( i in 1:lent ){
    if( type[i] == "B2" ) loci=grep("^B$", param$type)
    else if( type[i] == "B" ) loci=grep("^B[[:digit:]]*$", param$type)
    else if( type[i] == "Bn2" ) loci=grep("^Bn$", param$type)
    else if( type[i] == "Bn" ) loci=grep("^Bn[[:digit:]]*$", param$type)
    else if( type[i] == "G" ) loci=grep("^G$", param$type)
    else loci=which( param$type == type[i] )
    loc=c(loc,loci)
   }

   if( length(loc) <= 0 ){
    if( printerror == 1 )
     cat(err, " There are no type '", type,"' items in the data frame.\n"
         , sep="")
    return(NULL)
   }

   param=param[loc,,drop=0]

  } # type is given

  # set unspecified c-parameter to zero
  if( is.null(param$p3) ) param$p3=0
  loc=grep("(^B$)|(^Bn$)", param$type)
  if( length(loc) > 0 ) param$p3[loc]=0

  # remove all NA columns
  # Print(apply( param, 2, function(x) !all(is.na(x)) ))
  param=param[,apply( param, 2, function(x) !all(is.na(x)) ), drop=0]

  if( tomat == 1 ){
   # convert data frame to matrix
   locparamnum=grep("^p[[:digit:]]",colnames(param))
   lenp=ncol(param)
   parammat=as.matrix( cbind(param$ncat ,param[,locparamnum]) )
   colnames(parammat)=c("ncat",colnames(param)[locparamnum])
   if( !is.null(param$name) ) rownames(parammat)=param$name
   else rownames(parammat)=param$giid
   return( parammat )
  }
  else return(param)

 } # end of data frame
 else{

  # param is a matrix: return as it is.
  return( param )

 } # end of mat



} # end of checkparam2


