#' Checking and Subsetting the Parameter Data Frame
#'
#' @param param Item Parameter Data Frame
#' @param type Item type
#' @param modulename name of the module which calls this function
#' @param tomat = 0 not to convert the result to a matrix.
#'
#' @details
#' This function checks if the input param is a valid parameter data frame
#' with name, type, ncat, p1, p2, ...,  variables.\cr
#' When type = (B | B3 | Bn | Bn3 | G | Gn | P |  PN | N) is given,
#' the items with the given type will be selected
#' and the parameter matrix (w/o name and type) will be returned.
#' \cr In this case, the item name will be stored in the rownames.
#' \cr
#' Otherwize, the input parame will be returned.
#' \cr
#' Note that type="B" means "B" or "B2" or "B3" and
#' type="B3" means "B3" only and that
#'  type="Bn" means "Bn" or "Bn2" or "Bn3" and
#' type="Bn3" means "Bn3" only.
#'
#' @examples
#' checkparam(paramS1, type="B" )
#' checkparam(paramS1[1,], type="B" )
#' checkparam(paramS1[3,], type="B" )
#'
#' @return
#' A parameter matrix or data frame.
#'
#' @export
#'

checkparam <- function( param, type=NULL, modulename="checkparam"
                      , tomat=1 ){
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
 #
 # Args:
 #
 #   param       parameter data frame name
 #   type        item type to be selected or NULL to select all
 #               B, B2, B3, G, Gn, P, PN or N  or ALL
 #   modulename  name of the module from which checkparam is called
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

 # return null if param is empty
 if( is.null(param) ){
  cat(err, " Empty data frame or matrix:", substitute(param),"\n")
  return(NULL)
 }

 # param given as data frame or not
 if( is.data.frame(param) ){
  #
  #  Here, param is a data frame
  #

  # error message root
  err=paste("\nerror(",modulename,")",sep="")

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
       " item param/weight data frame must have a variable named 'p1'.\n")
   return(NULL)
  }

  if( !is.null(type) ){

   if(  type == "ALL" ){
    # return the data frame as it is:
    # assuming checkparam is called from irf.
    return(param)
   }
   if( type == "B" ) loc=grep("^B[[:digit:]]*$", param$type)
   else if( type == "Bn" ) loc=grep("^Bn[[:digit:]]*$", param$type)
   else if( type == "G" ) loc=grep("^G", param$type)
   else              loc=which( param$type == type )
   if( length(loc) <= 0 ){
    cat(err, " There are no type '", type,"' items in the data frame.\n"
        , sep="")
    return(NULL)
   }

   param=param[loc,,drop=0]

   # set unspecified c-parameter to zero
   if( substr(type,1,1) == "B" ){
    if( is.null(param$p3) ) param$p3=0
    else param$p3[is.na(param$p3)]=0
   }

  } # type is given

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



} # end of checkparam









param=paramS1
str(param)
res=checkparam(param,  , "checkparam")
Print(res,class(res))
str(res)
checkparam(res,  , "checkparam")




comments(
'

plotchar=c(18,19)
plotchar=1:2
plotchar=c(16,17,18,15,21,24,23,22)
colors=c("black","blue")
linetype=1:2
matplot( res$theta,res$ICRF[,1:2],type="b", pch=plotchar, col=colors
        , lty=linetype )

legend( range(res$theta)[1], 0.8, 0:1, cex=0.9, col=colors, lty=linetype
                , pch=plotchar,  title="cat" )






loceq( letters, letters[c(2,3,5,2,3)] )

loceq( letters[c(2,3,5,2,3)], letters )


seed=1701
set.seed(seed)

long= sample(0:9,10,replace=1)
short=c(2,6)
Print(long,short,loceq(short,long), loceq(long,short))

short=c(6,2)
Print(long,short,loceq(short,long), loceq(long,short))


')










