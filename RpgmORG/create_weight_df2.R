#' Create Weight Data Frame
#'
#' @param param Parameter data frame with which the weight data frame be
#' associated.
#' @param w Item weight vector, or a scalar.
#' @param vvec Item category weight vector, or NULL
#'
#' @details
#' Using the information in param, the natural item category weights will be
#' generated.
#'
#' The non NA values of the item category weights must be stored in \code{vvec}.
#' The length of non NA item category weight of item \code{i} is given by
#' \code{param$ncat[i]}.
#'
#' \code{w} can be a scalar.
#'
#' @return A list of \cr
#' weight The weight data frame created. \cr
#' minmax A vector containing the minimum and the maximum score of the test.\cr
#' minmax_i A matrix containing the minimum and the maximum score of each item
#' before applying the item weight \code{w}..
#'
#' @examples
#' create_weight_df( paramS1 )
#' create_weight_df( paramS1, w=c(2,1,2,1) )
#'
#' create_weight_df( paramS1, vvec=c(0,1, -1,1, 0,1,2,3, -2,2,4,6) )
#'
#' @export
#'

create_weight_df <- function( param, w=1, vvec=NULL ){
 # create weight data frame
 # Shin-ichi Mayekawa
 # 20230715cot,16cot
 #

 # constants
 nitem=nrow(param)
 ncat=param$ncat
 iname=param$name
 itype=param$type

 # weight is a vector, scalar or NULL
 if( is.null(w) ) w=1
 if( is.vector(w) ) w=matrix(w,nitem)
 rownames(w)=iname; colnames(w)="w"
 v=matrix(0:(max(ncat)-1),nitem,max(ncat),byrow=1)
 rownames(v)=iname; colnames(v)=paste("v",0:(max(ncat)-1),sep="")
 for( i in 1:nitem ){
  if( ncat[i]+1 <= ncol(v) ) v[i,(ncat[i]+1):ncol(v)]=NA
 }

 # item category weight
 if( !is.null(vvec) ){
  if( is.vector(vvec) && length(vvec)==length(v[!is.na(v)]) ){
   toP <- cumsum(ncat)
   fromP <- toP - ncat + 1
   for (j in 1:nitem) {
    v[j, 1:ncat[j]] <- vvec[fromP[j]:toP[j]]
   }
  }
 }

 # output data frame
 weight=data.frame(name=iname,type=itype, ncat, w, v)

 # max score
 temp=find_minmax_score( weight )
 minmax=temp$minmax
 minmax_i=temp$minmax_i

 res=named_list( weight, minmax, minmax_i )
 return( res )


} # create_weight_df


create_weight_df( paramS1 )
create_weight_df( paramS1, w=c(2,1,2,1) )

create_weight_df( paramS1 )
create_weight_df( paramS1, w=c(1,2,1,2), vvec=c(0,1, -1,1, 0,1,2,3, -2,2,4,6) )

