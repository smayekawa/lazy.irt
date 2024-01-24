#' Create Weight Data Frame
#'
#' @param param Parameter data frame with which the weight data frame be
#' associated.
#' @param w Item weight vector
#'
#' @details
#' Using the information in param, the natural item category weights willbe
#' generated.
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
#' @export
#'

create_weight_df <- function( param, w=1 ){
 # create weight data frame
 # Shin-ichi Mayekawa
 # 20230715cot
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
  weight=cbind(ncat,w,v); rownames(weight)=iname
  colnames(weight)=c("ncat","w",colnames(v))
  weight=data.frame(name=iname,type=itype, w, v)


 # max score
 temp=find_minmax_score( weight )
 minmax=temp$minmax
 minmax_i=temp$minmax_i

 res=named_list( weight, minmax, minmax_i )
 return( res )


} # create_weight_df


