#' Dummmy Expantion of Categorical Data
#'
#' @details
#'  In compressed format, Uc[i,j]=k means that subject i answers to item j
#'  by selecting option k, k=0,1,...,ncat[j]. \cr
#'  This function expands this compressed format to uncompressed
#'  dummy expanded format: \cr
#'  That is, each item takes ncat[j] colums of U[, fromP[j]:toP[j]], \cr
#'  and the k+fromP[j]-th colum will be 1, and 0 elsewhere.
#' @param Uc Item Response matrix
#' @param ncat a vector consisting of # of categories of each item
#' @param type NOT used.
#' @param zero = 0 to exclude the zero-th category from output
#'
#' @return  A list of \cr
#' U, ncat, type, zero, fromP, toP
#'
#' @examples
#' # generate item response
#' # in uncompressed format
#' set.seed(1701)
#' res <- gendataIRT( 1, paramS1, npoints=7, thmin=-3, thmax=3 )
#' U=res$U
#' # in compressed format
#' set.seed(1701)
#' res1 <- gendataIRT( 1, paramS1, npoints=7, thmin=-3, thmax=3, compress=1 )
#' Uc=res1$U
#'
#' # from compressed to uncompressed
#' U1=dummy_expand(Uc)$U
#' # back to compressed
#' Uc1=dummy_compress(U)$Uc
#' Print(U-U1)
#' Print(Uc-Uc1)
#'
#' # zero=0 option
#' U2=dummy_expand(Uc, zero=0)$U
#' Uc2=dummy_compress(U2, zero=0, ncat=paramS1$ncat)$Uc
#' Print(Uc-Uc2)
#'
#'
#' @export
#'

dummy_expand <- function( Uc, ncat=NULL, type=NULL, zero=1 ){
 # dummy expand the categorical valued variable u
 # Shin-ichi Mayekawa
 # 20121022
 # matrix input: 121024
 # missing value: 121025
 # make the data integer when possible: 121126
 # vectorize fromP: 20151113
 # 20161210

 nitems=ncol(Uc)
 npoints=nrow(Uc)
 iname=colnames(Uc)

 if( is.null(ncat) ){
  ncat=matrix(apply(Uc,2,max,na.rm=1),,1)+1
 }
 toP=cumsum(ncat)
 fromP=c(toP-ncat+1)
 # Print(nitems,npoints,ncat,toP,fromP, iname)
 U=matrix(as.integer(0),npoints,toP[nitems])
 for( j in 1:nitems ){
  uj=Uc[,j]
  locobs=!is.na(uj)
  min=0; max=ncat[j]-1
  factuj=factor(uj, levels=min:max)
  Uj=model.matrix(~factuj-1)
  #Print(uj,factuj,Uj)
  #  attributes(Uj)$assign=NULL
  #  attributes(Uj)$contrasts=NULL
  #  colnames(Uj)=min:max
  U[locobs,fromP[j]:toP[j]]=as.integer(Uj)
 }
 if( is.null(iname) )
  iname=paste("Q",1:nitems,sep="")
 Pname=character(sum(ncat))
 catname=outer(paste(iname,"_",sep=""),0:max(ncat),paste,sep="")
 for( j in 1:nitems ){
  Pname[fromP[j]:toP[j]]=catname[j,1:ncat[j]]
 }
 colnames(U)=Pname
 rownames(U)=rownames(Uc)

 if( zero != 1 ){
  # remove the 0-th category
  U=U[,-fromP]
  fromP=fromP-(0:(nitems-1))
  toP=toP-(1:nitems)
 }

 return( list(U=U, ncat=ncat, type=type, zero=zero, fromP=fromP, toP=toP) )

} # end of dummy_expand




#' Compress the Dummmy Expanded Categorical Data
#'
#' @details
#'  In compressed format, Uc[i,j]=k means that subject i answers to item j
#'  by selecting option k, k=0,1,...,ncat[j]. \cr
#'  This function recovers this compressed format from the dummy expanded
#'  uncompressed format: \cr
#'  That is, each item takes ncat[j] colums of U[, fromP[j]:toP[j]], \cr
#'  and the k+fromP[j]-th colum will be 1, and 0 elsewhere.
#' @param U The Item Response matrix in Compressed Format
#' @param ncat a vector consisting of # of categories of each item
#' @param type NOT used.
#' @param zero = 0 if the input does not have the zero-th column.
#'
#' @return  A list of \cr
#' Uc, ncat, type, zero, fromP, toP
#'
#' @examples
#' # generate item response
#' # in uncompressed format
#' set.seed(1701)
#' res <- gendataIRT( 1, paramS1, npoints=7, thmin=-3, thmax=3 )
#' U=res$U
#' # in compressed format
#' set.seed(1701)
#' res1 <- gendataIRT( 1, paramS1, npoints=7, thmin=-3, thmax=3, compress=1 )
#' Uc=res1$U
#'
#' # from compressed to uncompressed
#' U1=dummy_expand(Uc)$U
#' # back to compressed
#' Uc1=dummy_compress(U)$Uc
#' Print(U-U1)
#' Print(Uc-Uc1)
#'
#' # zero=0 option
#' U2=dummy_expand(Uc, zero=0)$U
#' Uc2=dummy_compress(U2, zero=0, ncat=paramS1$ncat)$Uc
#' Print(Uc-Uc2)
#'
#'
#' @export
#'

dummy_compress <- function( U, ncat=NULL, type=NULL, zero=1 ){
 # compress the dummy_expanded data
 # Shin-ichi Mayekawa
 # 20161210
 #

 # recover ncat
 if( is.null(ncat) ){
  if( zero == 0 ){
   cat("\nerror1(dummy_compress): zero must be 1 when ncat is not given.\n")
   stop()
  }
  toP=which( apply(apply(U,1,cumsum),1,var) == 0 )
  if( length(toP) == 0 ){
   cat("\nerror1(dummy_compress): Input matrix is not a dummy expanded data."
       , "\n")
   cat(" zero-th category may be missing.\n")
   stop()
  }
  fromP=c(1,toP+1)[-length(toP)-1]
  ncat=toP-fromP+1
 }
 else{
  toP=cumsum(ncat)
  fromP=toP-ncat+1
  if( zero != 1 ){
   fromP=fromP-(0:(nitems-1))
   toP=toP-(1:nitems)
  }
 }

 # Print(zero, ncat, fromP, toP)
 # const
 npoints=nrow(U)
 nitems=length(ncat)
 iname=colnames(U)[fromP]
 iname=sub("_[[:digit:]]","",iname)

 # create Uc from U
 Uc=matrix(0,npoints,nitems)

 if( zero == 1 ){
  for( j in 1:nitems ){
   loc=apply(U[,fromP[j]:toP[j],drop=0], 1, which.max)
   Uc[,j]=loc-1
  }
 }
 else if( zero != 1 ){
  for( j in 1:nitems ){
   loc=apply(U[,fromP[j]:toP[j],drop=0], 1, which.max)
   Uc[,j]=loc
   loc=which( apply(U[,fromP[j]:toP[j],drop=0], 1, sum) == 0 )
   Uc[loc,j]=0
  }
 }

 colnames(Uc)=iname

 return( named_list(Uc, ncat, type, zero, fromP, toP) )

} # end of dummy_compress





# generate item response
# in uncompressed format
set.seed(1701)
res <- gendataIRT( 1, paramS1, npoints=7, thmin=-3, thmax=3 )
U=res$U
# in compressed format
set.seed(1701)
res1 <- gendataIRT( 1, paramS1, npoints=7, thmin=-3, thmax=3, compress=1 )
Uc=res1$U

# from compressed to uncompressed
U1=dummy_expand(Uc)$U
# back to compressed
Uc1=dummy_compress(U)$Uc
Print(U-U1)
Print(Uc-Uc1)

# zero=0 option
U2=dummy_expand(Uc, zero=0)$U
Uc2=dummy_compress(U2, zero=0, ncat=paramS1$ncat)$Uc
Print(Uc-Uc2)


