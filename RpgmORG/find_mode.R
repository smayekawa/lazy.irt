#' Find the mode of icrf
#'
#' @param param Item Paramter Data Frame
#' @param print = 1 to print the result
#'
#' @details
#' When item type is "P", the original step parameter b_k
#' is the intersection of P_{k-1} and P_k, k=1,2,...,ncat-1 \{cr
#' When item type is "G" , the mode of P_k is given as (b_k+b_{k+1})/2 .
#' \cr
#' Note that, b_k is stored as p_{k+1} in the param data frame
#' since p_1 is always the slope.
#' \cr
#' The modes of type "P" items and the intersections of type "G" or "B3"
#' items are found numerically.
#'
#' \cr\cr
#' The numerical root finding may result in an error.
#'
#' @return
#' A matrix of nitem x max(ncat) consisting of the modes
#' corresponding to the 0-th, 1st, 2nd, ..., and ncat-th category.
#' \cr Note that the first (0) and the last (ncat) columns are filled with NAs.
#'
#' @examples
#' mode <- find_mode( paramS1, print=1 )
#'
#' @export
#'

find_mode <- function( param, print=0 ){
 # find the mode of icrf
 # Shin-ichi Mayekawa
 # 20161122dnc
 #


 # mode function
 mm <- function( theta, Pj, k ){
  # The category number k:  0 <= k <= ncat-1
  return( irf( Pj, theta=theta, print=0 )$ICRF[k+1] )
 } # end of mm


 nitem=nrow(param)
 ncat=param$ncat
 type=param$type
 locp1=which(colnames(param)=="p1")

 # abp=param[,colnames(param) %in% paste("p",1:nitems,sep=""), drop=0]
 # abp=as.matrix(abp)

 # Print(ncat,nitem)

 mode=matrix(NA,nitem,max(ncat))
 rownames(mode)=param$name
 colnames(mode)=0:(max(ncat)-1)


 for( j in 1:nitem ){

  ncatj=ncat[j]; typej=type[j]
  paramj=param[j,]

  if( typej == "G" ){
   for( k in 1:(ncatj-2) ){
    mode[j,k+1]=(paramj[1,locp1+k]+paramj[1,locp1+k+1])/2
   }
  } # end of G
  else if( typej %in% c("P","PN", "PN0") ){
   for( k in 1:(ncatj-2) ){
    res=optimize( mm, c(-15,15), Pj=paramj, k=k, maximum=1 )
    mode[j,k+1]=res$maximum
   }
  } # end of P

 } # end if j for item

 if( print ){
  Print( mode )
 }

 return( mode)

} # end of find_mode




#' Find the intersection of icrf
#'
#' @param param Item Paramter Data Frame
#' @param print = 1 to print the result
#'
#' @details
#' When item type is "P", the original step parameter b_k
#' is the intersection of P_{k-1} and P_k, k=1,2,...,ncat-1 \{cr
#' When item type is "G" , the mode of P_k is given as (b_k+b_{k+1})/2 .
#' \cr
#' Note that, b_k is stored as p_{k+1} in the param data frame
#' since p_1 is always the slope.
#' \cr
#' The modes of type "P" items and the intersections of type "G" or "B3"
#' items are found numerically.
#'
#'
#' \cr\cr
#' The numerical root finding may result in an error.
#'
#' @return
#' A matrix of nitem x max(ncat)-1 consisting of the intersections of
#' 0-th and 1st, 1st and 2nd, ..., (ncat-1)-th and ncat-th categories.
#'
#' @examples
#' intersection <- find_intersection( paramS1, print=1 )
#'
#' @export
#'

find_intersection <- function( param, print=0 ){
 # find the mode of icrf
 # Shin-ichi Mayekawa
 # 20161122dnc
 #



 # intersection function
 ii <- function( theta, Pj, k1, k2=k1+1 ){
  # The category number k1:  0 <= k1 <= ncat-1
  dif=irf( Pj, theta=theta, print=0 )$ICRF[k1+1]-
   irf( Pj, theta=theta, print=0 )$ICRF[k2+1]
  return( dif )
 } # end of ii


 nitem=nrow(param)
 ncat=param$ncat
 type=param$type
 locp1=which(colnames(param)=="p1")

 # Print(ncat,nitem)

 intersect=matrix(NA,nitem,max(ncat)-1)
 rownames(intersect)=param$name
 colnames(intersect)=paste(0:(max(ncat)-2),1:(max(ncat)-1),sep="")

 for( j in 1:nitem ){

  ncatj=ncat[j]; typej=type[j]
  paramj=param[j,]


  if( typej == "PN" ){
   paramj=convPN2P( paramj )
  }

  if( typej %in% c("G", "B3") ){
   # G or B3
   modej=find_mode( paramj )
   modej[1]=-5; modej[ncatj]=5
   for( k in 0:(ncatj-2) ){
    res=uniroot( ii, c(modej[k+1],modej[k+2]), Pj=paramj, k1=k )
    intersect[j,k+1]=res$root
   }
  } # end of G
  else if( typej %in% c("P","PN", "PN0","B") ){
   # P
   for( k in 0:(ncatj-2) ){
    intersect[j,k+1]=paramj[1,locp1+k+1]
   }
  } # end of P

 } # end if j for item

 if( print ){
  Print( intersect )
 }

 return( intersect)

} # end of find_intersection










intersect=find_intersection( param )
Print(intersect, fmt="8.3")





comments(
'



intersect=find_intersection( paramS1 )
Print(intersect, fmt="8.3")


mode=find_mode( paramS1 )
Print(mode)






temp=c("
 name type ncat p1 p2 p3  p4 p5 p6
       Q1 G  6  1  -3 -1   0  1  4
       Q2 P  6  1  -3 -1   0  1  4
       Q0 B3 2  1   0  0.2,  ,  ,  ,
       ")




temp=c("
name type ncat p1 p2 p3 p4
Q1 G 4  1,  -1, 1, 2
Q2 P 4  1,  -1, 1, 2
")


param=cards(temp,header=1)
ncat=param$ncat[1]

mode=find_mode( param[1,])

temp=irf( param[1,], print=0, plot=1, npoints=51)
ypos=0; tck=0.03
for( k in 1:(ncat-1)+1 ){
 lines( rep(mode[k],2), c(ypos-tck, ypos+tck), type="l", col="red" )
}


'
)





