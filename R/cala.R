#' Calibration of Independently Estimated Sets of Item Parameters
#' by Minimizing the LS Criterion Defined in terms of Parameter Values
#'
#'
#' @param indata   data frame containing at least following: \cr
#'             giid  gfid  type ncat  a  p1 p2 .... \cr
#'  indata can be created from a parameter data frame
#'  by renaming name as giid and adding gfid.
#'  It is assumed that indata is sorted by gfid.
#' @param baseform    form number or 0
#' @param nsubj a vector consisting of the # of subjects in each form or 1.
#' \cr To be used when calculating the mean of theta of the combined group.
#' @param maxiter     max of iteration for PCA
#' @param eps         eos for PCA
#' @param init        not used
#' @param print       > 0 to print results
#' @param debug = 1 to print intermediate result
#'
#' @details
#'
#'   This program finds u,v and b such taht \cr
#'      bhata[j,g] = u[g] + b[j] v[g] + error \cr
#'      j=1,2,..., # of items and  g=1,2,..., # of forms, \cr
#'   using PCA with missing data, where bhat is the collection of estimated
#'   b-parameters. \cr
#'  \cr
#'    a[g] will be estimated as the average of ahat[j,g]/r[g] , \cr
#'    c[g] will be estimated as the average of chat[j,g] , \cr
#'   where r[g]=1/v[g]
#'
#' @return A list of: \cr
#'  param Data frame containing the equated item parameters \cr
#'  uv  Transformations from common scale to each scale\cr
#'  qr  Transformations from each scale to common scale\cr
#'  nfini Vector of # of forms including each item \cr
#'  niinf Vector of # of items included in each form \cr
#'  fini  List of the Forms including each item. \cr
#'  iinf  List of the Items included in each form.\cr
#'  iftable  Matrix of item x form indicating the inclusion pattern.\cr
#'  toP  Vector of stating locations in ParamHat \cr
#'  fromP  Vector of ending locations in ParamHat  \cr
#'  iterPCA  # of iterations needed for convergence. \cr
#'  nbPCA  3 of b-parameters processed by PCA.\cr
#'  nisolo # of items that are included in only one form. \cr
#'  ParamHat Item Paramter matrix analyzed by PCA \cr
#'  rmse LS criterion minimized.
#'
#' @references
#' Sayaka Arai & Shin-ichi Mayekawa (2011)
#' A comparison of equating methods and linking designs for developing
#' an item pool under item response theory.
#' Behaviormetrika, Vol. 38 No. 1, pp. 1-16 \cr
#'
#' Shin-ich Mayekawa (1991) Chapter 4. Parameter Estimation.
#' In Shiba Ed. Item Response Theory. pp87-129. (in Japanese)
#'
#' @examples
#' # Binary Items
#' res=cala( indata=paramCal1, baseform=1, print=2  )
#'
#' # Mixture of Item Types
#' res=cala( indata=paramCal2, baseform=1, print=2  )
#'
#' @export
#'

cala <- function( indata, baseform=0, nsubj=1
                  , maxiter=100, eps=1e-5, init=1
                  , print=1, debug=0  ){
 # cala: Calibration of Separately Estimated Item Parameters (pArameter)
 # Shin-ichi Mayekawa
 # fortran version: 870903 -- 080727
 # 120227,29,0301
 # solo items: 120303,04
 # baseform=0: 120304
 # when PCA returns negative v: 120908
 # type Bxx items: 120921
 # output many results: 120922
 # cala.core renamed as cala: 120926,1004,05
 # bugfix: 121030
 # output rmse: 121103
 # Nominal Response Items: 121130  ***************** not ready yet
 # nsubj and nsb: 20161121
 #
 #
 #
 # Args:
 #
 #   indata      data frame containing at least following
 #               giid  gfid  type ncat  a  p1 p2 ....
 #
 #               indata can be created from a parameter data frame
 #               by renaming name as giid and adding gfid
 #
 #               It is assumed that indata is sorted by gfid.
 #
 #
 #   baseform    form number of 0
 #
 #   maxiter     max # of iteration for PCA
 #   eps         eos for PCA
 #   init        not used
 #
 #   print       > 0 to print results
 #   debug
 #
 #
 #  This program finds u,v and b such taht
 #     bhata[j,g] = u[g] + b[j] v[g] + error
 #  using PCA with missing data, where bhat is the collection of
 #  "b" parameters.
 #
 #   a[g] will be estimated as the average of ahat[j,g]/r[g] ,
 #   c[g] will be estimated as the average of chat[j,g] ,
 #  where r[g]=1/v[g]
 #
 #
 #
 # subscript j is for item, g for form, q is for theta point
 #
 #

 if( debug ) Print("****** top of cala", indata)

 # const
 ne=nrow(indata)

 # add seqnum
 if( length(indata$seqnum) == 0 )
  indata=cbind(seqnum=1:ne,indata)

 # remove liid lfid
 indata$liid=NULL
 indata$lfid=NULL

 # IDI -> giid, IDT -> gfid

 # ids in original order
 if( length(indata$giid) == 0 ) indata$giid=indata$name
 if( length(indata$gfid) == 0 ) indata$gfid=indata$form
 giid=as.character(indata$giid)
 gfid=as.character(indata$gfid)

 # unique giid and gfid in the order of indata
 ugiid=unique(giid)
 ugfid=unique(gfid)
 nitems=length(ugiid)
 nforms=length(ugfid)

 if( length(nsubj) == 1 ){
  if( nsubj == 0 ) nsubj=1
 }
 nsubj=matrix(nsubj,1,nforms)
 rownames(nsubj)=""; colnames(nsubj)=ugfid

 #Print(nitems,nforms,ne, digits=3)


 # add local id (liid and lfid) to indata by merging

 # liid to giid table:   liid=1,2,...,max_unique_giid
 # unique item names in the original order
 lgitable=data.frame( liid=as.numeric(1:nitems), giid=ugiid
                      , stringsAsFactors=0  )
 # merge liid
 indata=merge(lgitable, indata, by="giid",sort=0)

 # lfid to gfid table:   lfid=1,2,...,nforms
 lgftable=data.frame( lfid=as.numeric(1:nforms), gfid=ugfid
                      , stringsAsFactors=0  )
 # merge lfid
 indata=merge(lgftable,indata,by="gfid",sort=0)


 # sort back to the original order
 indata=indata[order(indata$seqnum),]
 rownames(indata)=c(1:ne)


 # # of items per form, # of forms per item
 nfini=matrix(table(indata$liid)) # of forms including item j
 niinf=matrix(table(indata$lfid)) # of items included  in form g

 # location of items in indata and  form names including item j
 loci=NULL; fini=NULL
 for( j in 1:nitems ){
  loc=which(indata$liid == j)
  loci=append(loci,list(loc))
  fini=append(fini,list(indata$lfid[loc]))
 }

 # location of forms in indata and item names included in form g
 locf=NULL; iinf=NULL
 for( g in 1:nforms ){
  loc=which(indata$lfid == g)
  locf=append(locf,list(loc))
  iinf=append(iinf,list(indata$liid[loc]))
 }
 # Print(loci,locf)
 # Print("*",fini,iinf)

 # iftable: item-form inclusion table
 iftable=matrix(NA,nitems,nforms)
 for( j in 1:nitems ){
  iftable[ j,indata$lfid[ loci[[j]] ] ]=1
 }
 rownames(iftable)=ugiid
 colnames(iftable)=ugfid
 # Print(iftable)

 # ncat of each item (unique):  pick up the first element
 ncat=indata$ncat[ sapply( loci,"[[",1 ) ]

 # unique item type and ncat:
 # pick up the first list-elements of loci as subscript
 lgitable=data.frame(lgitable
                     , type=as.character( indata$type[ sapply( loci,"[[",1 ) ] )
                     , ncat=ncat
                     , stringsAsFactors=0 )


 if( debug ){
  Print(nitems,nforms,ne,"/",lgitable,lgftable
        , nfini,fini, niinf,iinf, loci, locf)
  Print(indata)
 }

 #
 # note
 #
 # lgitable[j,2] = the global item id of item j, j=1,2,...,nitems
 # lgftable[g,2] = the global form id of form g, g=1,2,...,nforms
 #
 # nfini[j]   = number of forms including item j, j=1,2,...,nitems
 # fini[j]    = local form ids including item j
 #                 formini[j] is a vector of length nformini[j]
 # loci[j]       = locations of item j in indata
 #                 indata$a[loci[j]] returns nformini[j] a parameters
 #
 # niinf[g]  = number of items included in form g, g=1,2,...,nforms
 # iinf[g]   = local item ids included in form g
 #                 iinf[g] is a vector of length niinf[g]
 # locf[g]       = locations of form g in indata
 #                 indata$a[locf[g]] returns nformini[j] a parameters
 #
 # iftable[j,g]  = 1 if item j is included in form g
 #



 # collect a-parameters in a matrix form like iftable
 atable=matrix(NA,nitems,nforms)
 for( j in 1:nitems ){
  atable[ j,indata$lfid[ loci[[j]] ] ]=indata$p1[loci[[j]]]
 }
 rownames(atable)=ugiid

 # collect c-parameters in a matrix form like iftable
 ctable=matrix(NA,nitems,nforms)
 for( j in 1:nitems ){
  if( length(grep( "^B[[:digit:]]*$", indata$type[ loci[[j]][[1]] ] )) > 0 )
   ctable[ j,indata$lfid[ loci[[j]] ] ]=indata$p3[ loci[[j]] ]
 }
 rownames(ctable)=ugiid


 # collect b-parameters for PCA
 # from and to without zero-th category
 toP=cumsum(ncat)
 fromP=toP-ncat+1
 fromP=fromP-0:(nitems-1)
 toP=toP-1:nitems


 # storage for b-parameter
 nbparam=sum(ncat-1)
 ParamHat=matrix(NA,nbparam,nforms)
 catname=outer(paste(ugiid,"_",sep=""),1:max(ncat),paste,sep="")

 # rnPH=unlist(
 #  mapply( function(x,y){ catname[x,y] }, 1:nitems, sapply(ncat-1,seq) )
 # )
 # Above does not work when all ncat elements are the same.
 # simplify=0 in sapply or mapply does not work.
 # make the result of sapply a list
 sss=sapply(ncat-1,seq,simplify=0)
 if( is.vector(sss) & !is.list(sss) ) sss=sapply(sss,list)
 else if( is.matrix(sss) ){
  if( ncol(sss) == 1 ) sss=list(as.vector(sss))
  else sss=unlist(apply(sss,2,list),recursive=0)
 }
 rnPH=unlist(
  mapply( function(x,y){ catname[x,y] }, 1:nitems, sss )
 )
 #Print(catname,sapply(ncat-1,seq))
 rownames(ParamHat)=rnPH
 colnames(ParamHat)=ugfid
 locp2=which(colnames(indata) == "p2")
 for( k in 1:ne ){
  liidk=indata$liid[k]
  lfidk=indata$lfid[k]
  ncatk=indata$ncat[k]
  ParamHat[fromP[liidk]:toP[liidk],lfidk]=
   unlist(indata[k,locp2:(locp2+ncatk-2)])
 }
 ParamHat=as.matrix(ParamHat)
 nf=rowSums(!is.na(ParamHat))
 # rnPH=rownames(ParamHat)

 # subset of ParamHat:  unique items for each form
 locsolo=which(nf == 1)
 loccommon=which(nf > 1)
 ParamHat1=ParamHat[locsolo,,drop=0]
 nsolo=nrow(ParamHat1)
 nisolo=length(which(nfini == 1))
 # Print(nfini,nisolo,nsolo,which(nfini==1,arr.ind=1),which(nfini==1))
 rownames(ParamHat1)=rownames(ParamHat)[locsolo]
 # common items
 ParamHat=ParamHat[loccommon,]


 # PCA of b-parameters of common items
 res=princ( ParamHat, ndim=1, estm=1, maxiter=maxiter, eps=eps, print=print-3 )
 iterPCA=res$iter
 rmse=res$rmse
 b=res$F; rownames(b)=rownames(ParamHat)
 u=res$mu; v=res$A
 if( any(v < 0) ){
  v=-v; b=-b
 }
 rm(res)
 q=-u/v; r=1/v

 if( nsolo > 0 ){
  # transform the b-parameters of unique items by q and r and append to b
  b1=matrix(0,nsolo,nforms)
  rownames(b1)=rownames(ParamHat1)
  for( j in 1:nsolo ){
   b1[j,]=ParamHat1[j,]*r+q
  }
  b1=matrix( rowSums(b1,na.rm=1), nsolo )
  rownames(b1)=rownames(ParamHat1)
  b2=matrix(NA,nbparam,1)
  b2[locsolo]=b1; b2[loccommon]=b
  b=b2
  rownames(b)=rnPH
 }

 # adjust parameter to base form
 if( baseform != 0 ){
  qb=q[baseform]; rb=r[baseform]
  b=(b-qb)/rb
  q=(q-qb)/rb
  r=r/rb
  v=1/r
  u=-q/r
 }
 else{
  # baseform == 0:   calculate the combined mean and std of theta
  nsb=as.vector(nsubj)
  nsb=nsb/sum(nsb)
  tmean=sum(q*nsb)
  tstd=sqrt( sum((r^2 + q^2)*nsb) - tmean^2 )
  # Print(tmean,tstd)

  # standardize all by the combined mean/std
  b=(b-tmean)/tstd
  u=u+v*tmean
  v=v*tstd
  q=-u/v; r=1/v
 }


 # final result
 rownames(b)=rnPH
 uv=cbind(u,v); qr=cbind(q,r)
 colnames(uv)=c("u","v"); rownames(uv)=ugfid
 colnames(qr)=c("q","r"); rownames(qr)=ugfid

 # convert b vector to a matrix
 bmat=matrix(NA,nitems,max(ncat-1))
 rownames(bmat)=ugiid
 for( j in 1:nitems ){
  bmat[j,1:(ncat[j]-1)]=b[fromP[j]:toP[j]]
 }

 # convert a and c paramters: cannot use atable%*%diag(1/r) because of NAs
 for( g in 1:nforms ){
  atable[,g]=atable[,g]/r[g]
 }
 a=as.matrix( rowMeans(atable,na.rm=1), nitems )
 c=as.matrix( rowMeans(ctable,na.rm=1), nitems )
 c[is.nan(c)]=NA
 rownames(a)=ugiid
 rownames(c)=ugiid

 param=as.data.frame(cbind(a,bmat))
 # 3PLM c-paramters
 locB=grep( "^B[[:digit:]]*$", lgitable$type )
 if( length(locB) > 0 ){
  param[locB,3]=c[locB]
 }
 colnames(param)=paste("p",1:ncol(param),sep="")
 param=cbind(giid=lgitable$giid,type=lgitable$type,ncat=ncat,param)
 nbPCA=nbparam-nsolo

 if( print > 0 ){
  cat("\n\nCalibration of Separately Estimated Item Parameters (pArameter)"
      , "\n\n")
  cat(" # of estimeats =", ne, ",  # of items =", nitems
      , ", # of forms =", nforms,"\n")
  cat(" basic form =", baseform, "\n")
  cat(" item/form inclusion pattern\n")
  ift=format(iftable);
  ift[which(ift == "NA")]="."
  print(ift,quote=0);
  cat(" # of subjects per form\n")
  print(nsubj)
  cat("\n # of items used for PCA (common items) =", nitems-nisolo,"\n")
  cat(" # of b-parameters used for PCA =", nbPCA,"\n")
  cat("\n max # of iterations for PCA =", maxiter," with eps =",eps,"\n")
  cat("  iterations required for convergence =",iterPCA,"\n")
  cat("  rmse =",rmse,"\n")
  if( print >= 3 ){
   cat("\n Item Paramter matrix analyzed by PCA\n")
   print(ParamHat, digits=2)
  }
  cat("\n Transformations Estimated\n")
  cat(" from common to each\n")
  print(uv, digits=3)
  cat("\n from each to common\n")
  cat("  (mean and std of theta for each form)\n")
  print(qr, digits=3)
  cat("\n Item Paramteres on Common Scale\n")
  Print(param, digits=2)
 }

 res=list( param=param, uv=uv, qr=qr
           , nfini=nfini, niinf=niinf, fini=fini, iinf=iinf, iftable=iftable
           , toP=toP, fromP=fromP
           , iterPCA=iterPCA, nbPCA=nbPCA, nisolo=nisolo, ParamHat=ParamHat
           , rmse=rmse
 )
 return( res )



} # end of cala






