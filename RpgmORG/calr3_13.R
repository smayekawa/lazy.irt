#' Calibration of Independently Estimated Sets of Item Parameters
#'
#' This function calibrates (equates) independently estimated sets of item
#' parameters (with common linking items) by minimizing the LS criterion 
#' defined in terms of item response functions (irf) \cr
#' Japanese help file: (\link[lazy.irt]{calr_JPH})
#' 
#' @param indata   data frame containing at least following: \cr
#'             giid  gfid  type ncat  a  p1 p2 .... \cr
#'  indata can be created from a parameter data frame
#'  by renaming name as giid and adding gfid.
#'  It is assumed that indata is sorted by gfid.
#' @param baseform  form number or 0
#' @param nsubj a vector consisting of the # of subjects in each form or 1.
#' \cr To be used when calculating the mean of theta of the combined group.
#' @param npoints # of discrete theta points to be used.
#' @param thmin Minimum value of theta points.
#' @param thmax Maximum value of theta points.
#' @param maxiter  max of iterations
#' @param eps Criterion of convergence for the relative improvement of
#' rss value.
#' @param print       > 0 to print results
#' @param debug = 1 to print intermediate result
#'
#' @details
#'
#'  This program finds q, r and item param such taht \cr
#'     P(theta_g | est_param_g) = P( q[g] + r[g] theta | param )  + error \cr
#'      j=1,2,..., # of items and  g=1,2,..., # of forms, \cr
#'  where param is the item parameters on the common scale. \cr
#' \cr
#'  First, those items which are included in at least two forms will be
#'  picked up and analyzed. \cr
#'  Then, the item parameters for the solo items will be converted
#'  using estimated q and r. \cr
#' \cr
#'  Subscript j is for item, g for form, q is for theta point.  \cr
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
#'  iter # of iteration required for convergence \cr
#'  rss # LS criterion minimized \cr
#'  rssimpr Final relative improvement of rss. \cr
#'  eps Convergence criterion for rssimpr \cr
#'  rmse rmse calculated from rss \cr
#'  nisolo # of items that are included in only one form. \cr
#'  gradI Final gradients of common item parameters. \cr
#'  gradF Final gredients or transformation constants. \cr
#'
#' @references
#' Sayaka Arai, Tomohiro Ohtani, & Shin-ichi Mayekawa (2013)
#' Simultaneous Equating of Separately Calibrated Polytomous IRT Test Items.
#' IMPS 2013, CITO, Arnhem, The Netherland. 23-26 July 2013 \cr
#'
#' Shin-ichi Mayekawa (2012)
#' Simultaneous equating of separately calibrated item parameters
#' under the common item design.
#' COMPSTAT 2012, Amathus Beach Hotel, Limassol, Cyprus. 27-31 August 2012\cr
#'
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
#' res=calr( indata=paramCal1, baseform=1, print=2  )
#'
#' # Mixture of Item Types
#' res=calr( indata=paramCal2, baseform=1, print=2  )
#'
#' @export
#'

calr <- function( indata, baseform=0, nsubj=0
                  , npoints=21, thmin=-3, thmax=3
                  , maxiter=100, eps=1e-5
                  , print=1, debug=0  ){
 # calr: Calibration of Separately Estimated Item Parameters (pRobability)
 # Shin-ichi Mayekawa
 # Fortran version for binary items: 870903 -- 080727
 # 120227,29,0301
 # solo items: 120303,04
 # baseform=0: 120304
 # when PCA returns negative v: 120908
 # cala -> calr: 120908 cot, 11
 # type Bxx items: 120921,22
 # initial by function cala.core: 120922
 # full matrix RSSI: 120925
 # work on non-solo items only: 120925
 # cala.core renamed as cala: 120926
 # calr.core renamed as calr: 120926
 # vertically stacked icrfhatj: 120927
 # final version which uses NLM w/o grad: 120928
 # GN: 121004,05
 # Jacobian by hand: 121007 cot j,08
 # get rid of NLM for all: 121008,10,11
 # full iftable: 121012,13
 # dirf_p: 121013
 # bugfix: 121029,30,31,1102
 # type="P" items by "PN": 121119
 # type="P" items as the are: 121120
 # sort param back to the original order iwth solo items: 130206
 # size of gradI0: 130629
 # constraints on c: 20150303
 # bugfix: 20150303
 # rmse: 20150615
 # nsubj and nsb: 20161121
 # roxygen2: 20231027,28cot
 #
 #
 # Args:
 #
 #   indata      data frame containing at least following
 #               giid (or name)  gfid  type ncat  a  p1 p2 ....
 #
 #               indata can be created from a parameter data frame
 #               by renaming name as giid and adding gfid
 #
 #               It is assumed that indata is sorted by gfid.
 #   baseform    base form number or 0
 #   nsubj       not used
 #   npoints     # of theta points in [thmin,thmax]
 #
 #   maxiter     max # of iterations
 #   eps         eps for the relative improvement of rss
 #   print
 #
 #
 #
 #  This program finds q, r and item param such taht
 #     P(theta_g | est_param_g) = P( q[g] + r[g] theta | param )  + error
 #  where param is the item parameters.
 #
 #  First, those items which are included in at least two forms will be
 #  picked up and analyzed.
 #  Then, the item parameters for the solo items will be converted
 #  using estimated q and r.
 #
 #  Subscript j is for item, g for form, q is for theta point
 #
 #
 # Requires:
 #  irf, icrfB, icrfG, icrfPN, icrfP
 #  dirf, dicrfB, dicrfG, dicrfPN, dicrfP, dirf_p
 #  Mh2Mv
 #  cala, princ
 #
 #

 RSSI <- function( paramnumj, paramj=NULL, icrfhatj, thg, thd ){
  # calculation of RSS for item j
  # 120908,22,24
  # full matrix version: 120924
  # vertically stacked icrfhatj: 120927
  # thg as argument: 120927
  # use of irf function restored: 120928
  # # of parameters reduced: 120928
  #
  # Args:
  #
  #  paramnumj    numerical part of item paramter data frame
  #               for item j (common scale)
  #                i.e.  a and b parameters (no ncat here)
  #               RSSI will be minimized w.r.t this vector
  #
  #  paramj       item paramter data frame for item j (common scale)
  #               This will be haded to IRF module
  #
  #  icrfhatj     npoints x nfini*(ncat[j]-1)  icrf data for item j
  #
  #  thg          qr[fini[[j]],1] %x% rep(1,length(theta))
  #                  + qr[fini[[j]],2] %x% theta
  #  thd          normalized pdf of N(0,1) at theta
  #                should be rep(thd,nfini)
  #

  # Print("**RSSI**", paramnumj, j, paramcharj,nfinij, finij, theta, thd)

  # replace the item paramter part of paramj
  paramj[,4:(4+length(paramnumj)-1)]=paramnumj

  icrfj=irf(paramj,thg,zero=0,print=0)$ICRF

  rssj=icrfhatj-icrfj
  rssj=sum(thd*rssj*rssj)

  return( rssj )

 } # end of RSSI


 RSSF <- function( qrg, icrfhatg, paramg, theta, thd ){
  # calculation of RSS for from g
  # 120909 cot, 24
  # # of parameters reduced: 120928
  #
  # Args:
  #
  #  qrg          c( q(g), r(g) )
  #               RSSF will be minimized w.r.t this vector
  #
  #  icrfhatg     npoints x sum(ncat[g]-1)  icrf data for form g
  #
  #  paramg       item paramter data frame contained in form g
  #
  #  theta        npoints x 1   theta values
  #  thd          normalized pdf of N(0,1) at theta
  #
  #

  # Print("**RSSF**", qrg, g, icrfhatg, paramg, theta, thd)
  thg=qrg[1]+qrg[2]*theta
  icrfg=irf(paramg,thg,zero=0,print=0)$ICRF
  # Print(g,icrfhatg,icrfg)
  rssg=icrfhatg-icrfg
  rssg=sum(thd*rssg*rssg)

  return( rssg )

 } # end of RSSF


 comments('
          # convert to PN format
          P2PN=0
          if( any(indata$type == "P") ){
          indata=convP2PN(indata)
          P2PN=1
          # Print(indata)
          }
          ') # end of comments


 # const
 ne=nrow(indata)

 # delete exessive columns
 locnotna=which( !apply( indata, 2, function(x) all(is.na(x)) ) )
 indata=indata[,locnotna]

 # add seqnum
 if( length(indata$seqnum) == 0 )
  indata=cbind(seqnum=1:ne,indata)

 # set the c-param of type B2 items equal to zero if NA
 indata$p3[which(indata$type == "B2" & is.na(indata$p3))]=0

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


 # add local id (liid and lfid) to indata by merging
 # liid to giid table:   liid=1,2,...,max_unique_giid
 # unique item names in the original order
 lgitable=data.frame( liid=as.numeric(1:nitems), giid=ugiid
                      , stringsAsFactors=0  )
 if( length(indata$liid) == 0 ){
  # merge liid
  indata=merge(lgitable, indata, by="giid",sort=0) # sort=0 does not work.
 }

 # lfid to gfid table:   lfid=1,2,...,nforms
 lgftable=data.frame( lfid=as.numeric(1:nforms), gfid=ugfid
                      , stringsAsFactors=0  )
 if( length(indata$lfid) == 0 ){
  # merge lfid
  indata=merge(lgftable,indata,by="gfid",sort=0) # sort=0 does not work.
 }

 # sort back to the original order
 indata=indata[order(indata$seqnum),]
 rownames(indata)=c(1:nrow(indata))


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

 # # of categories of each item (common and solo): pick up the first element
 ncat=indata$ncat[ sapply( loci,"[[",1 ) ]

 # original fitable for all the items (common and solo)
 iftable0=matrix(NA,nitems,nforms)
 for( j in 1:nitems ){
  iftable0[ j,indata$lfid[ loci[[j]] ] ]=1
 }
 iftable0=cbind(iftable0,nfini)
 rownames(iftable0)=ugiid
 colnames(iftable0)=c(ugfid,"form_count")
 iinf0=NULL
 for( g in 1:nforms )
  iinf0=append(iinf0,list(which(iftable0[,g]==1)))

 # merge nfini and keep the subset in which nfini > 1
 # cannot add nfini to indata because later its column will be
 # refferred by number.
 nfinitable=data.frame( liid=as.numeric(1:nitems), nfini=nfini
                        , stringsAsFactors=0  )
 temp=merge(nfinitable, indata, by="liid")
 # sort back to the original order
 temp=temp[order(temp$seqnum),]


 # store the  common items in indata1
 locsolo=which(temp$nfini == 1)
 loccommon=which(temp$nfini > 1)
 indata_solo=indata[locsolo,]
 indata1=indata[loccommon,]
 ne1=nrow(indata1)
 nisolo=nrow(indata_solo)
 rm(temp)

 nitems1=sum(nfini > 1)

 if( debug ) Print("**** before cala", indata1)
 # initial by cala
 res=cala( indata1, baseform=baseform, maxiter=100, eps=1e-8
           , print=min(0,print-3), debug=debug  )
 uv=res$uv; qr=res$qr
 param=res$param  # common item parameters
 ncat1=param$ncat # # of categories of the common items
 type=param$type  # item type of the common items
 nfini=res$nfini  # # of forms including each of the common items
 niinf=res$niinf  # # of common items included in each form
 fini=res$fini    # list of formms including each of the common items
 iinf=res$iinf    # list of common items included in each form
 iftable=res$iftable;
 iterPCA=res$iterPCA; nbPCA=res$nbPCA
 ParamHat=res$ParamHat; rmse=res$rmse
 rm(res)



 if( print > 0 ){
  cat("\n\nCalibration of Separately Estimated Item Parameters"
      , "  (pRobability)", "\n")
  cat(" # of estimeats =", ne, ",  # of items =", nitems
      , ", # of forms =", nforms,"\n")
  cat(" # of items to be analzed (common items) =", nitems1, "\n")
  cat(" basic form =", baseform, "\n")
  cat("\n item/form inclusion pattern\n")
  ift=format(iftable0)
  ift[which(ift == "NA")]="."
  print(ift,quote=0); cat("\n")
  cat(" # of subjects per form\n")
  print(nsubj)
  cat(" # of theta points =", npoints, " in (", thmin, ",", thmax,")","\n")
  cat("\n max # of iterations =", maxiter," with eps =",eps,"\n")
  cat(" print control =", print,"\n")
  if( print >= 2 ){
   cat("\nInitial Estimates by CAL-A\n")
   cat("\n # of items used for PCA (common items) =", nitems-nisolo,"\n")
   cat(" # of b-parameters used for PCA =", nbPCA,"\n")
   cat("  iterations required for convergence =",iterPCA,"\n")
   cat("\n Initial Transformations\n")
   cat(" from common to each\n")
   print(uv, digits=3)
   cat("\n from each to common\n")
   cat("  (mean and std of theta for each form)\n")
   print(qr, digits=3)
   cat("\n Initial Item Paramteres on Common Scale\n")
   param0=param
   rownames(param0)=NULL
   Print(param0, digits=2)
  }
 }

 ######################### main body starts here. #######################



 # testing
 # param$p1=1
 # param$p2=-1
 # param$p3=0.2


 rm( lgftable )


 # theta distribution
 theta=seq(thmin,thmax,length.out=npoints)
 thd=exp(-0.5*theta^2); thd=thd/sum(thd)


 if( debug ) Print(" in calr, after cala", indata1)


 # keep necessary variables
 locparam=which( colnames(indata1) %in% c("gfid", "giid", "type", "ncat") )
 locparam=c(locparam,grep("^p[[:digit:]]",colnames(indata1)))


 # create data icrf (icrfFhat) sorted by form
 # Do not sort the data. Use the original order.
 indata2=indata1
 res=irf(indata2[,locparam],print=0, zero=0
         ,  npoints=npoints, thmin=thmin, thmax=thmax)
 icrfFhat=res$ICRF

 # index in icrfFhat
 fromicrfF=rep(0,nforms); toicrfF=rep(0,nforms)
 fromicrfF[1]=1; toicrfF[1]=sum(ncat1[iinf[[1]]])-niinf[1]
 for( g in 2:nforms ){
  fromicrfF[g]=toicrfF[g-1]+1
  toicrfF[g]=fromicrfF[g]+sum(ncat1[iinf[[g]]])-niinf[g]-1
 }

 if( debug ) Print(fromicrfF,toicrfF,niinf,ncat, ncat1)
 # Here, icrf (data) of form g is stored in
 #   icrfFhat[,fromicrfF[g]:toicrfF[g]]
 #  where  toicrfF[g]-fromicrfF[g]=sum(ncat[iinf[[g]]])-niinf[g]
 # The item ids included form g is stored in iinf[[g]].
 #



 # create data icrf (icrfIhat) sorted by item
 # sort indata by item(liid order)
 indata2=indata1
 od=order(indata2$liid,indata2$lfid)
 indata2=indata2[od,]

 if( debug ) Print(indata2)
 res=irf(indata2[,locparam],print=0, zero=0
         ,  npoints=npoints, thmin=thmin, thmax=thmax)
 icrfIhat=res$ICRF; fromPI=res$fromP; toPI=res$toP
 if( debug ) Print(theta, icrfIhat)
 #Print(fromP, toP, ncat1, digits=2)

 # index in icrfIhat
 fromicrfI=rep(0,nitems1); toicrfI=rep(0,nitems1)
 fromicrfI[1]=1; toicrfI[1]=nfini[1]*(ncat1[1]-1)
 for( j in 2:nitems1 ){
  fromicrfI[j]=toicrfI[j-1]+1
  toicrfI[j]=fromicrfI[j]+nfini[j]*(ncat1[j]-1)-1
 }
 if( debug ) Print(fromicrfI,toicrfI,nfini,ncat,ncat1)
 # Here, icrfhat (data) of item j is stored in
 #   icrfIhat[,fromicrfI[j]:toicrfI[j]]
 #  where  toicrfI[j]-fromicrfI[j]=nfini*(ncat[j]-1)-1.
 # The form ids including item j is stored in fini[[j]].
 #

 rm(indata2)

 # initial RSS
 rssI=0
 for( j in 1:nitems1 ){
  icrfIhatj=icrfIhat[,fromicrfI[j]:toicrfI[j],drop=F]
  icrfIhatj=Mh2Mv( icrfIhatj, npoints, ncat1[j]-1, nfini[j] )
  paramj=param[j,]
  paramnumj=as.numeric( paramj[,4:(4+ncat1[j]-1)] )
  if( type[j] == "B3" ) paramnumj=as.numeric( paramj[,4:(4+ncat1[j])] )
  thg=matrix( qr[fini[[j]],1]%x%rep(1,length(theta))+qr[fini[[j]],2]%x%theta
              ,,1 )
  rssjnew=RSSI( paramnumj, paramj=paramj, icrfhatj=icrfIhatj
                , thg, thd )
  rssI=rssI+rssjnew
 } # end of j loop

 rssF=0
 for( g in 1:nforms ){
  icrfFhatg=icrfFhat[,fromicrfF[g]:toicrfF[g],drop=F]
  paramg=param[iinf[[g]],]
  qrg=qr[g,]
  rssgnew=RSSF( qrg, icrfhatg=icrfFhatg, paramg=paramg, theta=theta
                , thd=thd )
  rssF=rssF+rssgnew
 } # end of g loop

 if( print >= 2 ){
  Print("initial RSS values", rssI, rssF)
 }


 # main iteration

 rssp=99999;
 rssimp=1e9

 # matrix to store the gradients: dRSS/dparam
 gradI=matrix(NA,nitems1,max(ncat1))
 if( length(grep("^B[[:digit:]]*$",type)) == nitems1 )
  gradI=matrix(NA,nitems1,max(ncat1)+1)
 rownames(gradI)=rownames(param)
 colnames(gradI)=colnames(param)[4:ncol(param)]
 gradF=matrix(NA,nforms,2)
 rownames(gradF)=rownames(qr)
 colnames(gradF)=colnames(qr)




 for( llll in 1:maxiter ){

  # update item parameter
  rssI=0
  for( j in 1:nitems1 ){

   # data icrf for item j: vertically stacked matrix
   #  (nforms[j]*npoints) x (ncat[j]-1)
   #  ( not yet vectorized )
   icrfIhatj=icrfIhat[,fromicrfI[j]:toicrfI[j],drop=F]
   icrfIhatj=Mh2Mv( icrfIhatj, npoints, ncat1[j]-1, nfini[j] )

   paramj=param[j,] # 1:name, 2:type, 3:ncat, 4:a, 5:b1, 6:b2, .....
   ncatj=ncat1[j]
   np=ncatj
   if( type[j] == "B3" ) np=np+1

   paramnumj=as.numeric( paramj[,4:(4+ncat1[j]-1)] )
   if( type[j] == "B3" ) paramnumj=as.numeric( paramj[,4:(4+ncat1[j])] )

   a=paramj$p1
   b=as.matrix(paramj[,5:(5+ncatj-2)],1)
   if( type[j] == "G" ) b=cbind(b,9999)
   c=paramj$p3


   # this is thetas for all nfini[j] groups stacked vertically
   thg=matrix( qr[fini[[j]],1]%x%rep(1,length(theta))+qr[fini[[j]],2]%x%theta
               ,,1 )
   lenthg=nrow(thg)

   # calculate Pj and Jacobian
   temp=dirf_p( paramj, theta=thg )
   Jack=temp$Jack
   Pj=temp$Pj

   # gradient and Hessian
   g=-2*t(Jack)%*%(thd*matrix(icrfIhatj-Pj,,1))
   H=2*t(Jack)%*%(thd*Jack)

   # descent direction by GN
   d=solve(H)%*%g

   #Print(j,ncat[j],nfini[j],"/",Jack)
   #Print(g,d,H)

   rssj=RSSI( paramnumj, paramj, icrfIhatj, thg, thd )

   # update parameters along the direction d
   step=1
   ok=0
   for( lllll in 1:20 ){
    pnew=paramnumj-step*d

    # constraints
    if( pnew[1] <= 0 ) break
    pnew[pnew < -10]=-10
    pnew[pnew >  10]=10
    if( type[j] == "B3" && pnew[3] < .0 ) pnew[3]=0
    if( type[j] == "B3" && pnew[3] > .7 ) pnew[3]=.7

    rssjnew=RSSI( pnew, paramj, icrfIhatj, thg, thd )
    if( rssjnew <= rssj ){
     ok=1
     break
    }
    step=step/2
   } # end of halving
   if( ok ){
    # update parameters
    param[j,4:(4+length(pnew)-1)]=pnew
    rssj=rssjnew
   }
   rssI=rssI+rssj

   gradI[j,1:np]=g

  } # end of j loop

  if( print >= 4 ){
   Print(llll, rssI, digits=5)
  }

  maxagradI=max(abs(gradI),na.rm=1)




  # update transformation parameter
  rssF=0
  for( g in 1:nforms ){
   icrfFhatg=icrfFhat[,fromicrfF[g]:toicrfF[g],drop=F]
   paramg=param[iinf[[g]],]
   qrg=qr[g,]
   thg=qrg[1]+qrg[2]*theta
   icrfFg=irf( paramg, thg, zero=0, print=0 )$ICRF

   temp=matrix( dirf( paramg, thg, zero=0, print=0, debug=0 )$dICRF, , 1 )
   Jack=cbind( temp, theta*temp )
   # Print(g, niinf[g], paramg)

   # gradient and Hessian
   gg=-2*t(Jack)%*%(thd*matrix(icrfFhatg-icrfFg,,1))
   H=2*t(Jack)%*%(thd*Jack)

   # descent direction by GN
   d=solve(H)%*%gg

   #Print(j,ncat[j],nfini[j],"/",Jack)
   #Print(g,d,H)

   rssg=RSSF( qrg, icrfFhatg, paramg, theta, thd )

   # update parameters along the direction d
   step=1
   ok=0
   for( lllll in 1:20 ){
    pnew=qrg-step*d
    rssgnew=RSSF( pnew, icrfFhatg, paramg, theta, thd )
    if( rssgnew <= rssg ){
     ok=1
     break
    }
    step=step/2
   } # end of halving
   if( ok ){
    # update parameters
    qr[g,]=pnew
    rssg=rssgnew
   }
   rssF=rssF+rssg

   gradF[g,]=gg


  } # end of g loop
  v=1/qr[,2];  u=-qr[,1]/qr[,2];  uv=cbind(u,v)

  maxagradF=max(abs(gradF),na.rm=1)


  if( print >= 4 ){
   Print(llll, rssF, rssI-rssF, digits=5)
  }

  # check convergence
  rss=rssF
  rssimp=rssp-rss
  rssimpr=rssimp/rssp
  if( print >= 2 ){
   cat("llll=",llll
       ,", rss =", format(rss,nsmall=6, width=10)
       , ", rssimp =", format(rssimp,nsmall=6, scientific=1)
       , ", rssimpr =", format(rssimpr,digits=6, scientific=1)
       , "\n")
  }
  if( rssimpr <= eps ) break

  # next iteration
  rssp=rss


  # adjust parameter to base form



 } # end of llll loop

 iter=llll


 if( debug ) Print( "in calr: after iteration" )


 # adjust parameter to base form
 locB=grep("^B[[:digit:]]*$",param$type)
 q=qr[,1]; r=qr[,2]
 a=param[,4]; b=param[,5:ncol(param)]
 cparam=param[locB,6]

 if( baseform != 0 ){
  qb=q[baseform]; rb=r[baseform]
  b=(b-qb)/rb;   a=a*rb
  q=(q-qb)/rb;   r=r/rb
  v=1/r;   u=-q/r
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

 # final result: common item part only
 uv=cbind(u,v); qr=cbind(q,r)
 param[,4]=a; param[,5:ncol(param)]=b
 param[locB,6]=cparam
 colnames(uv)=c("u","v"); rownames(uv)=ugfid
 colnames(qr)=c("q","r"); rownames(qr)=ugfid

 if( debug ) Print( "in calr: before solo",indata_solo, param )

 # add solo items
 if( nisolo > 0 ){
  # Print(colnames(indata_solo), grep("^p[[:digit:]]$",colnames(indata_solo)))
  loc=grep("^p[[:digit:]]$",colnames(indata_solo))
  param_solo=cbind( indata_solo$giid, indata_solo$type, indata_solo$ncat
                    , indata_solo[,loc] )
  colnames(param_solo)=c( "giid","type","ncat"
                          , colnames(param_solo)[4:ncol(param_solo)] )
  rownames(param_solo)=param_solo$giid
  #Print(param, indata_solo)

  fini_solo=NULL
  for( j in 1:nisolo){   # nrow(param_solo)
   g=indata_solo$lfid[j]
   loc=which( colnames(param_solo) %in% paste("p",2:ncat[j],sep="") )
   param_solo$p1[j]=param_solo$p1[j]/qr[g,2]
   param_solo[j,loc]=qr[g,1]+param_solo[j,loc]*qr[g,2]
   if( length(grep("^B[[:digit:]]*$",param_solo$type[j])) > 0 )
    param_solo[j,6]=indata_solo$p3[j]

   fini_solo=append(fini_solo,list(g))
  }

  #  param=rbind(param, param_solo)

  param=merge(param,param_solo,all=1)
  nfini=rbind(nfini,matrix(1,nisolo))
  fini=append(fini,fini_solo)
  gradI0=matrix(NA,nrow(param),ncol(param)-3) #*****************************
  gradI0[1:nitems1,1:ncol(gradI)]=gradI
  gradI=gradI0
  gradI[is.na(param[,4:ncol(param)])]=NA
 } # solo items

 # make giid and type character, not factor
 param$giid=as.character(param$giid)
 param$type=as.character(param$type)

 rownames(nfini)=param$giid
 rownames(gradI)=param$giid

 # Print(indata_solo, param_solo, nfini, fini,iinf)

 # sort back by liid
 param=cbind(param,temp=1:nitems)
 param=merge(param,lgitable,by="giid", all=1)
 rownames(param)=NULL
 od=order(param$liid)
 param=param[od,]
 temp=param$temp
 param=param[,1:(ncol(param)-2)]
 nfini=nfini[temp]
 fini=fini[temp]
 iinf=iinf0
 gradI=gradI[nfini > 1,]
 param0=cbind(param,form_count=nfini)

 rmse=sqrt(rss/npoints/sum(ncat1-1))

 comments('
          # convert back to P format
          if( P2PN == 1 ){
          param=convPN2P(param)
          }
          ') # end of comments


 if( print > 0 ){
  cat("\n\nCalibration of Separately Estimated Item Parameters"
      , "  (pRobability)", "\n\n")
  cat(" # of estimeats =", ne, ",  # of items =", nitems
      , ", # of forms =", nforms,"\n")
  cat(" # of items analzed (common items) =", nitems1, "\n")
  cat(" basic form =", baseform, "\n")
  cat("\n item/form inclusion pattern\n")
  print(ift,quote=0); cat("\n")
  cat(" # of subjects per form\n")
  print(nsubj)
  cat("\n # of theta points =", npoints, " in (", thmin, ",", thmax,")","\n")
  cat(" \n print control =", print,"\n")
  cat("\n max # of iterations =", maxiter," with eps =",eps,"\n\n")
  cat("  # of iterations before termination =",iter,"\n")
  cat("  final RSS value =", rss," with relative imp. =", rssimpr," \n")
  cat("  final rmse =", rmse,"\n")
  cat("  max abs gradient (item and form) =",maxagradI," and"
      , maxagradF,"\n")
  cat("\n Transformations Estimated\n")
  cat(" from common to each\n")
  print(uv, fmt="8.4")
  cat("\n from each to common\n")
  cat("  (mean and std of theta for each form)\n")
  print(qr, fmt="8.4")
  cat("\n Item Paramteres on Common Scale\n")
  print(param0)
  if( print >= 2 ){
   cat(" final gradient values\n")
   Print(gradI, fmt="12.8")
   Print(gradF, fmt="12.8")
  }
 }



 return( list( param=param,uv=uv,qr=qr
               , nfini=nfini, niinf=niinf, fini=fini
               , iinf=iinf, iftable=iftable0
               , iter=iter, rss=rss, rssimpr=rssimpr, eps=eps, rmse=rmse
               , nisolo=nisolo
               , gradI=gradI, gradF=gradF ) )



} # end of calr










res=calr( indata=paramCal1, baseform=1, print=2  )










comments(
  '


# must follow gen_caldata4.R


etime(0)
res=calr( indata=paramCal11, baseform=1, print=2  )
etime(1)


  '
)
