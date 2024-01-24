#' Calculation of Test, Item and Item Category Information Functions
#'
#' @param param Item Parameter Data Frame.
#' @param theta Discrete theta values
#' @param npoints # of discrete points for theta.
#' @param thmin Minimum value of discrete thata value.
#' @param thmax Maximum value of discrete thata value.
#' @param numderiv = 1 to use numerical first derivatives of irf.
#' @param smallP Minimum value of probability in irf and dirf functions.
#' @param print = 1 to print the summary \cr
#' = 2 to print test information functions. \cr
#' = 3 to print item information functions. \cr
#' = 4 to print item category information functions.
#' @param plot = 1 to plot test information functions\cr
#' = 2 to plot item information functions. \cr
#' = 3 to plot the locally best item and item categorie weights.
#' @param legend = 0 to skip printing legend.
#' @param maxinfo = The maximum value of item information function for plot.
#' @param debug = 1 to print intemediate result.
#'
#'
#' @return
#' A list of \cr\cr
#'  theta:    theta points\cr
#'  fromP, toP: location of each item category in item info \cr
#'  TRF: test response function (tcc) \cr
#'  icrf: item category response functions\cr
#'  dicrf: derivative of icrf \cr
#'  info:   test information function (if)\cr
#'  info_item:   item information function (iif) \cr
#'  info_item_cat:   item category information function (icif)
#'
#' @details
#' The item category information function, icif,
#'  (item response information function)
#' is defined as the 2nd derivative of \eqn{log(P_{kj}(\theta))}
#' where \eqn{P_{kj}(\theta)} is the item category response function: \cr
#' \eqn{ I_{kj}(\theta) = (P'_{kj}(\theta))^2 / P_{kj}(\theta) - }
#'  \eqn{ P''_{kj}(\theta) } \cr
#' where \eqn{P'_{kj}(\theta)} is the first derivative of \eqn{P_{kj}(\theta)}.
#'  \cr
#' The second derivatives, \eqn{ P''_{kj}(\theta) }, will be calculated
#' numerically by \code{dirt_num} using \code{lazy.mat::JacobianMat}.
#'
#' The item information function, iif, is defined as \cr
#' \eqn{  I_j(\theta) = \sum_k I_{kj}(\theta) }
#'   \eqn{= \sum_k (P'_{kj}(\theta))^2 / P_{kj}(\theta) } \cr
#' The test information function, if, is the sum of the above:\cr
#'  \eqn{  I(\theta) = \sum_j I_j(\theta) }
#'  \eqn{ = \sum_j \sum_k (P'_{kj}(\theta))^2 / P_{kj}(\theta) }.
#'
#' Note that
#' \eqn{  \sum_k P_{kj}(\theta) = 1 } and \eqn{ \sum_k P'_{kj}(\theta) =}
#' \eqn{ \sum_k P''_{kj}(\theta) = 0 } .
#'
#' Above corresponds to the information functions from \code{info_fun}
#' associated with locally best item category weights (LO).
#'
#'
#'
#' @references
#' Birnbaum, A.(1968) Some Latent Traint Models.
#' In F. M. Lord and M. R. Novick, Statistical Theories of Mental Test Scores.
#'  Reading, Mass.: Addison-Wesley. \cr
#'
#' Samejima,  F.  (2010). The General Graded Response Model. (p79-80)
#' In Nering, M. L. and Ostini, R. Eds. Handbook of Polytomous
#' Item Response Theory Models. NY, NY: Routledge
#'
#' Samejima,  F.  (1969). Estimation  of  a  latent  ability  using  a
#' response  pattern  of  graded  scores.  (Eq 6-6 in p39)
#' Psychometrika  Monographs, 34
#' (Suppl. 4).
#'
#' @examples
#' resInfo <- iif( paramA1, plot=4, print=4 )
#'
#' @export
#'

iif <- function( param, theta=NULL, npoints=31, thmin=-4, thmax=4
                       , legend=1, maxinfo=0, numderiv=0, smallP=1e-9
                       , print=1, plot=0, debug=0 ){
 # calculation of information function
 # Shin-ichi Mayekawa
 # obscore modified: 20140329
 # legend: 20150616
 # bugfix for LOW: 20150616,17
 # plot item info: 20170627
 # v_LO and wLO: 20171209
 # trf with v_LO: 20171209
 # return trf etc: 20171219dnc
 # do not use v[1,1]=1: 20171229
 # output icrf and dicrf, etc: 20180109
 # max(trf_LO) not adjusted: 20180109
 # numderiv: 20180111
 # smallP: 20180129
 # theta as an arg: 20180201
 # basic functioin: 20180209
 # lwd=2: 20180211
 # roxgen: 20180212
 # average info: 20180216cot,18cot
 # optimal -> best: 20180218cot
 # row.names=NULL: 20180221
 # output stdx_t: 20180430cot
 # roxygen2: 20221205cot
 # info_func simplified as iif: 20221206
 # _LO removed: 20221206
 # return(res): 20230202
 # roxgen2: 20230329
 #
 # Args:
 #    param    parameter data frame name
 #    weight   item and category weight data frame name
 #             see the descriptions in read.param or read.weight.
 #             If weight==NULL a set of narural weights will be used.
 #
 #    npoints  # of theta points
 #    thmin    min value of theta
 #    thmax    max value of theta
 #    print    = 1 to verbose
 #    plot     = 1 to produce several plot
 #
 # Values:
 #
 #  list( theta, fromP, toP, TRF, icrf, dicrf, info_LO, info_item_LO ,
 #   info_item_cat_LO)
 #
 #
 #
 # Needs:
 #  irf,  icrfB, icrfG,  icrfPN
 #  dirf,  dicrfB, dicrfG,  dicrfPN
 #  sumsmnw,  sumsmnw12
 #  checkparam
 #

 pdfname=as.character(deparse(substitute(param)))

 param=checkparam( param, "ALL", "iif" )
 if( is.null(param) ){
  cat("\n\nInput Item Parameter ERROR.\n\n")
  return()
 }

 # constants
 nitems=nrow(param)
 iname=param$name
 nc=ncol(param)
 ncat=param$ncat
 type=param$type


 # generate thata and prior theta dist
 if( is.null(theta) ){
  theta=seq(thmin,thmax,length.out=npoints)
 }
 thname=format(theta,digits=2)

 if( print > 0 ){
  cat("\n\nCalculation of the Information Functions \n")
  cat(" parameter data frame name =", pdfname,"\n\n")
  cat(" # of item parameters =",nitems,"\n")
  cat(" use of numerical first derivatives =", numderiv,"\n")
  cat(" range of theta = [",thmin,",",thmax,"] with", npoints
      ,"discrete points\n")
  cat("\n item parameters\n")
  print( param )
 }


 # calculate icrf and trf=cond. mean of X given theta
 temp=irf( param, theta, print=0, debug=0, plot=0, smallP=smallP )
 icrf=temp$ICRF
 TRF=temp$TRF
 fromP=temp$fromP
 toP=temp$toP
 rm(temp)

 if( debug > 0 ) Print(icrf, fromP,toP)

 # first derivative of icrf
 temp=dirf( param, theta, print=0, numderiv=numderiv
            , smallP=smallP )
 dicrf=temp$dICRF
 slope_trf=temp$dTRF
 dirf=temp$dIRF
 rm(temp)


 # information function with the locally best category weights
 dd=(dicrf^2)/icrf
 info=rowSums(dd)
 info_item=mapply( function(x,y) rowSums(dd[,x:y])
                      , fromP, toP, SIMPLIFY=FALSE )
 info_item=matrix(unlist(info_item),,nitems)
 colnames(info_item)=iname
 rownames(info_item)=thname
 names(info)=thname

 # item category info
 dd2=dicrf_num( param, theta, log=0, second=1, eps=1e-6 )
 info_item_cat=dd-dd2
 colnames(info_item_cat)=colnames(dicrf)

 # checing
 info_item2=mapply( function(x,y) rowSums((info_item_cat)[,x:y])
                      , fromP, toP, SIMPLIFY=FALSE )
 info_item2=matrix(unlist(info_item2),,nitems)

 if( print ){
  cat("\n")
  Print(max(abs(rowSums(dd2))))
  Print(max(abs(info_item2-info_item)))
 }

 if( print ){
  if( print >= 2 ){
   cat( "Test Information function\n" )
   Print(theta,info, fmt="6.3")
   if( print >= 3 ){
    cat(" Item information functions.\n")
    Print(theta, info_item, fmt="6.3")
   }
   if( print >= 4 ){
    cat(" Item Category Information Function.\n")
    Print(theta, info_item_cat, fmt="6.3")
   }
  }
 }


 if( plot > 0 ){
  # trf
  matplot( theta, info, type="l", col=1, lwd=2, lty=c(1,3)
           , main="Test Information function" )

  # plot item info
  if( plot >= 2 ){
   matplot( theta,info_item, type="l"
            , col=c("black", "red", "green", "blue"), lty=c(1:3), lwd=2
            , xlim=c(min(theta),max(theta)), ylim=c(0,max(info_item))
            , main="Item Information Functions", ylab="information")
   if( legend > 0 ){
    legend( range(theta)[1], max(info_item)-.1, iname
            , cex=0.8, col=c("black", "red", "green", "blue")
            , pch=NA, lty=c(1:3), title="legend" )
   }
  }
  if( plot >= 3 ){
   for( j in 1:(nitems) ){
    title=paste("Item Category Information Function for"
                , iname[j]," ( type = ",type[j],", ncat=",ncat[j], " )")
    pp=param[j,4:(ncat[j]+3)]
    pp[is.na(pp)]=0
    sub=paste( "param = ", paste(format(pp,digits=3)
                                 , collapse=",  ") )
    if( length(grep("^B", param$type[j])) > 0 ){
     pp=param[j,4:(ncat[j]+4)]
     pp[is.na(pp)]=0
     sub=paste("param = ", paste(format(pp,digits=3), collapse=",  "))
    }
    maxi=max(info_item_cat)
    #Print(j,theta,v[,fromP[j]:toP[j]])
    matplot(theta,info_item_cat[,fromP[j]:toP[j]], main=title, sub=sub
            , lwd=2
            , xlim=c(min(theta),max(theta)), ylim=c(0,maxi)
            , type="l", ylab="information")
    if( legend > 0 ){
     legend( range(theta)[2]-1, range(info_item_cat)[2], 0:(ncat[j]-1)
             , cex=0.8, col=c("black", "red", "green", "blue")
             , pch=NA, lty=c(1:3), title="legend" )
    }
   }
  }
 }

 res=named_list( theta, fromP, toP, TRF, icrf, dicrf
                 , info, info_item, info_item_cat )

 return( res )

} # end of iif




res=iif( paramA1[8,], npoints=21, plot=3, print=1 )





comments(
 '

paramS0=paramS1[1:2,]
paramS0$type=c("Bn","Bn3")
paramS0$name=c("Q1n","Q2n")
paramS0=rbind(paramS1[1:2,],paramS0)
paramS0$p1=1
theta=seq(-4,4,length=121)
resInfo=info_func( paramS0, theta, plot=4, print=1 )





resInfo=info_func( paramS1, npoints=15, plot=1, print=1 )
resInfo=info_func( paramS1, weight=weightS11, npoints=15, plot=1, print=0 )
resInfo=info_func( paramS1, weight=weightS12, npoints=15, plot=1, print=0 )


resInfo=info_func( paramS2, plot=3, print=3, npoints=151 )



param=paramB2[c(1:3,7:9,13:15),]
param$p3=0.5
resInfo=info_func( param, plot=3, print=3, npoints=151 )


param=paramS2[c(1,3,5,7),]
param$p1=1
resInfo=info_func( param, plot=3, print=3, npoints=51 )


param=paramS2[c(1,3,5,7),]
param$p1=1
weight=weightS2[c(1,3,5,7),]
resInfo=info_func( param, weight, plot=3, print=3, npoints=51 )

param=paramS2[c(1,3,5,7),]
param$p1=1
weight=weightS2[c(1,3,5,7),]
weight[4,7]=4
resInfo=info_func( param, weight, plot=3, print=3, npoints=151 )


')








