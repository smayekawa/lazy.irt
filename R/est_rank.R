#' Estimation of the Latent Rank of LRT Model
#'
#' @param Uc Item Response Data in compressed format
#' @param U  Item Response Data
#' @param V  LRT Item Parameter Matrix:   nitem x nclass
#' @param   alpha = vector of Dirichle parameters of length nclass \cr
#' or a number A to set alpha=rep(A,nclass) \cr
#' or a negative number -A to set alpha=A*normal pdf
#' @param   vmin minimum value of V
#' @param   vmax maximum value of V
#' @param   rho initial value of probability vector of each latent class.
#' @param print = 1 to print the result
#' @param plot = 1 to plot the estimated theta distributions
#' @param title Title strings
#'
#'
#' @details
#' Note that uLRT returns the V  as the nclass x nitem data frame,
#' whereas in this function  V is defined as nitem x nclass. (sorry!)
#' \cr\cr
#' The core part of uLRT, namely, the E-step, is used in this function.
#' \cr\cr
#' The rank of a person i is defined as the location of the highest
#' posterior probability of H[i,].
#' \cr
#'
#' @return A list of: \cr
#' rank The estimated rank: nrow(Uc) x 1 \cr
#' H The posterior probability matrix: nrow(Uc) x nclass \cr
#' rho The prior probability of the class: \cr
#' alpha The prior of rho: \cr
#' nclass The number of classes \cr
#' method EAP
#'
#' @examples
#' #### In the following examples, maxiter is set to 20 which is
#' #### not large enough to obtain convergence.
#' ####
#' #
#' #
#' #
#' set.seed(1701)
#'
#' param=paramB1[c(1:3,7:9,13:15),]
#' thmin=-2; thmax=2; npoint=5
#' N=1000
#'
#' # discrete theta
#' # theta0=seq(thmin,thmax,length=npoint)
#' theta0=c(-2, -1, 0, 2, 3)
#' theta=unlist(lapply( theta0, rep, round(N/npoint) ))
#' res2 <- gendataIRT( 1, paramB1, theta=theta, compress=1 )
#' Uc=as.data.frame(res2$U)
#' ncat=res2$ncat
#' type=res2$type
#'
#' nclass=5
#' res1 <- uLRT( Uc, nclass=nclass, estrho=1, monotone=1, alpha=20
#' , maxiter=20, plot=1, print=1 )
#' V1=res1$V[,seq(2,2*res1$nitems,2)]
#' res <- est_rank( Uc=Uc, V=t(V1), print=1, alpha=-4, plot=1, title="Test" )
#'
#' @export
#'

est_rank <- function( Uc=NULL, U=NULL, V=NULL
                      , alpha=rep(1,ncol(V)), vmin=0.0001, vmax=1-vmin
                      , rho=NULL, print=0, plot=0, title=NULL ){
 # estimation of LRT rank
 # Shin-ichi Mayekawa
 # 20171107
 # as lazy.irt:  20171113
 #

 method="EAP"
 method=toupper(method)


 # error
 if( is.null(V) ){
  stop("param must be given.")
 }

 # from parameter data set
 nitem=nrow(V)
 ncat=rep(2,nitem)
 nclass=ncol(V)

 # make it a prob
 V[V<vmin]=vmin
 V[V>vmax]=vmax


 # item response data
 if( !is.null(Uc) ){
  # Uc is given.

  # error
  if( !is.null(U) ){
   cat("\nerror1(est_rank): Both U and Uc are given.\n")
   return()
  }
  if( nitem != ncol(Uc) ){
   cat("\nerror1:(est_rank) # of items in V = ", nitem
       , " but # of columns in ", Uc, " = ", ncol(Uc), "\n" )
   return()
  }
  # uncompress data
  U=dummy_expand( Uc )$U
  rm(Uc)

 }
 else if( is.null(U) ){
  cat("\nerror1(est_rank): Either U or Uc must be given.\n")
  return()
 }

 # Here, U is given.

 # make U a matrix
 if( is.data.frame(U) ) U=as.matrix(U)

 # replace NA by 0
 U[is.na(U)]=0


 N=nrow(U)

 if( sum(ncat) != ncol(U) ){
  cat("\nerror1:(est_rank) Sum of ncat in V = ", nitem
      , " but # of columns in U = ", ncol(U), "\n" )
  return()
 }


 # class distribution
 if( length(alpha) == 1 ){
  if( alpha >= 0 ){
   alpha=rep(alpha,nclass)
  }
  else{
   dd=dnorm(seq(-3,3,length=nclass))
   alpha=-alpha*dd/sum(dd)
  }
 }
 alpha[alpha<1]=1
 alpha1=alpha-1

 if( is.null(rho) ){
  temp=alpha1+1
  rho=matrix(temp/sum(temp),,1)
 }
 else if( is.vector(rho) ){
  rho=matrix(rho,,1)/sum(rho)
 }
 classname=paste("class",1:nclass,sep="")
 rownames(rho)=classname



 if( print ){
  cat("\nEstimation of Theta\n")
  Print(title)
  cat(" # of items =" , nitem, "\n")
  cat(" # of classes =" , nclass, "\n")
  cat(" # of subjects = ", N, "\n")
  cat("\n")
  cat(" estimatio method =", method, "\n")
  Print(alpha,rho)
 }

 if( method == "EAP" ){

  # assuming V is binary
  V=interleave( 1-t(V), t(V) )

  P=V

  logP=log(P)

  # Print(dim(U),dim(V),dim(P))

  # H=exp( logP%*%t(U)+log(thd) )
  # H=H%*%diag(1/colSums(H))
  # thetahat=t( t(theta)%*%H )

  ULP=U%*%t(logP)
  # For each row, make the minimum value -745 hoping that
  #  -745 < ULP < 709
  rowmin=apply(ULP,1,min)
  rowmax=apply(ULP,1,max)
  ULP=ULP-rowmin  - 600       # instead of 745 -> 720 -> 600
  # rowmin1=apply(ULP,1,min)
  # rowmax1=apply(ULP,1,max)
  # Print(rowmin,rowmax,rowmin1,rowmax1)
  ULP[ULP > 709]=709
  eULP=exp(ULP)*matrix(1,nrow(ULP))%*%t(rho)
  rseULP=rowSums(eULP)
  rseULP[rseULP < 1e-307]=1e-307
  H=eULP/rseULP
  rank=max.col(H)


 }


 freq=table(rank)


 if( print >= 1 ){
  Print(freq)
 }

 if( plot > 0 ){
  # brk=0:(nclass+1)+0.5
  # hist( rank, breaks=brk )
  title1="Frequency Distribution of the Estimated Ranks"
  if( !is.null(title) ){
   title1=paste(title,"\n",title1,sep="")
  }
  barplot(height=freq, names.arg=classname
          , main=title1, space=0, las=2)
 }


 res=named_list( rank, H, rho, alpha, nclass, method )

 return( res )


} # end of est_rank


