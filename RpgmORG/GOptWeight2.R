#' Estimation of the Globally Optimal Item Category Weights
#'
#' This program estimates a set of globally optimal item category weights
#' which maximizes the expected test imformation function.
#'
#' @param param Item paramter data frame
#' @param i.weight A vector of item weights: w
#' @param c.weigth A vector of item category weights: v
#' @param npoints # of discrete theta points
#' @param thmin Minimum value of theta points
#' @param thmax Maximum value of theta points
#' @param th.m Mean of theta dist
#' @param th.sd Standard deviation of theta dist
#' @param method = "nlm"
#' @param L max # of iterations
#' @param eps eps for gradient
#' @param print = 1 to print result
#'
#'
#' @examples
#' res1 <- GOptWeight( paramS1 )
#' create_weight_df( paramS1, vvec=res1 )
#'
#' res2 <- GOptWeight( paramS2 )
#' create_weight_df( paramS2, vvec=res2 )
#'
#' @author Contributed by Dr. Sayaka Arai of DNC.
#'
#' @export
#'

#
# weight.df replaced by create_weight_df in the examples above: 20230716cot
#

GOptWeight <- function( param, i.weight=NULL, c.weight=NULL
                       , npoints=31, thmin=-4, thmax=4, th.m=0, th.sd=1
                       , method="nlm", L=100, eps=1e-6, print=0 ){

 nitems <- nrow(param)
 iname <- param$name
 ncat <- param$ncat
 Q <- npoints

 # irf & dirf
 temp1 <-  irf(param, weight=NULL, npoints=Q, thmin=thmin, thmax=thmax
               , print=print)
 temp2 <- dirf(param, weight=NULL, npoints=Q, thmin=thmin, thmax=thmax
               , print=print)
 Pmat  <- t(temp1$ICRF)    # sum(ncat) x Q matrix
 dPmat <- t(temp2$dICRF)
 toP   <- temp1$toP
 fromP <- temp1$fromP

 # weight(simple weight)
 if (is.null(i.weight)) {
  i.weight <- rep(1, nitems)
 }
 if (is.null(c.weight)) {
  c.weight <- make.simple.weight(ncat, nitems)
 }
 ws <- c.weight

 # theta
 theta <- seq(thmin, thmax, length.out=npoints)
 th.dist <- dnorm(theta, m=th.m, sd=th.sd)/ sum(dnorm(theta, m=th.m, sd=th.sd))

 m <- 0; mdif <- 100

 if( method != "nlm" ){
  for (l in 1:L){
   if (mdif < 1e-06) break
   if( print ) cat("--------------------------------- \n")
   for (J in 1:nitems){
    wsj <- ws[(fromP[J]+1):toP[J]]
    ans <- maxNR(fn = EIj(param, Pmat, dPmat, J, ws, th.dist),
                 grad = EIjD1(param,Pmat, dPmat, J, ws, Q, th.dist),
                 hess = EIjD2(param,Pmat, dPmat, J, ws, Q, th.dist),
                 start = wsj, print.level=0, tol = eps, iterlim = L)
    ws[(fromP[J]+1):toP[J]] <- ans$estimate
    if( print )
     cat("item=",iname[J]," maximum=",ans$maximum," code=",ans$code, "\n")
   }
   mdif <- abs((ans$maximum-m)/m)
   m <- ans$maximum
   if( print ){
    cat("--------------------------------- \n")
    cat("Iteration=",l,",   maximum=",ans$maximum,",   mdif=",mdif, "\n")
    Print(ws)
   }
  }
 }
 else{
  for (l in 1:L){
   if (mdif < 1e-06) break
   if( print ) cat("--------------------------------- \n")
   for (J in 1:nitems){
    wsj <- ws[(fromP[J]+1):toP[J]]
    ans <- nlm(f = mEIj(param, Pmat, dPmat, J, ws, th.dist)
               , p = wsj, print.level=0, gradtol = eps, iterlim = L)
    ws[(fromP[J]+1):toP[J]] <- ans$estimate
    if( print )
     cat("item=",iname[J]," maximum=",-ans$minimum," code=",ans$code, "\n")
   }
   mdif <- abs((-ans$minimum-m)/m)
   m <- -ans$minimum
   if( print ){
    cat("--------------------------------- \n")
    cat("Iteration=",l,",   maximum=",-ans$minimum,",   mdif=",mdif, "\n")
    Print(ws)
   }
  }
 }

 return(ws)

} # end of GOptWeight




 #-----------------------------------------
 # make simple weights
 #-----------------------------------------
 make.simple.weight <- function(ncat, nitems){
  ws <- numeric(0)
  for (j in 1:nitems){
   ws <- append(ws, c(0:(ncat[j]-1)))
  }
  return(ws)
 }



 #-----------------------------------------
 ### Small parts of functions           ###
 #-----------------------------------------
 # numerator of Eq(66) squared for item J
 SumdP <- function(param, dPmat, J, ws, wj){
  ncat <- param$ncat
  toP <- cumsum(ncat)
  fromP <- toP - ncat + 1

  w <- ws
  w[fromP[J]:toP[J]] <- c(0, wj[1:(ncat[J]-1)])
  return(w %*% dPmat)
 }

 # denominator of Eq(66) for item J
 VarX <-function(param, Pmat, J, ws, wj){
  nitems <- nrow(param)
  ncat <- param$ncat
  toP <- cumsum(ncat)
  fromP <- toP - ncat + 1

  w <- ws
  w[fromP[J]:toP[J]] <- c(0, wj[1:(ncat[J]-1)])
  frac2 <- 0
  for (j in 1:nitems){
   w1 <- w[fromP[j]:toP[j]]
   P1 <- Pmat[fromP[j]:toP[j], ]
   frac2 <- frac2 + (w1^2) %*% P1 - (w1 %*% P1)^2
  }
  return(frac2)
 }

 #==========================================
 # Unknown category weight
 #==========================================
 #-----------------------------------------
 # Expected Information (function of wj) Eq(66)
 EIj <- function(param, Pmat, dPmat, J, ws, th.dist){
  return(
   function(wj){
    ((SumdP(param, dPmat, J, ws, wj)^2) / VarX(param, Pmat, J, ws, wj)) %*%
     th.dist
   }
  )
 }

 mEIj <- function(param, Pmat, dPmat, J, ws, th.dist){
  return(
   function(wj){
    res=
     -((SumdP(param, dPmat, J, ws, wj)^2) / VarX(param, Pmat, J, ws, wj)) %*%
     th.dist
    return( res )
   }
  )
 }


 #-----------------------------------------
 # Expected Information (function of wj) with d1 Eq(68)
 EIjD1 <- function(param, Pmat, dPmat, J, ws, Q, th.dist){
  return(function(wj){

   ncat <- param$ncat
   toP <- cumsum(ncat)
   fromP <- toP - ncat + 1

   k <- ncat[J] - 1  # length of wj
   d1 <- numeric(k)

   frac1 <- as.vector(SumdP(param, dPmat, J, ws, wj))  #(1, Q)
   frac2 <- as.vector(VarX(param, Pmat, J, ws, wj))    #(1, Q)
   #term1
   dPJk <- matrix(dPmat[(fromP[J]+1):toP[J], ], ncol=Q) #(k,Q)matrix for itemJ
   term1 <- dPJk %*% diag(frac1) %*% solve(diag(frac2))
   #term2
   PJ  <- Pmat[(fromP[J]  ):toP[J], ]  #(k+1, Q)matrix for itemJ
   PJk <- matrix(Pmat[(fromP[J]+1):toP[J], ], ncol=Q)  #(k,Q)matrix for itemJ
   Wj  <- matrix(c(0, wj), nrow=1, ncol=k+1)  #k+1 vector form itemJ
   oneQ <- matrix(1, nrow=Q, ncol=1)
   onek <- matrix(1, nrow=k, ncol=1)
   term2 <- (PJk * (wj %*% t(oneQ) - onek %*% Wj %*% PJ)) %*%
    diag(frac1^2) %*% solve(diag(frac2^2))

   d1 <- (term1 - term2) %*% th.dist
   t(d1)
  })
 }



 #-----------------------------------------
 # Eq(73)
 EIjD2 <- function(param, Pmat, dPmat, J, ws, Q, th.dist){
  return(function(wj){

   ncat <- param$ncat
   toP <- cumsum(ncat)
   fromP <- toP - ncat + 1

   k <- ncat[J] - 1  # length of  wj
   frac1 <- as.vector(SumdP(param, dPmat, J, ws, wj))  #(1, Q)
   frac2 <- as.vector(VarX(param, Pmat, J, ws, wj))    #(1, Q)
   Wj  <- matrix(c(0, wj), nrow=1, ncol=k+1)  #k+1 vector form itemJ
   PJ   <-  Pmat[(fromP[J]  ):toP[J], ]  #(k,Q)matrix for itemJ
   PJk <- matrix(Pmat[(fromP[J]+1):toP[J], ], ncol=Q)  #(k,Q)matrix for itemJ
   dPJk <- matrix(dPmat[(fromP[J]+1):toP[J], ], ncol=Q) #(k,Q)matrix for itemJ
   oneQ <- matrix(1, nrow=Q, ncol=1)
   onek <- matrix(1, nrow=k, ncol=1)

   #hess4
   hm4 <- matrix(0, nrow=k, ncol=k)
   for (q in 1:k){
    hm4[, q] <- dPJk %*% diag(as.vector(dPJk[q, ])) %*%
     solve(diag(frac2)) %*% th.dist
   }

   #hess3
   hm3 <- matrix(0, nrow=k, ncol=k)
   term3L1 <- wj %*% t(oneQ) - onek %*% Wj %*% PJ
   #[]
   for (q in 1:k){
    term3L <- PJk %*% diag(as.vector(dPJk[q, ])) * term3L1
    term3R1 <- wj[q] * onek %*% t(oneQ) - onek %*% Wj %*% PJ
    term3R <- dPJk %*% diag(as.vector(PJk[q, ])) * term3R1
    hm3[, q] <- ((term3L + term3R) %*% diag(frac1) %*%
                  solve(diag(frac2^2))) %*% th.dist
   }

   #hess2
   hm2 <- matrix(0, nrow=k, ncol=k)
   for (q in 1:k){
    term3R1 <- wj[q] * onek %*% t(oneQ) - onek %*% Wj %*% PJ
    term2 <- PJk %*% diag(as.vector(PJk[q, ])) * term3L1 * term3R1
    hm2[, q] <- term2 %*% diag(frac1^2) %*% solve(diag(frac2^3)) %*%
     th.dist     }

   #hess1
   hm1 <- matrix(0, nrow=k, ncol=k)
   for (q in 1:k){
    term1 <- PJk %*% diag(as.vector(PJk[q,]))
    term1[q,] <- term1[q,] - as.vector(PJk[q,])  # q
    hm1[,q] <- term1 %*% diag(frac1^2) %*% solve(diag(frac2^2)) %*%
     th.dist
   }

   hm <- matrix(0, nrow=k, ncol=k)
   m2 <- matrix(2, nrow=k, ncol=k)
   m4 <- matrix(4, nrow=k, ncol=k)

   hm <- hm1 + m4*hm2 - m2*hm3 + hm4
   hm
  })
 }



