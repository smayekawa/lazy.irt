#-------------------------------------------------------------
# (1)Esimation of Globally Optimal Weights with  maxNR function
#-------------------------------------------------------------
# L: maximum number of iterations
# [thmin, thmax]の範囲にnpoints個の区分点を作成する。
# 期待情報量を計算するときに，N(th.m, th.sd^2)の正規分布を用いる。
# i.weight：項目全体に掛かる重み(w_j)
# c.weight：カテゴリの重み(v_kj) 


GOptWeight <- function(param, i.weight=NULL, c.weight=NULL, npoints=31, thmin=-4, thmax=4, th.m=0, th.sd=1, L=100, print=0){  

  nitems <- nrow(param)
  iname <- param$name
  ncat <- param$ncat
  Q <- npoints
  
# irf & dirf
  temp1 <-  irf(param, weight=NULL, npoints=Q, thmin=thmin, thmax=thmax, print=print)
  temp2 <- dirf(param, weight=NULL, npoints=Q, thmin=thmin, thmax=thmax, print=print)
  Pmat  <- t(temp1$ICRF)    #sum(ncat)行Q列の行列
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

#theta
  theta <- seq(thmin, thmax, length.out=npoints)
  th.dist <- dnorm(theta, m=th.m, sd=th.sd)/ sum(dnorm(theta, m=th.m, sd=th.sd))
  
  m <- 0; mdif <- 100
  for (l in 1:L){
    if (mdif < 1e-06) break
  	 cat("--------------------------------- \n")
    for (J in 1:nitems){
      wsj <- ws[(fromP[J]+1):toP[J]]
      ans <- maxNR(fn = EIj(param, Pmat, dPmat, J, ws, th.dist),
                   grad = EIjD1(param,Pmat, dPmat, J, ws, Q, th.dist),
                   hess = EIjD2(param,Pmat, dPmat, J, ws, Q, th.dist), 
                   start = wsj, print.level=0, tol = 1e-01, iterlim = L)
      ws[(fromP[J]+1):toP[J]] <- ans$estimate
      cat("item=",iname[J]," maximum=",ans$maximum," code=",ans$code, "\n")
    }
  	 mdif <- abs((ans$maximum-m)/m)  #変化した割合でチェックする
  	 m <- ans$maximum
  	 cat("--------------------------------- \n")
  	 cat("Iteration=",l,",   maximum=",ans$maximum,",   mdif=",mdif, "\n")
  	 Print(ws)
  }
  return(ws)
}



#-----------------------------------------
#シンプルな重みを作り，ベクトル形式で返す関数。
#例）ncat=c(2,3); nitems=2 ならば，c(0,1,0,1,2)が返る。
#-----------------------------------------
make.simple.weight <- function(ncat, nitems){
  ws <- numeric(0)
  for (j in 1:nitems){
    ws <- append(ws, c(0:(ncat[j]-1)))
  }
  return(ws)
}



#-------------------------------------------------------------
# weights (vector -> data.frame)
#-------------------------------------------------------------
weight.df <- function(param, ws){
  nitems <- nrow(param)
  iname <- param$name
  ncat <- param$ncat
  type <- param$type
  
  toP <- cumsum(ncat)
  fromP <- toP - ncat + 1
        
  w <- matrix(1, nitems)
  rownames(w) <- iname
  colnames(w) <- "w"

  v <- matrix(NA, nitems, max(ncat))
  rownames(v) = iname
  colnames(v) = paste("v", 0:(max(ncat) - 1), sep = "")
  for (j in 1:nitems) {
    v[j, 1:ncat[j]] <- ws[fromP[j]:toP[j]]
  }
  
  weight <- data.frame(iname, type, ncat, w, v)
  rownames(weight) <- iname
  colnames(weight) <- c("name", "type", "ncat", "w", colnames(v))
  
  return(weight)
}


#-----------------------------------------
### Small parts of functions           ###
#-----------------------------------------
# 式66の分子の二乗の中身。
# 1×Qのベクトルを返す
# wは長さsum(ncat)のカテゴリ重みベクトル。
# 項目Jの重みに相当する部分：wj（ベクトル）だけが変数。
# dPmatはsum(ncat)行Q列の行列
SumdP <- function(param, dPmat, J, ws, wj){
  ncat <- param$ncat
  toP <- cumsum(ncat)
  fromP <- toP - ncat + 1
    
  w <- ws
  w[fromP[J]:toP[J]] <- c(0, wj[1:(ncat[J]-1)])
  return(w %*% dPmat)
}

##式66の分母
# wは長さsum(ncat)のカテゴリ重みベクトル。
# 項目Jの重みに相当する部分：wj（ベクトル）だけが変数。
# Pmatはsum(ncat)行Q列の行列
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
# カテゴリの重みが未知の場合
#==========================================
#-----------------------------------------
# Expected Information (function of wj) 式(66)
EIj <- function(param, Pmat, dPmat, J, ws, th.dist){
  return(function(wj){
    ((SumdP(param, dPmat, J, ws, wj)^2) / VarX(param, Pmat, J, ws, wj)) %*% th.dist
  })
}



#-----------------------------------------
# Expected Information (function of wj) with d1 式(68)
EIjD1 <- function(param, Pmat, dPmat, J, ws, Q, th.dist){
  return(function(wj){

  ncat <- param$ncat
  toP <- cumsum(ncat)
  fromP <- toP - ncat + 1

    k <- ncat[J] - 1  #着目している項目Jのカテゴリ数-1。探したいwjの長さ。
    d1 <- numeric(k)

    frac1 <- as.vector(SumdP(param, dPmat, J, ws, wj))  #(1, Q)
    frac2 <- as.vector(VarX(param, Pmat, J, ws, wj))    #(1, Q)
  #term1
    dPJk <- matrix(dPmat[(fromP[J]+1):toP[J], ], ncol=Q)  #(k, Q)matrix for itemJ
    term1 <- dPJk %*% diag(frac1) %*% solve(diag(frac2))
  #term2
    PJ  <- Pmat[(fromP[J]  ):toP[J], ]  #(k+1, Q)matrix for itemJ
    PJk <- matrix(Pmat[(fromP[J]+1):toP[J], ], ncol=Q)  #(k,   Q)matrix for itemJ
    Wj  <- matrix(c(0, wj), nrow=1, ncol=k+1)  #k+1 vector form itemJ
    oneQ <- matrix(1, nrow=Q, ncol=1)
    onek <- matrix(1, nrow=k, ncol=1)
    term2 <- (PJk * (wj %*% t(oneQ) - onek %*% Wj %*% PJ)) %*% diag(frac1^2) %*% solve(diag(frac2^2))

    d1 <- (term1 - term2) %*% th.dist
    t(d1)
  })
}

 

#-----------------------------------------
#式(73)
EIjD2 <- function(param, Pmat, dPmat, J, ws, Q, th.dist){
  return(function(wj){
 
  ncat <- param$ncat
  toP <- cumsum(ncat)
  fromP <- toP - ncat + 1

    k <- ncat[J] - 1  #着目している項目Jのカテゴリ数-1。探したいwjの長さ。
    frac1 <- as.vector(SumdP(param, dPmat, J, ws, wj))  #(1, Q)
    frac2 <- as.vector(VarX(param, Pmat, J, ws, wj))    #(1, Q)
    Wj  <- matrix(c(0, wj), nrow=1, ncol=k+1)  #k+1 vector form itemJ
    PJ   <-  Pmat[(fromP[J]  ):toP[J], ]  #(k, Q)matrix for itemJ
    PJk <- matrix(Pmat[(fromP[J]+1):toP[J], ], ncol=Q)  #(k,   Q)matrix for itemJ
    dPJk <- matrix(dPmat[(fromP[J]+1):toP[J], ], ncol=Q)  #(k, Q)matrix for itemJ
    oneQ <- matrix(1, nrow=Q, ncol=1)
    onek <- matrix(1, nrow=k, ncol=1)

  #hess4
    hm4 <- matrix(0, nrow=k, ncol=k)
    for (q in 1:k){
      hm4[, q] <- dPJk %*% diag(as.vector(dPJk[q, ])) %*% solve(diag(frac2)) %*% th.dist
    }

  #hess3
    hm3 <- matrix(0, nrow=k, ncol=k)
    term3L1 <- wj %*% t(oneQ) - onek %*% Wj %*% PJ
    #[]
    for (q in 1:k){
      term3L <- PJk %*% diag(as.vector(dPJk[q, ])) * term3L1
      term3R1 <- wj[q] * onek %*% t(oneQ) - onek %*% Wj %*% PJ 
      term3R <- dPJk %*% diag(as.vector(PJk[q, ])) * term3R1
      hm3[, q] <- ((term3L + term3R) %*% diag(frac1) %*% solve(diag(frac2^2))) %*% th.dist
    }

  #hess2
    hm2 <- matrix(0, nrow=k, ncol=k)   #空行列を作る
    for (q in 1:k){
      term3R1 <- wj[q] * onek %*% t(oneQ) - onek %*% Wj %*% PJ 
      term2 <- PJk %*% diag(as.vector(PJk[q, ])) * term3L1 * term3R1
      hm2[, q] <- term2 %*% diag(frac1^2) %*% solve(diag(frac2^3)) %*% th.dist     }

  #hess1
    hm1 <- matrix(0, nrow=k, ncol=k)   #空行列を作る
    for (q in 1:k){
      term1 <- PJk %*% diag(as.vector(PJk[q,]))
      term1[q,] <- term1[q,] - as.vector(PJk[q,]) #q番目の項の時のみ，PqJを引く
      hm1[,q] <- term1 %*% diag(frac1^2) %*% solve(diag(frac2^2)) %*% th.dist
    }

    hm <- matrix(0, nrow=k, ncol=k)   #空行列を作る
    m2 <- matrix(2, nrow=k, ncol=k)   # 2の行列
    m4 <- matrix(4, nrow=k, ncol=k)   # 4の行列

    hm <- hm1 + m4*hm2 - m2*hm3 + hm4
    hm
  })
}

