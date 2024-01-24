uIRT <- function( Uc, groupvar=NULL, idvar=NULL, ncat=NULL, type=NULL, DinP=1
                , param=NULL, msn=NULL, fixeditems=NULL, baseform=1
                , theta=NULL, thd=NULL, estmu=0, estsigma=0
                , npoints=21, thmin=-4, thmax=4
                , maxiter=200, eps=1e-7, epsd=1e-6
                , maxiter2=20, eps2=1e-4, nstrict=9, maxiter22=5
                , maxabsparam=20
                , SQUAREM=3, minalpha=-6, nSQUAREM=1
                , print=2, plot=0, smallP=0, debug=0
                ){
 # unidimensional IRT parameter estimation
 # Shin-ichi Mayekawa
 # 121024,25
 # multi-group: 121026
 # group and id as a column of Uc: 121028,29
 # baseform: 121031
 # print: 121031,1101
 # output lmlh: 121103
 # DinP: 121110(London)
 # epsd, location of conv. check, plot: 121110(London)
 # initial msn: 121112
 # calculate H only in getlmlh: 121112
 # smallP to ordinal_reg: 121115
 # limits to args of exp or log: 121117
 # initial p3=c=0.3 for type="B3" items: 121118
 # group names etc: 121119
 # nominal response model: 121123,24
 # fixed parameters: 121124,25,26
 # inline SQUAREM: 121126,27
 # epsd for maxadp: 121128
 # positive minalpha: 121129
 # nstrict, nSQUAREM: 121130
 # check if the input dataframe exists: 121215
 # estimate baseform theta dist if fixeditems: 121216
 # avoid 3d array and exclude fixed items from SQUAREM: 121216
 # include ms in SQUAREM: 121217
 # calculation of R in getlmlh: 121217
 # new print: 121230
 # maxabsparam for ordinal_reg: 130204
 #
 #
 # Args:
 # 
 #  Uc       n x nitems+1 or +2    compressed item resopnse data frame
 #           in BILOG-MG's expanded format
 #  groupvar name of the grouping varible contained as a column of Uc
 #  idvar    name of the id varible contained as a column of Uc
 #
 #  ncat     nitems x 1  # of categories for each item
 #  type     nitems x 1  item type as   "Bn" | "G" | "PN" | "P"
 #
 #  param    nitems x max(ncat)  initial value data frame for item param
 #           Only the subset of the items can be given.
 #           param$name and the colname(Uc) will be used to match the items.
 #  msn      nG x 3  initial value matrix of (mean, std, n) for each group
 #
 #  fixeditems list of items whose item parameters are to be fixed 
 #           to the values given in the param data frame.
 #           item numbers or item names
 #
 #  baseform base form number whose theta distribution is fixed at the values
 #           given in msn[baseform,]
 #
 #           If there are fixed items, the values of the item paramters
 #           given in the param data frame must be on the baseform scale.
 #
 #  estmu    =1  to estimate multigroup theta means
 #  estsigma =1  to estimate multigroup theta std
 #
 #  theta    npoints x 1   discrete theta points
 #  npoints, # of discrete theta points between (thmin, thmax)
 #
 #  maxiter  max # of iterations
 #  eps      eps for the relative improvement of lmlh
 #  epsd     eps for the max. abs. diff. of msn
 #           Seems this is important.
 #  maxiter2 max # of iterations for ordinal_reg and smn when llll <= nstrict
 #  maxiter22  max # of iterations for ordinal_reg and smn when llll > nstrict
 #  eps2     eps for the relative improvement of llh in ordinal_reg and smn
 #  nstrict  maxiter2 will be reduced to maxiter22 after nstrict iterations
 #
 #  SQUAREM  = 1 to emply SQUAREM to accelerate
 #  minalpha minimum value of alpha parameter:  default is -6.
 #           If minalpha > 0,  alpha will be fixed to minalpha.
 #  nSQUAREM iteration number from which SQUAREM update begins
 #
 #
 #
 #  Values:
 #   list of
 #     param   estimated parameter data frame
 #     msn     mu, sigma, n  of theta distriution for each group
 #     (theta, thd) theta distribution
 #     H       n x npoints   posterior of theta given data
 #     (EAP, poststd)   n x 2    estimated theta and its std
 #   
 #
 #  Needs:
 #   dummy_expand, irf, dirf_p, ordinal_reg, smn
 #
 #
 #
 #  Comments:
 #   When theta distribution is estimated, mu and sigma for the baseform group
 #   is fixed to (0,1).
 #
 #   When some of the item parameters are fixed, 
 #   theta distribution of the baseform group will also be estimated 
 #   if estmu == 1.
 #
 #
 #   Method SqS3 is implemented for SQUAREM.
 #     alpha > minalpha is enforced.
 #
 #
 #
 #
 #
 #
 #   If SQUAREM == 3 and minalpha == -6,  it seems maxiter2=1 seems to be OK.
 #   However, for the first few iterations, use maxiter2=10 or so.
 #   maxiter2=1 may cause problem when n or nitems are large.
 #
 #
 
 
 
 
 getlmlh <- function( U, group, param, theta, thd ){
  # calculation of log marginal likelihood and the posterior expectation 
  # Shin-ichi Mayekawa
  # 20121024
  # output H: 121112
  # limits of exp or log: 121117,18,19
  # return R: 121217
  #
  # Args:
  #   U      n x sum(ncat)  expanded data
  #          U[i,fromP[j]:toP[j]]=0 if Uc[i,j] is missing
  #   group  n x 1
  #   param  parameter data frame: name, type, ncat, p1, p2, ...
  #   theta  npoint x 1   discrete theta points
  #   thd    npoint x nG  prob at theta points for each group (sums to unity)
  #
  # Values:
  #   lmlh   log marginal likelihood
  #   logP   npoints x  log of IRF
  #   H      n x npoints  posterior expectation of theta given obs
  #   R      npoints x sum(nitems)
  #
  #  Note
  #     min arg to 1/1e-x is 308
  #     min arg to log is 1e-308:
  #     max arg to exp is 709.782: 
  #
  
  nG=length(unique(group))
  P=irf( param, theta, DinP=DinP, print=0 )$ICRF
P[P<10e-324]=10e-324
  logP=log( P ); rm(P)
  H=matrix(0,nrow(U),length(theta))
  lmlh=0
  R=array(NA,c(ncol(H),ncol(U),nG))
  for( g in 1:nG ){
   locg=which(group==g)
   # eULP=exp(U[locg,,drop=0]%*%t(logP)) * matrix(1,length(locg))%*%t(thd[,g])
   ULP=U[locg,,drop=0]%*%t(logP)
   ULP[ULP > 709]=709
   eULP=exp(ULP)*matrix(1,length(locg))%*%t(thd[,g])
   rseULP=rowSums(eULP)
   rseULP[rseULP < 1e-307]=1e-307
   lmlh=lmlh+sum( log( rseULP ) )
   H[locg,]=eULP/rseULP   
   R[,,g]=t(H[locg,])%*%U[locg,,drop=0]
  }
  rm(ULP,eULP,rseULP)
  
  return( list(lmlh=lmlh, logP=logP, H=H, R=R) )
  
 } # end of getlmlh
 
 
 
 
 # chech if dataset exists
 dfname=deparse(substitute(Uc))
 dfname=unlist( strsplit(dfname,"[",fixed=1) )[1]
 if( !exists(dfname) ){
  cat("\n\nerror1(uIRT) **** data frame ", dfname, " does not exist. ****\n\n")
  return(NULL)
 }
 else if( !is.data.frame(Uc) ){
  cat("\n\nerror1(uIRT) **** ", dfname, " is not a data frame. ****\n\n")
  return(NULL)
}

 
 # grouping variable
 if( is.null(groupvar) ){
  groupvar="none"
  group=rep(1,nrow(Uc))
  baseform=1
 }
 else{
  locgvar=which(colnames(Uc)==groupvar)
  if( length(locgvar) == 0 ){
   cat("\n\n error1(uIRT): group variable does not exist.\n\n\n")
   Print(groupvar)
   return( NULL )
  }
  # make group an integer variable containing 1,2,...,nG
  group=Uc[,locgvar]
  Uc=Uc[,-locgvar,drop=0]
  ugroup=unique(group)
  nG=length(unique(group))
  names(ugroup)=paste("group",1:nG,sep="")
  group0=group
  for( g in 1:nG ) group0[group0 == ugroup[g]]=g
  group=group0
  rm(group0)
 }
 nG=length(unique(group))
 if( nG == 1 ) baseform=1
 
 # id variable
 if( is.null(idvar) ){
  idvar="none"
  id=paste("s",1:nrow(Uc),sep="")
 }
 else{
  locidvar=which(colnames(Uc)==idvar)
  if( length(locidvar) == 0 ){
   cat("\n\n error1(uIRT): ID variable does not exist.\n\n\n")
   Print(idvar)
   return( NULL )
  }
  #
  id=Uc[,locidvar]
  Uc=Uc[,-locidvar,drop=0]
 }
 
 # const
 n=nrow(Uc)
 nitems=ncol(Uc)
 itemname=colnames(Uc)
 nobs=table(group)

 # missing data count
 pmiss=sum(is.na(Uc))/(n*nitems)
 
 
 # dummy expand the compressed data
 #   U[i,fromP[j]:toP[j]]=0 if Uc[i,j] is missing
 temp=dummy_expand( Uc )
 rm(Uc)
 U=temp$U
 fromP=temp$fromP; toP=temp$toP; ncat=temp$ncat
 rm(temp)

 # item type
 if( is.null(type) ){
  type=rep("B",nitems)
  type[which(ncat>2)]="P"
 }
 else{
  type=type[1:nitems]
  if( !( (length(grep("^B",type[ncat == 2])) == sum(ncat == 2))
         &&  (length(grep("^B",type[ncat == 3])) == 0)   ) ){
   cat("\n\n error1(uIRT): Item type does not match the # of categories.\n\n")
   Print(param)
   return(NULL)
  }
 } 
 
 # of parameters of the model:  
 np=ncat
 np[which(type == "B3")]=3 
 np[which(type == "N")]=2*(ncat[which(type == "N")]-1)
 maxnp=max(np)
 if( length(grep("^B[2]*$",type)) == nitems ) maxnp=maxnp+1

 if( debug > 0 ) 
  Print(  nitems, nG, npoints, "/", type, ncat, np )
 
 # initial parameter data frame
 paramORG=param
 # create an empty data frame
 paramhead=as.data.frame( cbind( name=itemname , type=type )
                        , stringsAsFactors=0 ) 
 paramnum=matrix(NA,nitems,maxnp+1
             , dimnames=list(itemname,c("ncat", paste("p",1:maxnp,sep=""))))
 param=cbind(paramhead, paramnum)
 rownames(param)=paste("item",1:nitems,sep="")
 if( debug > 0 ) Print(" initial param with NA", param) 


 # insert the defined part of paramORG into param
 initemlist=NULL
 if( !is.null(paramORG) ){
  for( j in 1:nrow(paramORG) ){
   paramj=paramORG[j,]
   namej=paramj$name
   typej=paramj$type
   npj=paramj$ncat
   if( length(grep("^B[[:digit:]]*$",typej)) > 0 ) npj=3
   else if( typej == "N" ) npj=2*(npj-1)
   if( !any(is.na(paramj[3:(3+npj)])) ){
    locj=which( itemname == namej )
    if( length( locj ) > 0 ){
      param[locj,1:(3+npj)]=paramj[1:(3+npj)]
      initemlist=c(initemlist,locj)
    }
   }   
  }
 }
 rm(paramORG)
  if( debug > 0 ) 
  Print(  "after inserting paramORG", type, ncat, np, param )

 
 # default initial parameter values
 for( j in 1:nitems ){
  if( is.na(param$ncat[j]) ){
   # initial parameter values or fixed values are not given for item j.
   param$ncat[j]=ncat[j]; param$p1[j]=1
   if( type[j] %in% c("G", "PN", "P") ){
    # equally spaced b parameters
    param[j,5:(np[j]+3)]=(1:(np[j]-1))-np[j]*(np[j]-1)/2/(np[j]-1)
   }
   else if( type[j]  %in% c("B", "B2") ){
    param$p2[j]=0; param$p3[j]=0
   }
   else if( type[j] == "B3" ){
    param$p2[j]=0; param$p3[j]=0.3   # do not use 0 for 3PLM.
   }
   else if( type[j] == "N" ){
    #Print("****",j)
    param[j,4:(4+ncat[j]-1-1)]=1          # a paramters
    param[j,(4+ncat[j]-1):(4+np[j]-1)]=0  # b paramters
   }
  } # end of initial not given for item j  
 }
 
 # item paramters fixed or not
 if( length(fixeditems) > 0 && length(initemlist) > 0 ){
  if( is.character(fixeditems) )
        fixeditems=which( itemname %in% fixeditems )
  fixeditemlist=intersect(fixeditems,initemlist)
  if( length(fixeditemlist) > 0 ){
   FixedItems=1
  }
 } 
 else{
  FixedItems=0; fixeditemlist=NULL 
 }
 nonfixeditemlist=setdiff(1:nitems,fixeditemlist)
 
 # location of non-fixed params excluding the c-parameters of B2 items.
 pp=param[nonfixeditemlist,4:(4+maxnp-1)]
 ppt=param$type[nonfixeditemlist]
 pp$p3[which(ppt == "B" | ppt == "B2")]=NA
 locSQEM=which(!is.na(pp))
 # if( debug ) Print(pp,locSQEM)
 rm(pp,ppt)
 # param[nonfixeditemlist,4:(4+maxnp-1)][locSQEM]

 
 # c-parameters of B2 items : all items, not restrocted to nonfixeditems
 locB2=which(param$type == "B" | param$type == "B2")
 cB2=param$p3[locB2]
 
 # # of paramters to be estimated
 nparam1=sum(np)-sum(np[fixeditemlist])
 nparam2=0
 if( estmu == 1 ) nparam2=nparam2+nG-1
 if( estsigma == 1 ) nparam2=nparam2+nG-1
 if( FixedItems ){
  if( estmu == 1 ) nparam2=nparam2+1  # for the baseform
  if( estsigma == 1 ) nparam2=nparam2+1  # for the baseform
 } 
 nparam=nparam1+nparam2
 
 
 # generate theta
 if( is.null(theta) ){
  theta=matrix(seq(thmin,thmax,length.out=npoints),,1)
 }
 else{
  npoints=length(as.vector(theta))
 }
 
 # theta distribution
 thd=matrix(0,npoints,nG)
 rownames(thd)=format(theta,2)
 if( is.null(msn) ){
  # make all N(0,1)
  msn=matrix(c(0,1,1),nG,3,byrow=1)
 }
 for( g in 1:nG ){
  temp=exp(-0.5*( (theta-msn[g,1])/msn[g,2] )^2); 
  thd[,g]=temp/sum(temp)
 }
 msn[,3]=nobs
 colnames(msn)=c("mu","sigma","n")
 rownames(msn)=paste("group",1:nG,sep="")

 # what to estimate
 estms=matrix(0,nG,2)
 if( (estmu == 1  ||  estsigma == 1) )
  for( g in 1:nG )
   if( FixedItems == 1  ||  g != baseform )
    estms[g,]=c(estmu,estsigma)
 locms=which(estms == 1)

 
 if( print > 0 ){
  cat("\n\n\nuIRT: parameter estimation of unidimensional IRT models\n")
  cat("\n id varaible =", idvar, ",  group varaible =", groupvar,"\n")
  cat(" # of observations =", n, ",  # of items =", nitems,"\n")
  cat(" # of groups =", nG, "\n")
  if( nG > 1 ){
   cat("  ");  print(ugroup)
  }
  cat(" percentage of mising data =", pmiss,"\n")
  cat(" # of theta points =", npoints, " in [", thmin, ",", thmax,"]","\n")
  cat(" estimation of theta distribution\n")
  cat("   estmu =", estmu, ",  estsigma =", estsigma, "\n")
  cat("   baseform =", baseform, "\n")
  if( FixedItems ){
   cat(" The item paramters of the following items are fixed. \n")
   print(fixeditemlist)
   }
  cat(" total # of parameters to be estimated =", nparam
      ," (",nparam1," + ",nparam2 ,") \n")
  cat(" max abs value of estimated parameter =", maxabsparam,"\n")
  cat(" max # of iterations =", maxiter," with eps =",eps
      , " and epsd =", epsd,"\n") 
  cat(" max # of iterations in m-steps =", maxiter2,"/",maxiter22
      ," with eps =",eps2, "and nstrict =", nstrict,"\n") 
  if( SQUAREM){
   cat(" SQUAREM =",SQUAREM," will be used after", nSQUAREM, "iterations ")
   cat(" with min(alpha) =", minalpha,"\n")
  }
  cat(" print level =", print,"\n")
  cat(" D in P =", DinP,"\n")
  if( length(initemlist) > 0 ){
    cat(" initial item parameter values of the following items are given.\n")
    print(initemlist)
   }
  if( print >= 2 ){
   cat(" inital parameter value \n")
   Print(param)
   Print(msn)
   if( print >= 3 ){
    Print(thd,digits=1)
   }
  }
 }

 #initial
 temp=getlmlh( U, group, param, theta, thd )
 lmlhp=temp$lmlh
 logP=temp$logP
 H=temp$H
 R=temp$R
 rm(temp)

 if( print > 0 ) Print("initial", lmlhp)
 converged=0
 msnp=msn
 alpha=-1
 
 # define paramnum as a matrix consisting of non fixed items pnly
 paramnum=as.matrix(param[nonfixeditemlist,4:(4+maxnp-1)])
 paramnump=paramnum
 
 NN=matrix(0,npoints+1,nG); rownames(NN)=c(format(theta,2),"total")
 
 if( SQUAREM ){
  # SQUAREM
  paramnum_save=matrix(NA,length(locSQEM)+length(locms),3)
  # if( debug ) Print(nparam1,length(locSQEM),maxnp,ncol(param))
  paramnum_save[,1]=c(
                   as.matrix(param[nonfixeditemlist,4:(4+maxnp-1)])[locSQEM]
                 , msn[locms] )
  paramnum_save[,2]=paramnum_save[,1]
  paramnum_save[,3]=paramnum_save[,1]

  param0=param
  msn0=msn
  thd0=thd
  llll2=0
  nSQEM=0
  nSQEMs=0
 } 
 
 iter_hist=matrix(c(0,lmlhp,0,0,-1),1,5)
 colnames(iter_hist)=c("llll","lmlh","maxadp","maxadms","alpha")
 

 
 # main iteration
 for( llll in 1:maxiter ){
  
  
  # e-step
  # expected value of U at theta point
  # Rall is npoints x sum(ncat)
  Rall=apply(R,c(1,2),sum)
  
  if( llll <= nstrict ) maxiter21=maxiter2
  else maxiter21=maxiter22
 

  if( debug > 0 ) Print("*top of m-step for items")
  
  # m-step for item parameters
  for( j in 1:nitems ){
   if( !( j %in% fixeditemlist) ){
    Rj=Rall[,fromP[j]:toP[j]]
    paramj=param[j,]
    eps21=eps2
    printx=print-10
#Print(j, Rj, theta, type[j], paramj)
    temp=ordinal_reg( Rj, theta, type=type[j], param=paramj, DinP=DinP
                 , maxiter=maxiter21, eps=eps21, maxabsparam=maxabsparam
                 , print=printx, smallP=smallP )
    param[j,1:ncol(temp$param)]=temp$param
   }
  } # end of j loop
  rm(temp)
  
  if( debug > 0 ) Print("*top of m-step for theta dist")  
  
  # m-step for the theta distribution
  locg=which(group==baseform)
  # NN[1:npoints,baseform]=rowSums( t(H[locg,,drop=0])%*%U[locg,,drop=0] )
  NN[1:npoints,baseform]=rowSums( R[,,baseform] )
  if( (estmu == 1  ||  estsigma == 1) && llll >= 2 ){
   for( g in 1:nG ){
    if( FixedItems == 1 ||  g != baseform ){
     locg=which(group==g)
     Rg=R[,,g]
     Nj=rowSums(Rg)
     NN[1:npoints,g]=Nj
     temp=smn( theta, Nj, maxiter=maxiter2, eps=eps2, print=print-10
             , estmu=estmu, estsigma=estsigma
             , mu=msn[g,1], sigma=msn[g,2] )
     thd[,g]=temp$P
     msn[g,1]=temp$mu
     msn[g,2]=temp$sigma
     rm(temp)
    }
   }
   if( nG > 1) NN[npoints+1,]=colSums(NN[1:npoints,])
  } # end of theta dist
  
  if( debug > 0 ) Print("*top of convergence check")  
  # convergence:  moved to this place from above 121110
  temp=getlmlh( U, group, param, theta, thd )
  lmlh=temp$lmlh; H=temp$H; R=temp$R; logP=temp$logP; rm(temp)
  lmlhimpr=(lmlh-lmlhp)/abs(lmlhp)
  maxadp=max(abs(paramnump-param[nonfixeditemlist,4:(4+maxnp-1)]), na.rm=1)
  maxadms=max(abs(msnp-msn))
  if( print >= 2 ) Print( llll, lmlh, lmlhp, lmlhimpr,maxadp,maxadms
                         , fmt=c("i4",".6") )
  if( lmlhimpr < eps  && maxadp < epsd  && maxadms < epsd ){
   converged=1
   break
  }
  
  
  iter_hist=rbind(iter_hist,c(llll,lmlh,maxadp,maxadms,alpha))
  
  
  # inline SQUAREM
  if( llll >= nSQUAREM  &&  SQUAREM ){
  if( debug > 0 ) Print("*top of SQUREM")  
   # save history
   #
   # Better save the non-fixed part only.  not yet ready.
   # May be with mands.
   #
   llll2=llll2+1
   paramnum_save[,3]=paramnum_save[,2]   # 3 = theta_n (oldest)
   paramnum_save[,2]=paramnum_save[,1]   # 2 = F(theta_n)
   paramnum_save[,1]=c(
                    as.matrix(param[nonfixeditemlist,4:(4+maxnp-1)])[locSQEM]
                  , msn[locms] )
                     # 1 = F( F(theta_n) )

   # update every three iterationss
   if( llll2 %% 3 == 0 ){

    if( debug > 0 ) Print("*top of SQUREM mod 3 == 0")  
    
    # calculate alpha
    nSQEM=nSQEM+1
    r=paramnum_save[,2] - paramnum_save[,3]
    v=paramnum_save[,1] - 2*paramnum_save[,2] + paramnum_save[,3] 

    if( SQUAREM == 3 ) alpha=-sqrt(sum(r*r)/sum(v*v))
    else if( SQUAREM == 1 ) alpha=sum(v*r)/sum(v*v)
    if( alpha < minalpha ) alpha=minalpha
    if( alpha > -1 ) alpha=-1
    if( minalpha > 0 ) alpha=-minalpha    

    # update parameter
    newp=paramnum_save[,3]-2*alpha*r + (alpha^2)*v
    paramnum[locSQEM]=newp[1:nparam1]
    param0[nonfixeditemlist,4:(4+maxnp-1)]=paramnum
    if( length(locms) > 0 ){
     msn0[locms]=newp[(nparam1+1):nparam]
     thd0=matrix(0,npoints,nG)
     for( g in 1:nG ){
      temp=exp(-0.5*( (theta-msn[g,1])/msn[g,2] )^2); 
      thd0[,g]=temp/sum(temp)
     }
    }
    
    # evaluate lmlh
    temp=getlmlh( U, group, param0, theta, thd0 )
    lmlh0=temp$lmlh; H0=temp$H; R0=temp$R; logP0=temp$logP; rm(temp)

    # check if improved
    if( lmlh0 > lmlh ){
     # improved
     param=param0
     thd=thd0; msn=msn0
     lmlh=lmlh0; logP=logP0; H=H0; R=R0; 
     paramnum_save[,1]=c(
                as.matrix(param0[nonfixeditemlist,4:(4+maxnp-1)])[locSQEM]
              , msn[locms] )
     llll2=1   # counting this update as the first one.
     nSQEMs=nSQEMs+1
     maxadp=max(abs(paramnump-param[nonfixeditemlist,4:(4+maxnp-1)]), na.rm=1)
     iter_hist=rbind(iter_hist,c(llll,lmlh,maxadp,maxadms,alpha))
     if( print >= 2 ) {
      cat("SQUAREM success ")
      Print( lmlh0, lmlh, alpha,maxadp, fmt=".6" )
     }
    }
    else{
     # not improved
     llll2=2  # persistent setting
    }
    rm(logP0,H0,R0)
   } # end of mod=3
   
   
  } # end of SQUAREM

  
  # next
  lmlhp=lmlh   
  paramnump=param[nonfixeditemlist,4:(4+maxnp-1)]
  msnp=msn
  

  
 } # end of llll loop
 
 
 
 
 iter=llll
 if( print >= 1 && converged )
  cat("\nIteration converged with ", iter, " iterations.  eps= ",eps,"\n")

 if( print > 0 ){
  cat("\nuIRT: parameter estimation of unidimensional IRT models\n")
  cat("\n id varaible =", idvar, ",  group varaible =", groupvar,"\n")
  cat(" # of observations =", n, ",  # of items =", nitems,"\n")
  cat(" # of groups =", nG, "\n")
  if( nG > 1 ){
   cat("  ");  print(ugroup)
  }
  cat(" percentage of mising data =", pmiss,"\n")
  cat("\n # of theta points =", npoints, " in [", thmin, ",", thmax,"]","\n")
  cat(" estimation of theta distribution\n")
  cat("   estmu =", estmu, ",  estsigma =", estsigma, "\n")
  cat("   baseform =", baseform, "\n")
  if( FixedItems ){
    cat(" initial item parameter values of the following items are given.\n")
   print(fixeditemlist)
  }
  cat(" total # of parameters to be estimated =", nparam
      ," (",nparam1," + ",nparam2 ,") \n")
  cat(" max abs value of estimated parameter =", maxabsparam,"\n")
  cat("\n iteration terminated with", iter," iterations\n")
  cat(" max # of iterations in m-steps =", maxiter2,"/",maxiter22
      ," with eps =",eps2, "and nstrict =", nstrict,"\n") 
  if( SQUAREM){
   cat(" SQUAREM =", SQUAREM," will be used after", nSQUAREM, "iterations ")
   cat(" with min(alpha) =", minalpha,"\n")
   cat("    # of successes =", nSQEMs," out of", nSQEM, "\n")
  }
  cat("\n log marginal likelihood maximized =", lmlh,"\n\n") 
  cat(" max # of iterations =", maxiter," with eps =",eps
    , " and epsd =", epsd,"\n") 
  cat(" max # of iterations in m-steps =", maxiter2," with eps =",eps2
      , "and nstrict =", nstrict,"\n") 
  cat(" print level =", print,"\n")
  cat("\n estimated item parameters (DinP=", DinP,") \n")
  Printdf( param, fmt="s5 s5 i3 7.3" )
  cat("\n estimated mean and std of theta distribution \n")
  Print(msn)
  if( print >= 3 ){
   cat(" estimated # of persons at each of theta points\n")
   Print(NN)
  }
 }
  
 if( plot >= 1 ){
  maxthd=max(thd)
  title="theta distribution for each group"
  sub=NULL
  for( g in 1:(nG-1) ){
   plot(theta,thd[,g]
        , xlim=c(min(theta),max(theta)), ylim=c(0,maxthd), type="l", ylab="")
   par(new=1)
  }
  plot(theta,thd[,nG], main=title, sub=sub
       , xlim=c(min(theta),max(theta)), ylim=c(0,maxthd), type="l", ylab="prb")
  par(new=0)
 }
 
 EAP=H%*%theta
 poststd=sqrt(H%*%theta^2-EAP^2)
 return( 
  list( param=param, msn=msn, theta=theta, thd=thd
      , converged=converged, lmlh=lmlh,  iter_hist=iter_hist
      , H=H, NN=NN, EAP=EAP, pstdtd=poststd, id=id ) 
        )

 
} # end of uIRT



ordinal_reg <- function( R, x, type="G", param=NULL
                       , maxiter=100, eps=1e-8, epsd=1e-5, smallP=0, DinP=1
                       , minabsparam=1e-3, maxabsparam=20
                       , print=1 ){
 # ordinal regression of (R,N) on scalar x
 # Shin-ichi Mayekawa
 # 20121022,23,24
 # valid parameter value: 121025
	# bugfix: 121104
	# pass Pj to dirf_p: 121105
	# DinP: 121110(London)
 # check too large values of abs(paramnum): 121115
 # bigfix: 121118
 # dc=logit(c): 121118
 # initial p3=c=0.3 for type="B3" items: 121118
 # nominal model: 121123
 # matSwp: 130204 ????
 # maxabsparam: 130204
 #
 #  Args:
 # 
 #    R       npoints x ncat   response data matrix 
 #    x       npoints x 1      regressor matrix
 #    type    = B | B2 | B3 | G | P | PN | N   item type
 #    param  initial parameter data.frame: name, type, ncat, p1, p2, ...
 #
 #
 #  Values:    
 #   list of param, P etc
 #
 # Needs:
 #  irf, dirf_p
 #
 # This program finds the parameter which maximizes the following likelihood
 #    llh = sum( R*log(P(x | param)) )
 #  where R is the npoints x ncat dummy expanded responce matrix,
 #  x is the npoints x 1 regressor vector, and
 #  P(x | param) is the npoints x ncat probability vector calculated as
 #  the IRT binary or polytomous item category response function at theta=x.
 #  The Fisher Scoring method is used to estimate the parameters.
 # 
 
 # Print("**** Top of Ordinal Reg") 
 
 
 # const
 if( is.vector(x) ){
  nobs=length(x); nq=1
 }
 else{
  nobs=nrow(x); nq=ncol(x)
 }
 ncat=ncol(R)  # This is the # of categories.  Tha max score is m=ncat-1. 
 ncat1=ncat-1
 N=matrix(rowSums(R),,1)
 ntot=sum(N); nmiss=sum(N==0)
 nitems=1
 
  
 # error
 if( nq > 1 ){
  cat("\n\nerror:(ordinal_reg) *** Only one regressor is allowed. ******\n\n")
  return( NULL )
 }
 if( !substr(type,1,1) %in% c("B","G","P","N") ){
  cat("\n\nerror:(ordinal_reg) *** Incorrect type spec. **********\n\n")
  Print(type)
  return( NULL )
 }
 if( length(grep("^B[[:digit:]]*$",type)) > 0  && ncat >= 3 ){
  cat("\n\nerror:(ordinal_reg) *** Incorrect type spec. **********\n\n")
  Print(type, ncat)
  return( NULL )
 }
 
 # of parameters of the model:  
 if( type == "B3" ){
  np=3; ncatx1=3
 }
 else if( type == "N" ){
  np=2*ncat1; ncatx1=np
 }
 else{
  np=ncat; ncatx1=ncat
 } 

 # initial values
 if( is.null(param) ){
  # not given
  paramhead=as.data.frame( 
   cbind( name="Q1" , type=type )
   , stringsAsFactors=0 ) 
  paramnum=matrix(0,1,np+1
                  , dimnames=list(NULL,c("ncat", paste("p",1:ncatx1,sep=""))))
  paramnum[1,1]=ncat
  paramnum[1,2]=1
  # increasing b-parameters with mean=0
  if( length(grep("^B[[:digit:]]*$",type)) == 0 ){
   paramnum[,3:(np+1)]=(1:(np-1))-np*(np-1)/2/(np-1)
  }
  else if( type == "B3" ){
   paramnum[,3]=0
   paramnum[,4]=0.3
  }
  else{
   paramnum=cbind(paramnum,0)
   colnames(paramnum)[4]="p3"
  }  
  param=cbind(paramhead, paramnum)
 }
 else{
  # initial given
  # delete exessivecolumns
  locnotna=which( !apply( param, 2, function(x) all(is.na(x)) ) )
  param=param[,locnotna]
  param$type=type
 }
 
 if( is.null(param)  &&  type == "N" ){
  paramhead=as.data.frame( 
   cbind( name="Q1" , type=type )
   , stringsAsFactors=0 ) 
  paramnum=matrix(0,1,np+1
                  , dimnames=list(NULL,c("ncat", paste("p",1:ncatx1,sep=""))))
  paramnum[1,1:(np+1)]=c(ncat, 1:ncat1, (1:(ncat-1))-ncat*(ncat-1)/2/(ncat-1))
  paramnum[1,1:(np+1)]=c(ncat, rep(1,ncat1), rep(0,ncat1))
  param=cbind(paramhead, paramnum)
  Print(param)
 }
 
 if( print > 0 ){
  cat("\nOrdinal Regression\n")
  cat(" # of observations =", nobs, ",  # of regressors =", nq,"\n")
  cat(" total # of responses =", ntot, " with missing =", nmiss,"\n")
  cat(" Model Type =", type, " with the # of categories =",ncat,"\n")
  cat(" # of parameters to be estimated =", np, "\n")
  cat(" max # of iterations =", maxiter," with eps =",eps
      ," and epsd =", epsd,"\n") 
  cat(" maximum value of the absolute value of the parameters ="
     , maxabsparam, "\n")
  cat(" print level =", print,"\n")
  cat(" inital parameter value \n")
  Print(param)
  if( print >= 3 ){
   Print(x,R,N)
  }
 }
 

 # initial llh
 # P is npoints x ncat irf
 P=irf( param, x, smallP=smallP, DinP=DinP, print=0 )$ICRF
P[P<10e-324]=10e-324
 llh=sum(R*log(P))
 #Print("INITIAL", llh, R,N,P)
 
 # main iteration loop
 llhp=llh
 paramnump=param[,4:(4+np-1)]
 converged=0
 P1=NULL
 

 for( llll in 1:maxiter ){
  
  # log-Jacobian matrix by expanding the category first
  Jac=dirf_p( param, x, smallP=smallP, DinP=DinP, Pj=P1
              , log=1, cat.first=1, zero=1, print=0 )$Jack
  
  # simple d llh / d param
  #  gx=t( t(matrix(t(R),,1))%*%Jac )
  
  paramnum=param[,4:(4+np-1)]
  
  # expectation of  t(R[i,])%*%R[i,]
  Dr=NULL
  for( q in 1:nobs ){
   Drk=N[q]*( Diag(P[q,])-t(P[q,,drop=0])%*%P[q,,drop=0] )
   Dr=b_diag( Dr, Drk )
  }
  
  # type B3 or type PO or others
  if( type == "B3" ){
   
   # logit transformation of c=p3:    d logit(c) / d x = c*(1-c)
   locc=which(colnames(param)=="p3")
   c=param[1,locc]
   if( c == 0 ) c=1e-323
   Jac[,3]=Jac[,3]*c*(1-c)
   if( c > 0 ) dc=log(c/(1-c))
   else dc=-709
   paramnum[3]=dc
   
   # gradient etc
   delta=matrix(t((R-as.vector(N)*P)),,1)
   g=t(Jac)%*%delta
   H=t(Jac)%*%Dr%*%Jac
   
   # use sweep 
   td=t( solve(H)%*%g )
     
   paramnew=param
   llhnew=llhp
   
   step=1; ok=0; llhnew=llhp
   for( lllll in 1:20 ){
    paramnumnew=paramnum + step*td
   if( paramnumnew[1] <=  0.1 ) paramnumnew[1]=0.1  
#    if( paramnumnew[3] < -709 ) paramnumnew[3]=-709
    if( paramnumnew[3] < -20 ) paramnumnew[3]=-20
    paramnumnew[3]=1/(1+exp(-paramnumnew[3]))    
    paramnew[,4:(4+np-1)]=paramnumnew
    P=irf( paramnew, x, smallP=smallP, print=0 )$ICRF
P[P<10e-324]=10e-324
    llhnew=sum( R*log(P) )
    if( llhnew >= llhp ){
     ok=1
     break
    }    
    step=step/2
   } # end of lllll   
   
  } # end of B3
  else if( type == "PO" ){
   
   comments('***')
   
  } # end of PO
  else{
   
   # Nominal, 2PLM, Graded, or GPCM w/o restrictions
   
   # gradient etc
   delta=matrix(t((R-as.vector(N)*P)),,1)
   g=t(Jac)%*%delta
   H=t(Jac)%*%Dr%*%Jac
   
   td=t( matSwp(H)%*%g )
   
   paramnew=param
   llhnew=llhp
   
   step=1; ok=0; llhnew=llhp
   for( lllll in 1:20 ){
    paramnumnew=paramnum + step*td
    paramnumnew[paramnumnew > maxabsparam]=maxabsparam
    paramnumnew[paramnumnew < -maxabsparam]=-maxabsparam
    if( type != "N"  &&  paramnumnew[1] <=  0.1 ) paramnumnew[1]=0.1
    paramnew[,4:(4+np-1)]=paramnumnew
P=irf( paramnew, x, smallP=smallP, print=0 )$ICRF
P[P<10e-324]=10e-324
    llhnew=sum( R*log(P) )
    if( llhnew >= llhp ){
     ok=1
     break
    }    
    step=step/2
   } # end of lllll   
   
   
  } # end of others
  
  
  if( ok ){
   # update parameters
   param=paramnew
   llh=llhnew
  }
  else{
   if( print >= 0 )
    cat("\n\n(ordinal_reg)**** Iteration did not converge. ****\n\n")
   llhimpr=0
   maxag=9999
   break
  }
  
  # convergence
  maxag=max(abs(g))
  maxadp=max(abs(paramnump-paramnumnew))
  llhimpr=(llh-llhp)/abs(llhp)
  if( print >= 2 ) Print(llll,llh,llhp,llhimpr, maxadp, maxag, digits=3)
  if( llhimpr <= eps  &&  maxadp <= epsd ){
   converged=1;
   break
  }
  
  # next iteration
  llhp=llh
  paramnump=paramnumnew
  P1=P[,2:ncat,drop=F]
  
 } # end of llll loop 
 
 
 
 iter=llll
 if( print >= 1 && converged )
  cat("\nIteration converged with ", iter, " iterations.  eps= ",eps,"\n")
 
 if( print > 0 ){
  cat("\nOrdinal Regression\n")
  cat(" # of observations =", nobs, ",  # of regressors =", nq,"\n")
  cat(" total # of responses =", ntot, " with missing =", nmiss,"\n")
  cat(" Model Type =", type, "\n")
  cat(" # of parameters to be estimated =", np, "\n")
  cat(" # of iterations before termination =",iter)
  cat(" with eps =", eps, " and epsd =", epsd,"\n")
  cat(" final log likelihood value =", llh," with relative imp. ="
      , llhimpr," \n")
  cat(" max value of the gradient =", maxag, "\n")
  cat(" print level =", print,"\n")
  cat("\n estimated parameter value  (DinP =",DinP,")\n")
  print(param)
  cat("\n final gradient value\n")
  print(t(g))
 }
 
 return( list(param=param, llh=llh, P=P, x=x, g=g, R=R, maxag=maxag) )
 
} # end of ordinal_reg




