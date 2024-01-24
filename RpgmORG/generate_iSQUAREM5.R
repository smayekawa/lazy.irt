#' Inline Squared Extrapolation Methods for accelerating fixed-point iterations
#'
#' This function generates iSQUAREM function as a closure
#' in whose environment the past values of the parameters and related
#' information are stored.
#'
#' @param param A vector of parameters
#' @param SQUAREM = 3 or 1
#' @param debug = 1 to print the genereted function
#'
#' @details
#' See the help of iSQUAREM function.
#'
#' @references
#' R Varadhan and C Roland (2008),
#' Simple and globally convergent numerical schemes for acceler-
#' ating the convergence of any EM algorithm,
#' Scandinavian Journal of Statistics, 35:335-353.\cr
#' C  Roland,  R  Varadhan,  and  CE  Frangakis  (2007),
#' Squared  polynomial  extrapolation  methods
#' with cycling: an application to the positron emission tomography problem,
#' Numerical Algorithms 44:159-172.
#'
#' @export
#'

generate_iSQUAREM <- function( param, debug=0 ){
 # function to generate iSQUAREM  function
 # Shin-ichi Mayekawa
 # 20151219(W9-607)
 # 20151220,21,27,28dnc,29
 # badness0 etc renamed: 20160103
 # list from badness_of_fit: 20160103   bf NOT USED.
 #

 # variables in this environment to be accessed from the function
 # generated in this function.


 #############################################################(1)
 nparam=length(param)
 param_save=matrix(NA,nparam,3)
 param_save[,1]=param
 param_save[,2]=param
 param_save[,3]=param
 llll2=0
 nSQEM=0
 nSQEMs=0
 badness0=999999999999
 badness1=NA
 alpha=NA
 #############################################################(1)




 ###################################################################(2)
 iSQUAREM=function( param, ...
                    , enforce_constraints=NULL, badness_of_fit=NULL
                    , SQUAREM=3, minalpha=-999, maxalpha=-1, bof_value=NULL
                    , always=0, reset1=1, reset2=2
                    , print=0, debug=0 ){
  # inline SQUAREM
  # This function stores the current param vector in param_save
  # and returns the updated param vector in every 3 call.
  #
  # Do NOT copy and pasete this function:
  # This function must be created by generate_iSQUAREM as a closure.
  #

  if( debug ){
   Print("top of iSQUAREM",param)
  }

  # save history
  #
  llll2 <<- llll2+1
  param_save[,3] <<- param_save[,2]          # 3 = theta_n (oldest)
  param_save[,2] <<- param_save[,1]          # 2 = F(theta_n)
  param_save[,1] <<- param                   # 1 = F(F(theta_n))

  if( !is.null(bof_value) ) badness0 <<-bof_value

  # update every three iterationss
  if( llll2 %% 3 == 0 ){

   # calculate alpha
   nSQEM <<- nSQEM+1
   r=param_save[,2] - param_save[,3]
   v=param_save[,1] - 2*param_save[,2] + param_save[,3]
   if( SQUAREM == 3 ) alpha=-sqrt(sum(r*r)/sum(v*v))
   else alpha=sum(v*r)/sum(v*v)

   if( alpha < minalpha ) alpha=minalpha
   if( alpha > maxalpha ) alpha=maxalpha

   # update parameter
   newparam=param_save[,3] -2*alpha*r + (alpha^2)*v

   # constraints
   if(!is.null(enforce_constraints) )
    newparam=enforce_constraints( newparam, ... )

   # evaluate badness of fit
   if( !always ){
    bf=badness_of_fit( newparam, ... )
    if( is.list(bf) ) badness1=bf$critval
    else badness1=bf
    rm(bf)
   }

   # Print(badness0,  badness1)

   # check if improved
   # badness0 is the value of he badness_of_fit function
   # evaluated with the latest SQUAREM update in this function,
   # or it is supplied by the bof_value argument.

   if( always  ||  badness1 <= badness0 ){

    # improved
    param_save[,1] <<- newparam
    badness0 <<- badness1

    # llll2 <<- 0      # counting next update as the first one.
    llll2 <<- reset1   # counting this update as the first one.
    nSQEMs <<- nSQEMs+1

    if( print ){
     if( !always ) cat("SQUAREM success. ")
     else  cat("SQUAREM employed. ")
     Print( badness1, alpha, fmt=c(".7"))
    }

    res=list( param=newparam, critval=badness1 )
    return( res )

   } # end of success
   else{

    # not improved: do nothing
    # llll2 <<- 1  #
    llll2 <<- reset2  # persistent setting

   } # failure

   res=list( param=param, critval=badness1 )
   return( res )


  } # end of mod=3
  else{

   res=list( param=param, critval=badness1 )
   return( res )

  }

 } # end of iSQUAREM function

 if( debug ){
  func=deparse(iSQUAREM, width.cutoff=79, control="useSource")
  cat("\nThe following function was generated:\n")
  print(func,quote=0)
 }

 return( iSQUAREM )

} # end of generate_iSQUAREM







#' iSQUAREM
#' A function to be used to update parameters using inline SQUAREM: \cr
#' Inline Squared Extrapolation Methods for accelerating fixed-point iterations
#'
#' This function must be generated using generate_iSQUAREM function.
#'
#'
#' @param param A parameter vector to be updated
#' @param ... Additional parameters to the functions enforce_constraints
#' and badness_of_fit
#' @param enforce_constraints A function to enfoce the constraints, if any,
#' on the parameter vector. \cr
#' This will be called as enforce_constraints(param, ...) from
#' iSQUAREM function.
#' @param badness_of_fit A function to evaluate the badness of fit of
#' the current parameter vector.
#' This will be called as badness_of_fit(param, ...) from
#' iSQUAREM function.
#' @param SQUAREM = 1 or 3
#' @param minalpha The minimul value of alpha parameter of iSQUAREM
#' @param maxalpha <= 1 The maximum value of alpha parameter of iSQUAREM
#' @param bof_value current value of the badness of fit prior to iSQUAREM
#' \cr Usually, it is not necessary to specify this.
#' @param always = 1 to update parameters regardless of the badness of fit.
#' @param reset1 = 0 or 1 counter reset value after update
#' @param reset2 = 1 or 2 counter reset value after failure when always=0
#' @param print = 1 to print the update process
#' @param debug = 1 to print the generated function
#'
#'
#' @details This is a dummy version of iSQUAREM function only to
#' display its help. \cr
#' The real one must be generated using generate_iSQUAREM as a closure
#' in the function which uses the accelaration process.
#' \cr\cr
#'  How to use iSQUAREM: \cr
#'  Try always=1 with maxalpha=-1 or less first. \cr
#'  Changing to reset1=0 and reset2=1 or increasing nSQUAREM may help. \cr
#'  If it seems not working, use always=0 with maxalpha=1 or less. \cr
#'  Changing to reset1=0 and reset2=1 may help. \cr
#'  If all of the above fail, be patient and use SQUAREM=0.
#' \cr\cr
#' When always=1 the badness of fit will not be evaluated in this function.
#' \cr However, note
#' that, even if always=0, the badness of fit evaluated in
#' the calling function may not improve
#' if iSQUAREM update is accepted with badness1 < badness
#' where badness1 is the badness of fit evaluated at new update.
#' \cr
#' This happens because badness in this function
#' is the value evaluated at the last  SQUAREM update
#' and not the one evaluated with the latest param value
#' in the calling function.
#'
#' The following is an example of a function which uses iSQUAREM.
#'
#' \preformatted{
#'
#' my_algorithm <- function( data, param, etc ){
#'
#'  # function definitions
#'  const <- function( param, additional_param_list ){
#'   # enforce constraints on param
#'   # All the variables except those listed as formal arguments
#'   # reffer to the ones defined in the environmend in which
#'   # this function is defined: namely environment(my_algorithm).
#'   return( param )
#'  } # end of const
#'  fit <- function( param, additional_param_list ){
#'   # calculate the badness of fit criterion of param
#'   # All the variables except those listed as formal arguments
#'   # reffer to the ones defined in the environmend in which
#'   # this function is defined: namely environment(my_algorithm).
#'   # This function can return a list whose first element is named as critval.
#'   return( crit )
#'  } # end of fit
#'
#'  # initial value of parameter vector
#'  param <- initial_value
#'  converged <- 0
#'
#'  ################ new block of code added 1 #########################
#'  # generate iSQUREM function as a closure
#'  iSQUAREM <- generate_iSQUAREM( param )
#'  ####################################################################
#'
#'  for( iter in 1:maxiter ){
#'
#'   # usual parameter update process satisfying constraints.
#'   param <- update_param( param )
#'
#'   ################ new block of code added 2 #########################
#'   if( SQUAREM > 0  &  nSQUAREM <= iter ){
#'    # acceleration
#'    temp <- iSQUAREM( param, additional_param_list
#'                  , enforce_constraints=const, badness_of_fit=crit
#'                  , SQUAREM=SQUAREM )
#'    param <- temp$param
#'    critval <- temp$critval
#'   } # end of iSQUAREM
#'   ####################################################################
#'
#'   # check convergence in terms of critval or changes of param value.
#'   if( converged ) break
#'
#'  } # end of iterations
#'
#'  return( param )
#'
#' } # end of my_algorithm
#' }
#'
#' @return
#' A list of \cr
#' param  A vector or parameters (updated or same as the input) \cr
#' critval The value of criterion evaluated at param.
#'
#'
#' @references
#' R Varadhan and C Roland (2008),
#' Simple and globally convergent numerical schemes for acceler-
#' ating the convergence of any EM algorithm,
#' Scandinavian Journal of Statistics, 35:335-353.\cr
#' C  Roland,  R  Varadhan,  and  CE  Frangakis  (2007),
#' Squared  polynomial  extrapolation  methods
#' with cycling: an application to the positron emission tomography problem,
#' Numerical Algorithms 44:159-172.
#'
#' @return
#' A list of \cr
#' param The updated parameter value with constraints if any. \cr
#' critval The value of badness_of_fit function when always=0 or NA.
#'
#' @export
#'

iSQUAREM <- function( param, ...
                      , enforce_constraints=NULL, badness_of_fit=NULL
                      , SQUAREM=3, minalpha=-999, maxalpha=-1, bof_value=NULL
                      , always=0, reset1=1, reset2=2
                      , print=0, debug=0 ){
 # inline SQUAREM
 # This function stores the current param vector in param_save
 # and returns the updated param vector in every 3 call.
 #
 # Do NOT copy and pasete this function:
 # This function must be created by generate_iSQUAREM as a closure.
 #

 if( debug ){
  Print("top of iSQUAREM",param)
 }

 # save history
 #
 llll2 <<- llll2+1
 param_save[,3] <<- param_save[,2]          # 3 = theta_n (oldest)
 param_save[,2] <<- param_save[,1]          # 2 = F(theta_n)
 param_save[,1] <<- param                   # 1 = F(F(theta_n))

 if( !is.null(bof_value) ) badness0 <<-bof_value

 # update every three iterationss
 if( llll2 %% 3 == 0 ){

  # calculate alpha
  nSQEM <<- nSQEM+1
  r=param_save[,2] - param_save[,3]
  v=param_save[,1] - 2*param_save[,2] + param_save[,3]
  if( SQUAREM == 3 ) alpha=-sqrt(sum(r*r)/sum(v*v))
  else alpha=sum(v*r)/sum(v*v)

  if( alpha < minalpha ) alpha=minalpha
  if( alpha > maxalpha ) alpha=maxalpha

  # update parameter
  newparam=param_save[,3] -2*alpha*r + (alpha^2)*v

  # constraints
  if(!is.null(enforce_constraints) )
   newparam=enforce_constraints( newparam, ... )

  # evaluate badness of fit
  if( !always ){
   bf=badness_of_fit( newparam, ... )
   if( is.list(bf) ) badness1=bf$critval
   else badness1=bf
   rm(bf)
  }

  # Print(badness0,  badness1)

  # check if improved
  # badness0 is the value of he badness_of_fit function
  # evaluated with the latest SQUAREM update in this function,
  # or it is supplied by the bof_value argument.

  if( always  ||  badness1 <= badness0 ){

   # improved
   param_save[,1] <<- newparam
   badness0 <<- badness1

   # llll2 <<- 0      # counting next update as the first one.
   llll2 <<- reset1   # counting this update as the first one.
   nSQEMs <<- nSQEMs+1

   if( print ){
    if( !always ) cat("SQUAREM success. ")
    else  cat("SQUAREM employed. ")
    Print( badness1, alpha, fmt=c(".7"))
   }

   res=list( param=newparam, critval=badness1 )
   return( res )

  } # end of success
  else{

   # not improved: do nothing
   # llll2 <<- 1  #
   llll2 <<- reset2  # persistent setting

  } # failure

  res=list( param=param, critval=badness1 )
  return( res )


 } # end of mod=3
 else{

  res=list( param=param, critval=badness1 )
  return( res )

 }

} # end of iSQUAREM function

