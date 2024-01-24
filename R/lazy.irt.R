#' lazy.irt: Some IRT functions for lazy boys and girls
#'
#' The following are the description of the item parameter data frame and
#' the item weight data frame, and the list of functions.\cr
#' Japanese help file: (\link[lazy.irt]{lazy.irt_JPH})
#'
#'
#' @section Item Parameter Data Frame and Item Weight Data Frame:
#' \itemize{
#'  \item Structure of Item Parameter Data Frame \cr
#'  Each row of the Item Parameter Data Frame corresponds to an item.
#' The data frame must have the following variables.
#'  \preformatted{
#'  name item name
#'  type item type  B | B3 | Bn | Bn3 | G | Gn | P 
#'  ncat # of categories
#'  p1 item parameter 1 (discrimination)
#'  p2 item parameter 2 (difficulty)
#'  p3 item parameter 3 (diffuculty or asymptote)
#'  ::
#'  ::
#'  }
#'  For \code{type = "B"  or "B2"  or  "B3"   or   "Bn" or "Bn3"} items:
#'  \preformatted{
#'     p1 discrimination \eqn{a_j}
#'     p2 difficulty \eqn{b_j}
#'     p3 lower asymptote \eqn{c_j} or 0
#'     }
#'  For \code{type = "B"  or "B2"  or  "B3"   or   "Bn" or "Bn3"} items:
#'  \preformatted{
#'     p1 discrimination \eqn{a_j}
#'     p2 difficulty \eqn{b_j}
#'     p3 lower asymptote \eqn{c_j} or 0
#'     }
#' For \code{type = "G" or "Gn"} items:
#' \preformatted{
#'     p1 discrimination \eqn{a_j}
#'     p2 threshold 1 \eqn{b_{j1}}
#'     p3 threshold 2  \eqn{b_{j2}}
#'     ::       ::
#'     p_ncat thresold ncat-1  \eqn{b_{j,[ncat-1]}}
#'     }
#' For \code{type = "P"} items:
#' \preformatted{
#'     p1 discrimination \eqn{a_j}
#'     p2 step 1 \qen{b_{j1}}
#'     p3 step 2 \qen{b_{j2}}
#'      ::     ::
#'     p_ncat step ncat-1  \eqn{b_{j,[ncat-1]}}
#'     }
#' For \code{type = "PN"} items:
#' \preformatted{
#'     p1 slope \eqn{a_j}
#'     p2 intercept 1 \eqn{b_{j1}}
#'     p3 intercept 2 \eqn{b_{j2}}
#'      ::     ::
#'     p_ncat intercept ncat-1 \eqn{b_{j,[ncat-1]}}
#'     }
#' For \code{type = "N"} items:
#' \preformatted{
#'     p1through p_(ncat[j]-1) slope parameters
#'     p_ncat[j] through p_2*(ncat[j]-1) intercept parameters
#'     }
#' Binary items have three parameters, namely, discrimination, difficulty
#' and asymptote which is set equal to 0 for B or Bn items.
#' Ordered polytomous items (P, G, Gn) have \code{ncat[j]} parameters, and,
#' Nomimal Response itms (N) have \code{2*(ncat[j]-1)} parameters.
#' \cr
#'  See the test data section below.
#'
#'  \item Structure of Item Weight Data Frame \cr
#'  Each row of the Item Weight Data Frame corresponds to an item.
#' The data frame must have the following variables.
#'  \preformatted{
#'  name item name
#'  type item type B | B3 | Bn | Bn3 | G | Gn | P
#'  ncat # of categories
#'  w weight to the item
#'  v0 category weight 0
#'  p1 category weight 1
#'  p2 category weight 2
#'  ::
#'  ::
#'  }
#'  See the test data section below.
#'
#' }
#'
#'
#'
#' @section Parameter Estimation and Equating:
#' \itemize{
#'  \item \link[lazy.irt]{uIRT}:	Item Parameter Estimation of Unidimensional IRT
#'  \cr\cr
#'  \item \link[lazy.irt]{est_theta}: Estimation of Theta
#'  \cr\cr
#'  \item \link[lazy.irt]{cala}:	Calibration of Independently Estimated Sets of Item Parameters
#'  by Minimizing the LS Criterion Defined in terms of Parameter Values
#'  \item \link[lazy.irt]{calr}:	Calibration of Independently Estimated Sets of Item Parameters
#'  by Minimizing the LS Criterion Defined in terms of Probability Values
#'  \cr\cr
#'  \item \link[lazy.irt]{tseq}: IRT True Score Equating
#'  \item \link[lazy.irt]{oseq}: IRT Observed Score Equating
#'  \cr
#'  \item \link[lazy.irt]{coseq}: Classical Observed Score Equating
#'  \cr\cr
#'  \item \link[lazy.irt]{ordinal_reg}: Oridinal Regression (used in the m-step of \code{uIRT}.)
#'  \item \link[lazy.irt]{smn}: Parameter Estimation of Parametric (Normal pdf)
#'  Scored Multinomial Distributions
#'  \item \link[lazy.irt]{invtrf}: The inverse function of trf (used in \code{tseq}.)
#'  }
#'
#' @section Item Response Functions:
#' \itemize{
#' \item \link[lazy.irt]{irf}:	Calculation of Item Response Function
#' \cr\cr
#' \item \link[lazy.irt]{icrfB}: Calculation of Item Response Function of Binary Logistic items
#' \item \link[lazy.irt]{icrfG}: Calculation of Item Response Function of Graded Response items
#' \item icrfN: Calculation of Item Response Function of Nominal items (N.A.)
#' \item \link[lazy.irt]{icrfP}: Calculation of Item Response Function of Partial Credit items
#' \item \link[lazy.irt]{icrfPN}: Calculation of Item Response Function of Partial Credit items
#' in Nominal Model Format
#' \item \link[lazy.irt]{icrfPN0}: Calculation of Partial Credit ICRF in Nominal format
#' }
#'
#'
#' @section Derivative of Item Response Functions:
#' \itemize{
#'  \item \link[lazy.irt]{dirf}: Calculation of Derivative of Item Response Function
#'  \cr\cr
#'  \item \link[lazy.irt]{dicrfB}: Calculation of the Derivative of
#'   Item Response Function of Binary Logistic Items
#'  \item \link[lazy.irt]{dicrfG}: Calculation of the Derivative of
#'   Item Response Function of Graded Response Items
#'  \item \link[lazy.irt]{dicrfN}: Calculation of the Derivative of
#'   Item Response Function of Nominal Items  (N.A.)
#'  \item \link[lazy.irt]{dicrfP}: Calculation of the Derivative of
#'   Item Response Function of Partial Credit Items
#'  \item \link[lazy.irt]{dicrfPN}: Calculation of the Derivative of
#'   Item Response Function of Partial Credit Items in Nominal Format
#'  \item \link[lazy.irt]{dicrfPN0}: Calculation of the Derivative of
#'   Partial Credit ICRF in Nominal format
#'   \cr
#'  \item \link[lazy.irt]{dicrf_num}: Numerical Derivative of Item Response Function
#'     using JacobianMat
#'   \cr\cr
#'  \item \link[lazy.irt]{dirf_p}: Calculation of the Derivative of Item Response Function
#'   with respect to Item Parameters
#' }
#'
#' @section Observed Score Distributions:
#' \itemize{
#'  \item \link[lazy.irt]{obscore}: Calculation of Observed Score Distribution and Posterior
#'   Distribution of Theta with Various Information Functions
#'  \item \link[lazy.irt]{obscore_s}: Calculation of Observed Score Distribution
#'  (simple version of \code{obscore}.)
#'  \cr\cr
#'  \item \link[lazy.irt]{sumsmnw}: Distribution of the Weighted Sum of
#'  Several Independent Scored Multinomial Distributions
#'  \item \link[lazy.irt]{sumsmnw12}: Distribution of the Weighted Sum of
#'  Two Independent Scored Multinomial Distributions
#'  \cr\cr
#'  \item \link[lazy.irt]{rel_irt}: Calculation of Test Reliability and Average SEM
#'  under IRT model
#'
#' }
#'
#'
#' @section Information Functions:
#' \itemize{
#'  \item \link[lazy.irt]{iif}: Calculation of Information Functions
#'  \cr
#'  \item \link[lazy.irt]{info_func}: Calculation of Various Information Functions
#'  \item \link[lazy.irt]{graded_info}: Calculation of the Information Function
#'  associated with the Graded Observed Score.
#'  \item \link[lazy.irt]{flatten_SEM}: Find the Transformation of the Observed Score
#'  such that the Resulting Score has a Flat SEM almost everywhere.
#'  \item \link[lazy.irt]{flatten_SEM_theta}: Find the Transformation of
#'  the Thetahat Based Observed Score
#'  such that the Resulting Score has a Flat SEM almost everywhere.
#'  \cr\cr
#'  \item \link[lazy.irt]{GOptWeight} Estimation of the Globally Optimal Item Category Weights
#'  }
#'
#'
#' @section Conversion Functions:
#' \itemize{
#'  \item \link[lazy.irt]{fitG2P}: Approximate Conversion of Partial Credit Items to
#'  Graded Response Items Using Logit Transformation
#'  \item \link[lazy.irt]{fitP2G}: Approximate Conversion of Graded Response Items to
#'  Partial Credit Items Using Logit Transformation
#'  \item \link[lazy.irt]{conv2G}: Conversion to  Logistic Graded Response Model
#'  \item \link[lazy.irt]{conv2Gn}: Conversion to Normal Graded Response Model
#'  \item \link[lazy.irt]{conv2P}: Conversion to Partial Credit Model
#'  }
#'
#'
#' @section Reading Item Paramter File:
#' \itemize{
#'  \item \link[lazy.irt]{read.param}: Reading Parameter File
#'  \item \link[lazy.irt]{read.weight}: Reading Weight File
#'  \cr
#'  \item \link[lazy.irt]{read_blg_par}: Read Bilog .par file.
#'  }
#'
#'
#' @section Utility Functions:
#' \itemize{
#'  \item \link[lazy.irt]{create_weight_df}: Creation of Weight Data Frame
#'  \item \link[lazy.irt]{find_minmax_score}: Find the min and max score from Weight Data Frame
#'  \item \link[lazy.irt]{checkparam2}: Checking Parameter Data Frame (New Version)
#'  \item \link[lazy.irt]{gendataIRT}: Generation of Simulated Item Response Data
#'  \item \link[lazy.irt]{find_mode}: Find the Mode of icrf
#'  \item \link[lazy.irt]{find_intersection}: Find the Intersection of icrfs
#'  \item \link[lazy.irt]{gen_icrfnames}: Generate the row names of vec(icrf)
#'  \item \link[lazy.irt]{graded_prob}; Calculation of Probability Contents
#'  associated with Graded Score
#'  \cr\cr
#'  \item \link[lazy.irt]{convP2N}: Convert Standard Partial Credit Item Parameters
#'  (step parameters) to in Nominal Format
#'  \item \link[lazy.irt]{convP2PN}: Convert Standard Partial Credit Item Parameters
#'  (step parameters) to in Nominal Format
#'  Standard Format to Partial Credit Parameters in Nominal Format
#'  \item \link[lazy.irt]{convPN2P}: Convert Partial Credit Item Parameters in
#'  Nominal Format to Standard Format (step parameters)
#'  }
#'
#'
#' @section Data:
#' \itemize{
#'  \item \link[lazy.irt]{paramA1}: Set of All Types of Items # 1. (10 items)
#'  \cr
#'  \item \link[lazy.irt]{paramB1}: Binary Item Parameter Data Frame # 1.
#'   (18 2PLM items)
#'  \item \link[lazy.irt]{paramB2}: Binary Item Parameter Data Frame # 2.
#'   (18 3PLM items)
#'  \cr\cr
#'  \item \link[lazy.irt]{paramS1}: Small Item Parameter Data Frame # 1.
#'  (3 mixed type items)
#'  \item \link[lazy.irt]{paramS2}: Small Item Parameter Data Frame # 2.
#'  (8 mixed type items)
#'  \item \link[lazy.irt]{paramS3}: Small Item Parameter Data Frame # 2.
#'  (20 mixed type items)
#'  \cr\cr
#'  \item \link[lazy.irt]{weightA1}: Item Weight Data Frame for paramA1.
#'  Natural category weight.
#'  \cr\cr
#'  \item \link[lazy.irt]{weightB1}: Item Weight Data Frame for
#'  paramB1 or paramB2 1. Natural weight.
#'  \item \link[lazy.irt]{weightB11}: Item Weight Data Frame for
#'  paramB1 or paramB2 2. Large weight for difficult items.
#'  \item \link[lazy.irt]{weightB12}: Item Weight Data Frame for
#'   paramB1 or paramB2 3. Large weight for easy items.
#'  \item \link[lazy.irt]{weightB13}: Item Weight Data Frame for
#'  paramB1 or paramB2 4. Large weight for high discriminative items.
#'  \item \link[lazy.irt]{weightB14}: Item Weight Data Frame for
#'   paramB1 or paramB2 5. Large weight for low discriminative items.
#'  \cr\cr
#'  \item \link[lazy.irt]{weightS1}: Small Item Weight Data Frame # 1
#'  To be used in conjunction with paramS1.
#'  \item \link[lazy.irt]{weightS11}: Small Item Weight Data Frame # 2.
#'  To be used in conjunction with paramS1.
#'  \item \link[lazy.irt]{weightS12}: Small Item Weight Data Frame # 3.
#'  To be used in conjunction with paramS1.
#'  \item \link[lazy.irt]{weightS2}: Small Item Weight Data Frame # 4.
#'  To be used in conjunction with paramS2.
#'  \item \link[lazy.irt]{weightS21}: Small Item Weight Data Frame # 5.
#'  To be used in conjunction with paramS2.
#'  \item \link[lazy.irt]{weightS22}: Small Item Weight Data Frame # 6.
#'  To be used in conjunction with paramS2.
#'  \item \link[lazy.irt]{weightS3}: Small Item Weight Data Frame # 7.
#'  To be used in conjunction with paramS3.
#'  \cr\cr
#'  \item \link[lazy.irt]{paramCal1}: Small Item Parameter Data Frame for
#'  Testing cala and calr #1.
#'  \item \link[lazy.irt]{paramCal2}: Small Item Parameter Data Frame for
#'  Testing cala and calr #2.
#'  }
#'
#'
#' @section LRT:
#' \itemize{
#'  \item \link[lazy.irt]{uLRT}: Item Parameter Estimation of Unidimensional LRT \cr
#'  with strict monotonicity restrictions.
#'  \item \link[lazy.irt]{est_rank}: Estmation of LRT latent rank for each person.
#'  \item \link[lazy.irt]{fitI2L}: Approxiamte Conversion of LRT model to 2PLM IRT model
#'  \item \link[lazy.irt]{fitI2L_ls}: LS Conversion of LRT model to 2PLM IRT model
#'  }
#'
#' @docType package
#' @name lazy.irt
NULL
