#' lazy.irt: Some IRT functions for lazy boys and girls
#'
#'
#' The following are the classes of the functions.
#'
#'
#' @section Item Parameter Estimation and Equating:
#' \itemize{
#'  \item uIRT:	Item Parameter Estimation of Unidimensional IRT
#'  \cr\cr
#'  \item cala:	Calibration of Independently Estimated Sets of Item Parameters
#'  by Minimizing the LS Criterion Defined in terms of Parameter Values
#'  \item calr	Calibration of Independently Estimated Sets of Item Parameters
#'  by Minimizing the LS Criterion Defined in terms of Probability Values
#'  \cr\cr
#'  est_theta: Estimation of Theta
#'  }
#'
#' @section Item Response Functions:
#' \itemize{
#' \item irf	Calculation of Item Response Function
#' \cr\cr
#' \item icrfB: Calculation of Item Response Function of Binary Logisticitems
#' \item icrfG: Calculation of Item Response Function of Graded Responseitems
#' \item icrfN: Calculation of Item Response Function of Nominalitems (N.A.)
#' \item icrfP: Calculation of Item Response Function of Partial Credititems
#' \item icrfPN: Calculation of Item Response Function of Partial Credititems
#' in Nominal Model Format
#' \item icrfPN0: Calculation of Partial Credit ICRF in Nominal format
#' }
#'
#'
#' @section Derivative of Item Response Functions:
#' \itemize{
#'  \item dirf: Calculation of Derivative of Item Response Function
#'  \cr\cr
#'  \item dicrfB: Calculation of the Derivative of
#'   Item Response Function of Binary Logistic Items
#'  \item dicrfG: Calculation of the Derivative of
#'   Item Response Function of Graded Response Items
#'  \item dicrfN: Calculation of the Derivative of
#'   Item Response Function of Nominal Items  (N.A.)
#'  \item dicrfP: Calculation of the Derivative of
#'   Item Response Function of Partial Credit Items
#'  \item dicrfPN: Calculation of the Derivative of
#'   Item Response Function of Partial Credit Items in Nominal Format
#'  \item dicrfPN0: Calculation of the Derivative of
#'   Partial Credit ICRF in Nominal format
#'   \cr
#'  \item dicrf_num: Numerical Derivative of Item Response Function
#'     using JacobianMat
#'   \cr\cr
#'  \item dirf_p: Calculation of the Derivative of Item Response Function
#'   with respect to Item Parameters
#' }
#'
#' @section Observed Score Distributions:
#' \itemize{
#'  \item obscore: Calculation of Observed Score Distribution and Posterior
#'   Distribution of Theta with Various Information Functions
#'  \cr\cr
#'  \item smn: Parameter Estimation of Parametric (Normal pdf)
#'  Scored Multinomial Distributions
#'  \item sumsmnw: Distribution of the Weighted Sum of
#'  Several Independent Scored Multinomial Distributions
#'  \item sumsmnw12: Distribution of the Weighted Sum of
#'  Two Independent Scored Multinomial Distributions
#'  \cr\cr
#'  \item rel_irt: Calculation of Test Reliability and Average SEM
#'  under IRT model
#'
#' }
#'
#'
#' @section Information Functions:
#' \itemize{
#'  \item info_func: Calculation of Various Information Functions
#'  \item graded_info: Calculation of the Information Function
#'  associated with the Graded Observed Score.
#'  \item flatten_SEM: Find the Transformation of the Observed Score
#'  such that the Resulting Score has a Flat SEM almost everywhere.
#'  \item flatten_SEM_theta: Find the Transformation of
#'  the Thetahat Based Observed Score
#'  such that the Resulting Score has a Flat SEM almost everywhere.
#'  }
#'
#'
#' @section Conversion Functions:
#' \itemize{
#'  \item fitG2P: Approximate Conversion of Partial Credit Items to
#'  Graded Response Items Using Logit Transformation
#'  \item fitP2G: Approximate Conversion of Graded Response Items to
#'  Partial Credit Items Using Logit Transformation
#'  \item fitG2P_ls: Conversion of Partial Credit Items to
#'  Graded Response Items
#'  \item fitP2G_ls: Conversion of Graded Response Items to
#'  Partial Credit Items
#'  \item fit223_ls: Conversion of 3PLM Items to 2PLM Items
#'  \cr\cr
#'  \item convP2N: Convert Partial Credit Item Parameters to in Nominal Format
#'  \item convP2PN: Convert Partial Credit Item Parameters in
#'  Standard Format to Partial Credit Parameters in Nominal Format
#'  \item convPN2P: Convert Partial Credit Item Parameters in
#'  Nominal Format to Standard Format
#'  }
#'
#'
#' @section Reading Item Paramter File:
#' \itemize{
#'  \item read.param: Reading Parameter File
#'  \item read.weight: Reading Weight File
#'  }
#'
#'
#' @section Utility Functions:
#' \itemize{
#'  \item checkparam: Checking Parameter Data Frame
#'  \item gendataIRT: Generation of Simulated Item Response Data
#'  \item find_mode: Find the Mode of icrf
#'  \item find_intersection: Find the Intersection of icrfs
#'  \item ordinal_reg: Oridinal Regression (used in uIRT)
#'  }
#'
#'
#' @section Data:
#' \itemize{
#'  \item paramB1: Binary Item Parameter Data Frame # 1. (18 2PLM items)
#'  \item paramB2: Binary Item Parameter Data Frame # 2. (18 3PLM items)
#'  \cr\cr
#'  \item paramS1: Small Item Parameter Data Frame # 1. (3 mixed type items)
#'  \item paramS2: Small Item Parameter Data Frame # 2. (8 mixed type items)
#'  \item paramS3: Small Item Parameter Data Frame # 2. (20 mixed type items)
#'  \cr\cr
#'  \item weightS1: Small Item Weight Data Frame # 1
#'  To be used in conjunction with paramS1.
#'  \item weightS11: Small Item Weight Data Frame # 2.
#'  To be used in conjunction with paramS1.
#'  \item weightS12: Small Item Weight Data Frame # 3.
#'  To be used in conjunction with paramS1.
#'  \item weightS2: Small Item Weight Data Frame # 4.
#'  To be used in conjunction with paramS2.
#'  \item weightS21: Small Item Weight Data Frame # 5.
#'  To be used in conjunction with paramS2.
#'  \item weightS22: Small Item Weight Data Frame # 6.
#'  To be used in conjunction with paramS2.
#'  \item weightS3: Small Item Weight Data Frame # 7.
#'  To be used in conjunction with paramS3.
#'  \cr\cr
#'  \item paramCal1: Small Item Parameter Data Frame for
#'  Testing cala and calr #1.
#'  \item paramCal2: Small Item Parameter Data Frame for
#'  Testing cala and calr #2.
#'  }
#'
#'
#' @section LRT:
#' \itemize{
#'  \item uLRT: Item Parameter Estimation of Unidimensional LRT \cr
#'  with strict monotonicity restrictions.
#'  \item est_rank: Estmation of LRT latent rank for each person.
#'  }
#' \itemize{
#'  \item fitI2L: Approxiamte Conversion of LRT model to 2PLM IRT model
#'  \item fitI2L_ls: LS Conversion of LRT model to 2PLM IRT model
#'  }
#'
#' @docType package
#' @name lazy.irt
NULL
