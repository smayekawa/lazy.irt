#' テスト得点の分布の計算と特性値の事後分布の計算
#'
#' 項目反応理論に基づき、与えられた項目パラメタと特性値の分布から
#' テスト得点の分布と重み付きテスト得点を与えたときの
#' 特性値の事後分布を計算する。 \cr
#' English help file: (\link[lazy.irt]{obscore})
#'
#' @usage obscore( param, weight=NULL
#' , npoints=31, thmin=-4, thmax=4, thdist=1
#' , alpha=0.1, compress=0, print=1, plot=0, debug=0 )
#'
#'
#' @param param 項目パラメタデータフレーム
#' @param weight 項目重みデータフレーム \cr
#' または項目重みからなるベクトル
#' @param npoints 特性値 theta の離散点のベクトル
#' @param thmin theta の最小値
#' @param thmax theta の最大値
#' @param thdist 特性値の分布の種類 \cr
#' = 0 一様分布,  \cr = 1  N(0,1)
#' @param alpha 信頼区間の確率（1-信頼水準）
#' @param compress = 1 確率が 0 となる得点を除く
#' @param print = 1 結果の出力
#' @param plot = 1 グラフの出力
#' @param debug = 1 中間結果の出力
#'
#'
#' @return 以下のリスト
#' \cr
#'   list( theta_stat, obs_stat, Px_t, Pt_x,  etcetc )\cr
#' \cr
#'    ただし\cr
#' \cr
#'   theta_stat は以下の変数を持つデータフレーム:\cr
#'     theta:       theta points\cr
#'     Pt:          prior probability distribution of theta\cr
#'     TRF:         test response function\cr
#'     slope_TRF:   slope of TRF\cr
#'     stdx_t:      standard deviation of X (observed score) given theta\cr
#'     info:        information function defined as
#'     (slope_TRF)^2 / (stdx_t)^2\cr
#'     info_LOW:    information function with the locally optimal item weight
#'     given category weights\cr
#'     info_LO:     information function with the locally optimal category
#'      weight\cr
#'     qt_L:        upper quantile of X given theta\cr
#'     qt_U:        lower quantile of X given theta\cr
#'     poststd:     posterior std of theta given X\cr
#'                    as a function of posterior mean\cr
#' \cr
#'   obs_stat は以下の変数を持つデータフレーム:\cr
#'     score:      domain of X (observed score)\cr
#'     Px:         marginal probability  of X\cr
#'     post_mean:  posterior mean of theta given X\cr
#'     post_std:   posterior standard deviation of theta given X\cr
#'     ci_L:       upper limit of conficence interval for theta given X\cr
#'     ci_U:       lower limit of conficence interval for theta given X\cr
#'     ci_hwid:    half the width of CI\cr
#' \cr
#'   SigmaX2, SigmaT2, SigmaE2:　観測得点の分散、真値の分散、誤差の分散の平均
#' \cr
#'   aSEM=sqrt(SigmaE2): （平均的）測定の標準誤差, \cr
#'   rel: 信頼性係数
#' \cr
#'   SigmaT_theta2, SigmaE_theta2, aSEM_theta, rel_theta
#' \cr
#'   Px_t           score x theta  theta を与えたときの X の分布 \cr
#'   Pt_x           score x theta  X を与えたときの theta の（事後）分布 \cr
#'   Px             score x 1      X の周辺分布 \cr
#' \cr
#' pdfname: item parameter data frame name \cr
#' wdfname: item weight data frame name \cr
#' npoints: # of theta points \cr
#' thmin, thmax: the range of theta \cr
#' thdist: Type of theta distribution \cr
#' nitems: # of items \cr
#' minscore_t: minimum score \cr
#' maxscore_t: maximum score \cr
#' alpha: small probability value \cr
#'
#' @details
#' この関数は、重み付テスト得点 X の分布（特性値 theta が所与の場合の
#' 条件付き分布、特性値と X との同時分布、ならびに X の周辺分布）を計算する。
#' また、テスト得点 X を与えたときの特性値の事後分布とその平均ならびに標準偏差
#' を計算する。
#' ただし、特性値を与えたときの X の条件付分布は
#' \link[lazy.irt]{sumsmnw_JPH} を用いて計算される。
#' \cr
#' 同時にに、重み付得点から特性値を推定する際の情報関数と、２種類の局所的
#' 最適な重みを用いた際の情報関数を出力する。
#'
#' 項目パラメタデータフレームと項目重みデータフレームに関しては
#' \link[lazy.irt]{lazy.irt_JPH} を参照のこと。
#'
#' 項目カテゴリならびに項目への重みが所与の場合、重み付得点 X から特性値 theta
#' を推定する際の情報関数は以下のように定義される。\cr
#' \code{
#'  ( slope of trf at theta )^2 / (variance of X at theta)
#'  }\cr
#' ただし、trf は特性値 theta を与えたときの重み付得点 X の条件付き期待値、
#' すなわち Test Information Function (trf) である。
#' \cr
#' 重みを変化させることにより上記の方法関数の値は変化するが、一般的に、
#' 情報関数を最大とする重みは特性値の値に依存することが知られている。
#' したがって、このような最適重みを局所的最適重みと呼ぶ。
#'
#' 項目カテゴリの重みが所与の場合に上記の情報関数を最大とする項目への重みは \cr
#' \eqn{ w_j(\theta)= P^{*'}_j(\theta) / var( U_j^{*} | \theta) } \cr
#' ただし \cr
#' \eqn{ U_j^{*} = \sum_k v_{kj} U_{kj}} and
#' \eqn{ P^{*'}_j(\theta)=\sum_k v_{kj} P_{kj}(\theta) }.
#' \cr
#' で与えられ、局所的な最適項目重み Locally Optimal Item Weight (LOW)
#' と呼ばれる。
#'
#' また、情報関数を最大とする項目カテゴリの重みは局所的最適重み Locally
#' Optimal Weight (LO) と呼ばれるが、それは以下の式で与えられる。 \cr
#' \eqn{ v_{kj}(\theta)=P'_{kj}(\theta) / P_{kj}(\theta)
#'  - P'_{0j}(\theta) / P_{0j}(\theta) }.
#' \cr
#' この LO は Samejima(1969) の the basic function と関連している。
#'
#' LO を用いた最大化された情報関数は \cr
#' \eqn{ \sum_j \sum_k (P'_{kj}(\theta))^2 / P_{kj}(\theta) } \cr
#' ただし \eqn{P_{kj}(\theta)} は項目カテゴリ情報関数 icif であり
#' \eqn{P'_{kj}(\theta)} はその導関数であるが、
#' これを上回る情報関数は存在しない。
#' \cr
#' LOW を用いた情報関数は \cr
#' \eqn{ \sum_j  (P_j^{*'}(\theta))^2 / var(U_j^{*} | \theta) } \cr
#' ただし \eqn{U_j^{*} = \sum_k v_{kj} U_{kj}} は重み付き得点であり
#' \eqn{P_j^{*'}(\theta)} は特性値を与えたときの重み付得点の期待値、すなわち
#' 項目反応関数 irf、の導関数である。
#'
#' \cr
#' 信頼性係数は以下のように定義される。 \cr
#' \code{rel = (SigmaX^2-SigmaE^2)/SigmaX^2}, \cr
#'
#'
#' @references
#' Birnbaum, A.(1968) Some Latent Traint Models.
#' In F. M. Lord and M. R. Novick, Statistical Theories of Mental Test Scores.
#'  Reading, Mass.: Addison-Wesley. \cr
#'
#' Mayekawa, S. (2018) 重み付き合計点を用いた特性値 θ の
#' 推定方法と推定の標準誤差の最小化.
#' 大学入試センター研究開発部ルサーチノート RN-18-01. \cr
#'
#' Mayekawa, S., & Arai, S. (2008).
#' Distribution of the Sum of Scored Multinomial Random Variables
#' and Its Application to the Item Response Theory.
#' In K. Shigemasu, A. Okada, T.Imaizumi, & T. Hoshino (Eds.)
#' New Trends in Psychometrics. Tokyo: University Academic Press. \cr
#'
#' Samejima,  F.  (1969). Estimation  of  a  latent  ability  using  a
#' response  pattern  of  graded  scores. Psychometrika  Monographs, 34
#' (Suppl. 4).
#'
#'
#' @examples
#' # Define the observed raw score X
#' # and calculate the score distribution, information function, etc.
#' res <- obscore( paramS1, plot=1 )
#'
#' # Define X using the item weights w stored in weightS11
#' res <- obscore( paramS1, weight=weightS11, plot=1 )
#'
#' # Define X using the item and category weights w and v stored in weightS12.
#' res <- obscore( paramS1, weight=weightS12, plot=1 )
#'
#' # @keywords internal
#' @docType package
#' @name obscore_JPH
NULL
