#' 反応関数（特性曲線）の計算
#'
#' この関数は、テスト反応関数 (trf)、項目反応関数 (irf)、ならびに
#' 項目カテゴリ反応関数 (icrf) を計算する。\cr
#' English help file: (\link[lazy.irt]{irf})
#'
#' @usage irf( param, theta=NULL, weight=NULL, zero=1, smallP=1e-9
#' , thmin=-4, thmax=4, npoints=21, DinP=1
#' , print=1, debug=0
#' , plot=0, colors="black", plotchar=NULL, linetype=NULL )
#'
#' @param param 項目パラメタデータフレーム
#' @param theta 特性値 theta の離散点
#' @param weight 項目重みデータフレーム
#' @param zero = 0 第 0 番目のカテゴリの出力をしない。
#' @param smallP 確率の最小値
#' @param thmin 特性値 theta の最小値
#' @param thmax 特性値 theta の最大値
#' @param npoints 離散点の数
#' @param DinP = 1 ロジスティク関数に D=1.7 を用いる。
#' @param print = 1 結果を出力する
#' @param debug = 1 中間結果を出力する
#' @param plot = 1 グラフを出力する
#' @param colors = グラフの線の色
#' @param plotchar = グラフの点の種類
#' @param linetype = グラフの線の種類
#'
#' @details
#' 項目パラメタデータフレームと項目重みデータフレームに関しては
#' \link[lazy.irt]{lazy.irt_JPH} を参照のこと。
#'
#' 多値項目の項目反応関数 irf は以下のように定義されるが、それは、特性値 theta
#' を与えたときの重み付き項目得点の期待値である。 \cr
#' \eqn{ E(\sum_k v_{kj} U_{kj} | \theta) = \sum_{k} v_{kj} P_{kj}(\theta)}
#' \cr
#'
#'実際の計算に際しては以下の関数群が用いられている。\cr
#'  \code{icrfB, icrfBN, icrfG, idrfGn, icrfPN, icrfP, icrfN, checkparam }
#'
#' @return　以下の要素からなるリスト \cr
#'    \code{ list( ICRF, IRF, TRF, fromP, toP ) }
#' \cr
#'   ただし、
#'   \preformatted{
#'ICRF  npoints x sum(ncat) 項目カテゴリ反応関数 (icrf)
#'IRF npoints x nitems   項目カテゴリ重みで重みづけられた項目反応関数 (irf)
#'TRF npoints x 1 項目カテゴリ重みと項目重みで重みづけられたテスト反応関数 (trf)
#'fromP, toP  ICRF の列中の各項目（のカテゴリ）の位置
#'      }
#'
#' @examples
#' Cards="
#' name type ncat   p1   p2   p3   p4
#'   Q1    B    2  0.9  0.0,,
#'   Q2    B3   2  1.0  0.0  0.2,
#'   Q3    G    4  1.0 -2.0  0.0  2.0
#'   Q4    P    4  0.8 -2.0  0.0  2.0
#' "; param <- cards( Cards, header=1 )
#' irf( param, plot=1 )
#'
#' # @keywords internal
#' @docType package
#' @name irf_JPH
NULL
