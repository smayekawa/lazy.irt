#' 正規分布の密度関数に比例する確率を持つ多項分布の確率の推定
#'
#' 数値的カテゴリを持つ多項分布の確率を、それが正規分布の密度関数に比例する
#' という仮定の下に推定する。\cr
#' English help file: (\link[lazy.irt]{smn})
#'
#' @usage smn( x, f, mu=NULL, sigma=NULL
#' ,    estmu=1, estsigma=1, maxiter=100, eps=1e-8, print=1, plot=0 )
#'
#' @param x 数値的カテゴリのベクトル
#' @param f 各カテゴリの度数
#' @param mu 正規分布の平均の初期値
#' @param sigma 正規分布の標準偏差の初期値
#' @param estmu = 0 \code{mu} を推定しない。
#' @param estsigma = 0 \code{sigma} を推定しない。
#' @param maxiter 繰り返し数の上限
#' @param eps eps 収束判定基準
#' @param print = 1 結果を出力する。
#' @param plot = 1 結果をプロットする。
#'
#' @details
#' 度数分布の形のデータ \code{(x,f)} が与えたれたときに, この関数は
#' 多項分布の対数尤度関数 \cr
#'   \eqn{ lnL = \sum_{i} f_i \log{p(x_i|mu,sigma)} } \cr
#' を最大とするような \code{mu} と \code{sigma} を求める。
#' ただし \cr
#'  \eqn{ p(x_i|mu,sigma) = c \times dnorm( x_i, mu, sigma )  } \cr
#'  and  \eqn{ c = 1 / \sum_{i} dnorm( x_i, mu, sigma ) }  \cr
#' である。
#'
#' @return  以下の要素からなるリスト
#' \preformatted{
#' P: 推定された確率のベクトル（正規分布の密度関数に比例する）
#' x: 数値的カテゴリベクトル
#' mu: mu の最尤推定値
#' sigma: sigma の最尤推定値
#' estmu:
#' estsigma:
#' llh: 最大化された対数尤度
#' mean: 標本平均
#' std: 標本標準偏差
#' }
#'
#' @examples
#' set.seed(1701)
#' n <- 1000
#' npoints <- 21
#' resg <- list(midpoints=c(-2,-1,0,1,2), freq=c(5,6,8,6,2))
#' res <- smn( resg$midpoints, resg$freq, print=1, plot=1, mu=0, sigma=1 )
#'
#' # @keywords internal
#' @docType package
#' @name smn_JPH
NULL
