#' 古典的方法によるテスト得点の等化
#'
#' この関数は古典的方法による二つのテスト得点の等化（線形等化法と等百分位法）
#' を行う。 \cr
#' English help file: (\link[lazy.irt]{coseq})
#'
#' @usage coseq( score1, freq1, cdf1=NULL, score2, freq2, cdf2=NULL
#' , lim1=NULL, lim2=NULL
#' , smooth1=0, bandwid1=3, smooth2=0, bandwid2=3
#' , method=2, interpol_method="linear"
#' , title="", nolow=0, print=1, plot=1 )
#'
#' @param score1 テスト１の得点からなるベクトル
#' @param freq1  テスト１の得点の度数からなるベクトル
#' @param cdf1   テスト１の得点の相対累積度数からなるベクトル
#' @param score2 テスト２の得点からなるベクトル
#' @param freq2  テスト２の得点の度数からなるベクトル
#' @param cdf2   テスト２の得点の相対累積度数からなるベクトル
#'
#' @param lim1 テスト１の得点の理論的最小値と最大値
#' @param lim2 テスト２の得点の理論的最小値と最大値
#'
#' @param smooth1 テスト１の cdf の移動平均法による平滑化の階数
#' @param bandwid1 移動平均法の幅
#' @param smooth2 テスト２の cdf の移動平均法による平滑化の階数
#' @param bandwid2 移動平均法の幅
#'
#'
#' @param method = 1 線形等化法 \cr
#'               = 2 等百分位法
#'
#' @param interpol_method = 補間法のオプションで "constant", "linear" or "spline"
#'
#' @param title  タイトルに使う文字列
#' @param nolow  = 1 元の点数を下げない変換を計算する
#'
#' @param print = 1 結果を出力
#' @param plot = 1 変換表のプロット \cr
#'             = 2 cdf のプロット \cr
#'             = 3 平滑化された cdf のプロット
#'
#' @details
#' テスト１の得点 \code{x1} のテスト２の得点e \code{x2} への等百分位法に
#' よる等化法は以下の様に定義されている。
#' \preformatted{
#'  x21 = invF2( F1(x1) )
#' }
#' ただし、 \code{F1} は \code{x1} の分布関数、\code{invF2} は \code{x2} の
#' 分布関数の逆関数である。
#' この関数では、\code{invF2} は \code{( F2(x2), x2 )} の表、
#' もしくはその平滑化されたものを \code{F1(x1)} で補間することにより求めている。
#'
#' 等百分位法は所謂 QQプロットと同値であり、以下の二つはほぼ同じ結果を返す。
#' \preformatted{
#'  coseq( score1=score1, freq1=freq1, score2=score2, freq2=freq2 )
#'  qqplot( expand_freqdist( score1, freq1 )
#'              , expand_freqdist( score2, freq2 ), type="l" )
#' }
#' なお、R の \code{qqplot} は度数分布の形式のデータを取り扱えないので
#' \code{lazy.tools::expand_freqdist} を用いて度数分布を素データに変換している。
#'
#' \code{cdf} の指定が \code{freq} の指定に優先する。
#'
#'
#'
#' @return 以下の要素を持つリスト
#' \preformatted{
#' ctable: \code{(score, score21, freq1)} からなる整数化された変換表
#' ctable0: 整数化される前の変換表
#' mands:  変換された得点の統計量
#' newfreq: 変換された得点の度数分布表 \code{(score21, freq21)}
#' sdist1 and sdist2: テストの得点分布（平滑化の結果を含む）
#' cntr: この関数へのパラメタのリスト
#' }
#'
#' 出力される \code{ctable} は
#' テスト１の \code{score1[i]} 点がテスト２の \code{score21[i]} 点に
#' 等しいことを示している。
#'
#'
#' @examples
#'
#' set.seed(1701)
#'
#' x1 <- round( 10*rbeta(500,2,4) )
#' f1 <- table(x1)
#' x1 <- as.numeric(names(f1))
#' x2 <- round( 15*rbeta(1000,4,2) )
#' f2 <- table(x2)
#' x2 <- as.numeric(names(f2))
#'
#' reseq <- coseq( score1=x1, freq1=f1, score2=x2, freq2=f2
#'                 , lim1=c(0,10), lim2=c(0,15)
#'                 , smooth1=3, smooth2=3, method=2, plot=0 )
#'
#' # @keywords internal
#' @docType package
#' @name coseq_JPH
NULL
