#' 項目反応理論を用いた観測得点の等化 IRT Observed Score Equating
#'
#' この関数は、項目反応理論を用いて観測される重み付きテスト得点の分布を算出し、
#' それを用いてテスト得点の等化を行う。\cr
#' English help file: (\link[lazy.irt]{oseq})
#'
#' @usage oseq( param1=NULL, param2=NULL, weight1=NULL, weight2=NULL
#' , thmin=-4, thmax=4, npoints=31, thdist=1
#' , smooth1 = 0, bandwid1 = 3, smooth2 = 0, bandwid2=3
#' , round=0, print=0, plot=0 )
#'
#'
#' @param param1 テスト１ の項目パラメタデータフレーム
#' @param param2 テスト２ の項目パラメタデータフレーム
#' @param weight1 テスト１ の項目重みデータフレーム, or NULL
#' @param weight2 テスト２ の項目重みデータフレーム, or NULL
#' @param npoints 観測値の分布に用いる区間 \code{[thmin,thmax]} の特性値
#' theta の離散点の数
#' @param thmin 特性値 theta の最小値
#' @param thmax 特性値 theta の最大値
#' @param thdist = 1 特性値の分布に正規分布を仮定, = 0 一様分布を仮定.
#' @param smooth1	移動平均法による cdf1 の平滑化数.
#' @param bandwid1 cdf1 の平滑化に際する移動平均の幅.
#' @param smooth2	移動平均法による cdf2 の平滑化数.
#' @param bandwid2 cdf2 の平滑化に際する移動平均の幅
#' @param round = 1 結果を整数に丸める
#' @param print = 1 or 2 結果を出力
#' @param plot = 1 グラフを出力
#'
#'
#' @details
#' この関数は、項目パラメタと特性値 theta の分布が所与の場合に、
#' テスト１ の得点 x1 ならびにテスト２の得点 x2 の得点分布を計算し、
#' その後 x1 を x2 へ等百分位法を用いて等化する。
#'
#' 得点分布の算出には \link[lazy.irt]{obscore_JPH} の核となる部分を抽出した
#'  \link[lazy.irt]{obscore_s} 関数が用いられ、等百分位法は
#'   \link[lazy.irt]{coseq_JPH} 関数を用いて行われる。
#' また、得られた結果は \code{x2_1} と命名される。
#'
#' 項目パラメタデータフレームと項目重みデータフレームに関しては
#' \link[lazy.irt]{lazy.irt_JPH} を参照のこと。
#'
#' 以下の例が示すように、 \code{round=1} が指定された場合には、テスト１の
#' テスト２への等化と、テスト２のテスト１への等化は必ずしも対称とはならない。
#'
#'
#' @return 以下の要素からなるリスト
#' \preformatted{
#' theta \code{x1} に対応する特性値 theta の値
#' x1 テスト１の得点ベクトル
#' x2 テスト２の得点ベクトル
#' x2_1 テスト２に等化されたテスト１の得点ベクトル
#' p1 テスト１の確率分布
#' p2 テスト２の確率分布
#' }
#' \code{(x2_1,p1} が変換されたテスト得点の得点分布となる。
#'
#' @examples
#' res=oseq( paramS1, paramS2, print=2 )
#' res=oseq( paramS1, paramS2, weight1=weightS12, weight2=weightS21, print=2 )
#'
#' # The effect of item weight.
#' res=oseq( paramS1, paramS1, weight2=weightS12, print=3 )
#'
#' # checking if symmetric
#' # equate test1 to test2
#' res1 <- oseq( paramS1, paramS2 )
#' # equate test2 to test1
#' res2 <- oseq( paramS2, paramS1 )
#' # merge result
#' r1 <- data.frame(x1=res1$x1, y1=res1$x2_1)
#' r2 <- data.frame(y2=res2$x1, x2=res2$x2_1)
#' rx <- merge( r1, r2, by.x="x1", by.y="x2", all=TRUE)
#' rx$y <- ifelse( is.na(rx$y1), rx$y2, rx$y1 )
#' ry <- merge( r1, r2, by.x="y1", by.y="y2", all=TRUE)
#' ry$x <- ifelse( is.na(ry$x1), ry$x2, ry$x1 )
#' rxr <- round(rx,4); ryr=round(ry,4)
#' printm(rxr,ryr)
#'
#' # @keywords internal
#' @docType package
#' @name oseq_JPH
NULL
