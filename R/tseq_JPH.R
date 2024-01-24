#' 項目反応理論を用いた真値の等化 IRT True Score Equating
#'
#' この関数は、項目反応理論を用いて観測される重み付きテスト得点の
#' 条件付き期待値を算出し、それを用いてテスト得点の真値の等化を行う。\cr
#' English help file: (\link[lazy.irt]{tseq})
#'
#' @usage tseq( param1=NULL, param2=NULL, weight1=NULL, weight2=NULL
#' , method=0, interpol="spline"
#' , thmin=-5, thmax=5, npoints=121
#' , round=0, by_x=1, print=0, plot=0 )
#'
#' @param param1 テスト１ の項目パラメタデータフレーム
#' @param param2 テスト２ の項目パラメタデータフレーム
#' @param weight1 テスト１ の項目重みデータフレーム, or NULL
#' @param weight2 テスト２ の項目重みデータフレーム, or NULL
#' @param by_x テスト得点１の増分
#' @param round = 1 結果を丸める場合
#' @param method = 0 trf の逆関数を補間を用いて求める \cr
#' = 1 trf の逆関数を \code{uniroot} を用いて求める。
#' @param interpol \code{lazy.tools::interpol} の補間法のオプション
#' @param thmin trf の逆関数を求める際の特性値 theta の最小値
#' @param thmax trf の逆関数を求める際の特性値 theta の最大値
#' @param npoints 区間 \code{[thmin,thmax]} における特性値の離散点
#' @param print = 1 or 2 結果を出力
#' @param plot = 1 グラフを出力
#'
#'
#' @details
#' この関数は、項目パラメタが与えられた場合にテスト反応関数 trf を計算し
#' それらの関係を用いて、テスト１ の得点 x1 をテスト２の得点 x2
#' へ等化する。
#'
#' \code{trf1} をテスト１の、\code{trf2} をテスト２のテスト反応関数とすれば、
#' テスト得点２に等化されたテスト得点１は \code{x2_1} \cr
#' \code{x2_1 = trf2( inv_trf1(x1))}
#' で計算される。
#'
#' テスト１の得点としては、観測される整数値のみではなく
#' テスト得点の最小値と最大値の間の以下の点 \cr
#' \code{x1 <- seq(minscore,maxscore, by=by_x)} \cr
#' が用いられる。
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
#' x1 テスト１の得点ベクトル
#' x2_1 テスト２に等化されたテスト１の得点ベクトル
#' theta x1 に対応する特性値の値
#' locmin 特性値が負の極めて大きな値になる場所
#' locmax 特性値が負の極めて小さな値になる場所
#' minmax1 テスト１得点の最小値と最大値
#' minmax2 テスト２得点の最小値と最大値
#' }
#'
#'
#' @examples
#' res <- tseq( paramS1, paramS2, print=2 )
#' res <- tseq( paramS1, paramS2, weight1=weightS12, weight2=weightS21, print=2 )
#'
#' # The effect of item weight.
#' res <- tseq( paramB1, paramB1, weight1=weightB11, weight2=weightB12, print=3 )
#'
#' # checking if symmetric
#' # equate test1 to test2
#' res1 <- tseq( paramS1, paramS2, by_x=0.5 )
#' # equate test2 to test1
#' res2 <- tseq( paramS2, paramS1, by_x=0.5 )
#' # merge result
#' r1 <- data.frame(x1=res1$x1, y1=res1$x2_1)
#' r2 <- data.frame(y2=res2$x1, x2=res2$x2_1)
#' rx <- merge( r1, r2, by.x="x1", by.y="x2", all=TRUE)
#' rx$y <- ifelse( is.na(rx$y1), rx$y2, rx$y1 )
#' ry <- merge( r1, r2, by.x="y1", by.y="y2", all=TRUE)
#' ry$x <- ifelse( is.na(ry$x1), ry$x2, ry$x1 )
#' printm(rx,ry)
#'
#' # comparison of methods
#' res1 <- tseq( paramS1, paramS1, method=1, by_x=0.1, weight2=weightS12 )
#' res2 <- tseq( paramS1, paramS1, method=0, by_x=0.1, weight2=weightS12 )
#' # Print(res1$x1, res1$x2_1, res1$x2_1-res2$x2_1, fmt="6.1, 6.4 9.6")
#' summary( abs(res1$x2_1-res2$x2_1) )
#'
#'
#' # @keywords internal
#' @docType package
#' @name tseq_JPH
NULL
