#' 順序的カテゴリを持つ基準変数の１次元の説明変数への回帰
#'
#' 度数分布 \code{(R,x)} が与えたれたときに \code{R} を \code{x} へ回帰する \cr
#' English help file: (\link[lazy.irt]{ordinal_reg})
#'
#' @usage ordinal_reg( R, x, type="P", param=NULL
#'  , maxiter=100, eps=1e-8, epsd=1e-5, smallP=0, DinP=1
#'  , minp1=1e-2, maxabsparam=20, print=1 )
#'
#' @param R lenght(x) x カテゴリ数の基準変数、ただし \cr
#'   \code{R[,1]} は第 0 番目のカテゴリの度数のベクトル、
#'   \code{R[,2]} は第 1 番目のカテゴリの度数ベクトル, ...
#'   \code{R[,ncat]} は第 \code{ncat-1} 番目のカテゴリの度数ベクトル
#' @param x １次元の連続量の説明変数ベクトル
#' @param type 項目タイプで  "B", "B3", "G", "P", "PN", "N"
#' @param param 項目パラメタデータフレームの初期値
#' @param maxiter 繰り返し数の上限
#' @param eps 対数尤度の増加率を用いた収束判定基準
#' @param epsd パラメタの変化の絶対値差を用いた収束判定基準
#' @param smallP 確率の最小値
#' @param DinP = 0 ロジスティク関数で D=1.7 を用いない。
#' @param minp1 type != "N" の場合の p1 の最小値
#' @param maxabsparam パラメタの最大値
#' @param print = 0 結果を出力しない。
#'
#' @details
#' \code{R,x} を説明変数ベクトルの各値ごとの基準変数の各カテゴリの度数分布と
#' する。ただし、\code{R} は \code{n x ncat} の行列である。  \cr
#' この関数は以下の対数尤度関数を最大とするようなパラメタ \code{p1,p2.,,,} を
#' 求める。 \cr
#' \eqn{ llh = \sum_{i=1}^{n} \sum_{k=1}^{ncat} R[i,k] log( P_k(x[i] | param)}
#' \cr
#' ただし、\eqn{P_k(x[i] | param)} は \code{type} の指定で決まる項目反応関数
#' であり、 \eqn{param} はそのパラメタ  \code{p1,p2.,,,} である。
#'
#' たとえば \code{type="B"} の場合はロジスティク回帰に、
#' \code{type="P"} の場合は多項ロジット回帰に相当する。
#'
#' また、 \code{type="B3"} の場合には \eqn{p3=c} パラメタをロジット変化した
#' \eqn{y=logit(c)} がパラメタとして用いられるが、その場合 \cr
#' \eqn{d \; llh / d; y = d \; llh / d \; c \times c(1-c)}, \cr
#' また \eqn{d \; logistic(y) / d \; y = c(1-c)} である。
#'
#'
#' @return 以下の要素からなるリスト
#' \preformatted{
#' param 推定された項目パラメタデータフレーム
#' llh 最大化された対数尤度
#' P \code{x] における確率}
#' x 説明変数ベクトル
#' R 基準変数
#' g 対数尤度の１次微分
#' maxag １次微分の絶対値の最大値
#' }
#'
#' @examples
#' # binary data
#' set.seed(1701)
#' theta=seq(-4,4,length.out=51)
#' resg=gendataIRT( 500, paramS1[1,], theta=theta, thdist="NORMAL", thd=NULL
#'                , thmean=0, thstd=1, compress=0 )
#' testdata=as.matrix( resg$U )
#'
#' # 3PLM
#' res=ordinal_reg( testdata, theta, type="B3", param=NULL, print=1 )
#' # 2PLM
#' res=ordinal_reg( testdata, theta, type="B", param=NULL, print=1 )
#' res=ordinal_reg( testdata, theta, type="P", param=NULL, print=1 )
#' # 2PLM with nonzero c parameter
#' initp=paramS1[2,]; initp$type="B"
#' res=ordinal_reg( testdata, theta, type="B", param=initp, print=1 )
#'
#'
#' # polytomous data
#' set.seed(1701)
#' theta=seq(-4,4,length.out=51)
#' resg=gendataIRT( 500, paramS1[3,], theta=theta, thdist="NORMAL", thd=NULL
#'                  , thmean=0, thstd=1, compress=0 )
#' testdata=as.matrix( resg$U )
#'
#' # graded response model
#' res=ordinal_reg( testdata, theta, type="G", param=NULL, print=1 )
#' # partial credit model
#' res=ordinal_reg( testdata, theta, type="P", param=NULL, print=1 )
#' res=ordinal_reg( testdata, theta, type="PN", param=NULL, print=1 )
#'
#'
#' # @keywords internal
#' @docType package
#' @name ordinal_reg_JPH
NULL
