#' lazy.irt: 項目反応理論のための関数群 for lazy boys and girls
#'
#'
#' 項目パラメタならびに項目重みデータフレームの説明と機能別の関数のリスト \cr
#' English help file: (\link[lazy.irt]{lazy.irt})
#'
#' @section 項目パラメタデータフレームと項目重みデータフレーム:
#' \itemize{
#'  \item 項目パラメタデータフレームの構造 \cr
#'  項目パラメタデータフレームは各項目を行に配置し、以下の変数を持つこと。
#'  \preformatted{
#'  name 項目名
#'  type 項目タイプ B | B3 | Bn | Bn3 | G | Gn | P
#'  ncat 項目カテゴリ数（２値項目の場合は 2）
#'  p1 項目パラメタ 1 （項目識別力パラメタ）
#'  p2 項目パラメタ 2 （項目困難度パラメタ）
#'  p3 項目パラメタ 3 （B3 もしくは Bn3 の場合は c パラメタ、それ以外は困難度）
#'  ::
#'  ::
#'  }
#'  \code{type = "B"  or "B2"  or  "B3"   or   "Bn"　or "Bn3"} の場合
#'  \preformatted{
#'     p1 識別力パラメタ \eqn{a_j}
#'     p2 困難度パラメタ \eqn{b_j}
#'     p3 ゲッシングパラメタ \eqn{c_j} or 0
#'     }
#' \code{type = "G" or "Gn"} の場合
#' \preformatted{
#'     p1 識別力パラメタ \eqn{a_j}
#'     p2 閾値パラメタ \eqn{b_{j1}}
#'     p3 閾値パラメタ  \eqn{b_{j2}}
#'     ::       ::
#'     p_ncat カテゴリの閾値パラメタ  \eqn{b_{j,[ncat-1]}}
#'     }
#' \code{type = "P"} の場合
#' \preformatted{
#'     p1 識別力パラメタ \eqn{a_j}
#'     p2 ステップパラメタ \qen{b_{j1}}
#'     p3 ステップパラメタ \qen{b_{j2}}
#'      ::     ::
#'     p_ncat ステップパラメタ  \eqn{b_{j,[ncat-1]}}
#'     }
#' \code{type = "PN"} の場合
#' \preformatted{
#'     p1 傾きパラメタ \eqn{a_j}
#'     p2 切片パラメタ \eqn{b_{j1}}
#'     p3 切片パラメタ \eqn{b_{j2}}
#'      ::     ::
#'     p_ncat 切片パラメタ \eqn{b_{j,[ncat-1]}}
#'     }
#' \code{type = "N"} の場合
#' \preformatted{
#'     p1 から p_(ncat[j]-1) 傾きパラメタ
#'     p_ncat[j] から p_2*(ncat[j]-1) 切片パラメタ
#'     }
#' ２値項目の場合は項目パラメタの数は（ゲッシングを含むため）3 である。
#' 順序付の多値項目 (P, G, Gn) の場合は、項目パラメタの数は \code{ncat[j]}。
#' 名義反応項目 (N) の項目パラメタの数は \code{2*(ncat[j]-1)} である。
#' \cr
#'  具体例は、以下に示すテストデータを参照のこと。
#'  \item 項目重みデータフレームの構造 \cr
#'  項目重みデータフレーム各項目を行に配置し、以下の変数を持つこと。
#'  \preformatted{
#'  name 項目名
#'  type 項目タイプ B | B3 | Bn | Bn3 | G | Gn | P
#'  ncat 項目カテゴリ数（２値項目の場合は 2）
#'  w 項目への重み（配点）
#'  v0 項目カテゴリ 0 への重み
#'  p1 項目カテゴリ 1 への重み
#'  p2 項目カテゴリ 2 への重み （多値項目の場合は以下へ続く）
#'  ::
#'  ::
#'  }
#'  以下に示すテストデータを参照のこと。
#' }
#'
#' @section 項目パラメタの推定と等化のための関数:
#' \itemize{
#'  \item \link[lazy.irt]{uIRT_JPH}:	1次元の項目パラメタの推定
#'  \cr\cr
#'  \item \link[lazy.irt]{est_theta}: 特性値 theta の推定
#'  \cr\cr
#'  \item \link[lazy.irt]{cala}:	項目パラメタの等化（最小２乗基準 = 項目パラメタ）
#'  \item \link[lazy.irt]{calr_JPH}:	項目パラメタの等化
#'  （最小２乗基準 = 項目反応関数）
#'  \cr\cr
#'  \item \link[lazy.irt]{tseq_JPH}: IRT を用いた True Score Equating
#'  \item \link[lazy.irt]{oseq_JPH}: IRT を用いた Observed Score Equating
#'  \cr
#'  \item \link[lazy.irt]{coseq_JPH}: 観測得点の等化
#'  \cr\cr
#'  \item \link[lazy.irt]{ordinal_reg_JPH}: 順序尺度の基準変数の
#'  １次元の説明変数への回帰
#'  \item \link[lazy.irt]{smn_JPH}: 正規分布の密度関数に比例する
#' 数値カテゴリを持つ多項分布の確率の計算

#'  \item \link[lazy.irt]{invtrf}: テスト反応関数の逆関数
#'  }
#'
#' @section 反応関数（特性曲線）の計算と表示のための関数:
#' \itemize{
#' \item \link[lazy.irt]{irf_JPH}: 項目反応関数 (irf) や項目カテゴリ関数 (icrf)
#'  の計算と表示
#' \cr\cr
#' \item \link[lazy.irt]{icrfB}: 3PLM (binary logistic model)
#' \item \link[lazy.irt]{icrfG}: GRM (graded response model) ロジスティク関数
#' \item icrfN: NRM (not yet available)
#' \item \link[lazy.irt]{icrfP}: GPCM (generalized partial credit model)
#' \item \link[lazy.irt]{icrfPN}: GPCM (generalized partial credit model)
#' NRM による表現
#' \item \link[lazy.irt]{icrfPN0}: GPCM (generalized partial credit model)
#' NRM による表現 2
#' }
#'
#'
#' @section 導関数の計算と表示のための関数:
#' \itemize{
#'  \item \link[lazy.irt]{dirf}: irf や icrf の特性値に関する１次導関数
#'  \cr\cr
#'  \item \link[lazy.irt]{dicrfB}: 3PLM
#'  \item \link[lazy.irt]{dicrfG}: GRM ロジスティク関数
#'  \item \link[lazy.irt]{dicrfN}: NRM (not yet available)
#'  \item \link[lazy.irt]{dicrfP}: GPCM
#'  \item \link[lazy.irt]{dicrfPN}: GPCM NRM による表現
#'  \item \link[lazy.irt]{dicrfPN0}: GPCM NRM による表現 2
#'   \cr
#'  \item \link[lazy.irt]{dicrf_num}: lazy.irt::JacobianMat を用いた数値微分
#'   \cr\cr
#'  \item \link[lazy.irt]{dirf_p}: irf の項目パラメタに関する１次導関数の計算
#' }
#'
#'
#' @section 情報関数の計算と表示のための関数:
#' \itemize{
#'  \item \link[lazy.irt]{iif_JPH}: 項目情報関数 (iif) や
#'  項目カテゴリ情報関数 (icif) の計算
#'  \cr
#'  \item \link[lazy.irt]{info_func}: 重み付き特典から特性値を推定する際の
#'  情報関数の計算
#'  \item \link[lazy.irt]{graded_info}: 段階化されたテスト得点から
#'  特性値を推定する際の情報関数の計算
#'  \item \link[lazy.irt]{flatten_SEM}: 測定の標準誤差 (SEM)
#'   を定数とするようなテスト得点の変換
#'  \item \link[lazy.irt]{flatten_SEM_theta}: 特性値の測定の標準誤差 (SEM)
#'   を定数とするような特性値 (theta-hat) に基づくテスト得点の変換
#'  \cr\cr
#'  \item \link[lazy.irt]{GOptWeight} 大局的に最適な項目カテゴリへの重みの計算
#'  }
#'
#'
#' @section テスト得点の分布の計算のための関数:
#' \itemize{
#'  \item \link[lazy.irt]{obscore_JPH}: テスト得点の分布の計算と
#'  特性値の事後分布の計算
#'  \item \link[lazy.irt]{obscore_s}: \code{obscore} の核部分
#'  \cr\cr
#'  \item \link[lazy.irt]{sumsmnw_JPH}: 数値カテゴリを持つ複数の多項分布の
#'  重み付き和の計算
#'  \item \link[lazy.irt]{sumsmnw12}: 数値カテゴリを持つ二つの多項分布の
#'  重み付き和の計算
#'  \cr\cr
#'  \item \link[lazy.irt]{rel_irt}: irt の元でのテストの信頼性係数等の計算
#'
#' }
#'
#'
#' @section 項目反応モデル間の相互変換のための関数:
#' \itemize{
#'  \item \link[lazy.irt]{fitG2P}: GPCM から GRM への変換（簡易版）
#'  \item \link[lazy.irt]{fitP2G}: GRM から GPCM への変換（簡易版）
#'  \item \link[lazy.irt]{conv2G_JPH}: logistic GRM への変換
#'  \item \link[lazy.irt]{conv2Gn_JPH}: normal GRM への変換
#'  \item \link[lazy.irt]{conv2P_JPH}: GPCM への変換
#'  }
#'
#'
#' @section テキストファイルの読み込みのための関数:
#' \itemize{
#'  \item \link[lazy.irt]{read.param}: 項目パラメタファイルの読込み
#'  項目パラメタデータフレームを作成
#'  \item \link[lazy.irt]{read.weight}: 項目重みファイルの読込み
#'  項目重みデータフレームを作成
#'  \cr
#'  \item \link[lazy.irt]{read_blg_par}: Bilog-MG の .par ファイルを読み込み
#'  項目パラメタデータフレームを作成
#'  }
#'
#'
#' @section ユーティリティ関数:
#' \itemize{
#'  \item \link[lazy.irt]{create_weight_df}: 項目パラメタデータフレームから
#'  項目重みデータフレームを作成 natural category weight
#'  \item \link[lazy.irt]{find_minmax_score}: 項目重みデータフレームから
#'  最高・最低得点を探す
#'  \item \link[lazy.irt]{checkparam2}: 項目パラメタデータフレームの
#'  整合性のチェック
#'  \item \link[lazy.irt]{gendataIRT}: 項目反応の生成
#'  \item \link[lazy.irt]{find_mode}: icrf のモードの計算
#'  \item \link[lazy.irt]{find_intersection}: icrf の交点の計算
#'  \item \link[lazy.irt]{gen_icrfnames}: vec(icrf) の行名の作成
#'  \item \link[lazy.irt]{graded_prob}; 正規分布の区間確率の計算
#'  \cr\cr
#'  \item \link[lazy.irt]{convP2N}: GPCM (step parameters) の NRM への変換
#'  \item \link[lazy.irt]{convP2PN}: GPCM (step parameters) の NRM への変換
#'  \item \link[lazy.irt]{convPN2P}: GPCM (step parameters) の NRM への変換
#'  }
#'
#'
#' @section テストデータ:
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
#' @section LRT関係の関数:
#' \itemize{
#'  \item \link[lazy.irt]{uLRT}: LRT 単調増加制約を課した項目パラメタの推定
#'  \item \link[lazy.irt]{est_rank}: 特性値の推定.
#'  \item \link[lazy.irt]{fitI2L}: LRT から IRT への変換（簡易版）
#'  \item \link[lazy.irt]{fitI2L_ls}: LRT から IRT への変換
#'  }
#'
#' # @keywords internal
#' @docType package
#' @name lazy.irt_JPH
NULL
