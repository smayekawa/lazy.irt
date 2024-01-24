x <- 1:10
y <- c(2, 1, 3, 2, 4, 8, 14, 20, 22, 29)
z <- c(2, 1, 3, 5, 7, 7, 6, 5, 4, 6)

# 二つ目の y 軸を描くために余白を調整
# par(oma = c(0, 0, 0, 2))

# x と y のプロット
plot(x, y, xlim = c(0, 10), ylim = c(0, 30),
     xlab = "x", ylab = "y", type = "l", col = "orange", lwd = 2,  # いろいろなオプション
     axes = FALSE)                                                 # 座標軸を描かないようにする
axis(1)   # 下の横軸を表示
axis(2)   # 左の縦軸を表示

# 次のグラフを重ね合わせるための作業
par(new = TRUE)

# x と z のプロット
plot(x, z, xlim = c(0, 10), ylim = c(0, 10),
     xlab = "", ylab = "", type = "l", col = "cyan", lwd = 2,
     axes = FALSE)
mtext("z",side = 4, line = 3)  # 右の縦軸のラベル
axis(4)                        # 右の縦軸を表示

# プロット全体に枠線で囲む
box()

# 凡例
legend("topleft", legend = c("y", "z"), col = c("orange", "cyan"), lty = 1)
