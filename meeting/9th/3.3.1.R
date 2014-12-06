# 3.3.1.R
# 図3.7のy(x,w)=w0 + w1x のベイズ学習
# Start()で実行する
# Author: Takahiro Sasaki
###############################################################################

kNumberOfLearning <- 10
alpha <- 2.0
beta <- (1 / 0.2) ^ 2

TrueFunction <- function(x) {
  # 真の関数 Y = a0 + a1 X
  -0.3 + 0.5 * x
}

# 処理実行
Start <- function() {

  # 初期事前分布
  m <- c(0, 0)
  s <- diag(1 / alpha, 2)

  par(mfrow = c(1, 2))
  Draw(m, s)
  curve(TrueFunction, from = -1, to = 1, 
        col = "red", lwd = 2,
        main = "data", xlab = "x", ylab = "t")
  
  # 訓練データ(x, t)を1対ずつ流して事後分布を描画する
  x <- runif(kNumberOfLearning, min = -1, max = 1)
  t <- TrueFunction(x) + rnorm(length(x), mean = 0, sd = 0.2)
  for (i in 1:length(x)) {
    cat("Hit Return", "\n")
    readline()

    # 新しいデータでmとsを更新する
    cat("x=", x[i], ", t=", t[i], "\n")
    posterior <- GetPosterior(m, s, x[i], t[i])
    m <- posterior$m
    s <- posterior$s
    Draw(m, s)
    curve(TrueFunction, from = -1, to = 1, 
          col = "red", lwd = 2,
          main = "data", xlab = "x", ylab = "t")
    points(x = x[1:i], y = t[1:i])
    abline(m)
  }
}

p.w <- function(w, m, s) {
  # 2次元ガウス分布 p(w) = N(w | m, s)
  1 / (2 * pi * sqrt(det(s))) *
    exp(-1 / 2 * t(w - m) %*% solve(s) %*% (w - m))
}

GetPosterior <- function(prior.m, prior.s, x, t) {
  # 事前分布と新しい測定値(x, t)から事後分布を計算する
  # Args:
  #   prior.m: 事前分布の平均m
  #   prior.s: 事前分布の共分散s
  #   x: 入力変数x
  #   t: 目標変数t
  # Returns:
  #   事後分布の平均mと共分散sが格納されたリスト
  design.matrix <- matrix(c(1, x), nrow = 1)
  s.inverse <- solve(prior.s) + beta * t(design.matrix) %*% design.matrix
  s <- solve(s.inverse)
  m <- s %*% (solve(prior.s) %*% prior.m + beta * t(design.matrix) * t)
  list(m = m, s = s)
}

Draw <- function(m, s) {
  # 平均m、共分散sの2次元ガウス分布を描画する
  grid <- seq(-1, 1, 0.05)
  w0 <- grid
  w1 <- grid
  z <- outer(w0, w1,
             function(w0, w1) {apply(cbind(w0, w1), 1, p.w, m, s)})
  image(main = "w", w0, w1, z, col = rainbow(100))
  cat("m=", m, ", |s|=", det(s), "\n")
}
