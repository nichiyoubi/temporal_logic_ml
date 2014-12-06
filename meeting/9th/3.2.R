# 3.2.R
# 図3.5のbias-variance分解の実装
# Start()で実行する
# Author: Takahiro Sasaki
###############################################################################

library(MASS)

# 観測試行回数
L <- 100

# ダミー（φ0=1）含む基底関数の数
M <- 25

# 1試行あたりの観測点の数
N <- 25

# 正則項
lambda <- 0.001

# 使用する基底関数
BaseFunction <- NULL

Start <- function(base.function = Gauss) {
  # メイン処理. L回試行してグラフを描画する.
  # 引数で基底関数を指定する.
  # Args:
  #   base.function: 基底関数. GaussまたはPolynomial. 指定しなければGauss
  
  BaseFunction <<- base.function

  # 真の分布関数のプロット
  x.grid <- seq(0, 1, 0.01)
  
  curve(TrueFunction, from = 0, to = 1,
        type="l", col="Red", lwd=3, xlab="x", ylab="t")
  
  # L回試行
  y.average <- rep(0, length(x.grid))
  for (k in 1:L) {
    # ノイズ付き観測データ(x, t)をN個ランダム生成
    x <- runif(N)
    t <- TrueFunction(x) + rnorm(n = N, mean = 0, sd = 0.3)
    
    # 観測データからwを導出
    w <- GetW(x, t)
    
    # 得たwを用いてyを導出
    y <- GetY(x.grid, w)
    y.average <- y.average + y
    
    # 試行20回分までをプロットする
    if (k <= 20) {
      lines(x.grid, y, col = "Darkgray")
      points(x, t)
    }
  }
  
  # yの平均をプロット
  y.average <- y.average / L
  lines(x.grid, y.average, lwd = 2, col = "Blue")
  
  # 凡例出力
  legend("topright",
         legend = c("真の分布", "予測期待値"), 
         lwd = c(3, 2),
         col = c("red", "Blue"))
}

# 基底ガウス関数 i=0,1,..,M-1
Gauss <- function(x, i) {
  if (i == 0) {
    return (1)
  }
  sigma <- 0.5
  return (exp(-(x - (i - 1) / (M - 2)) ^ 2 / (2 * sigma ^ 2)))
}

# 基底線形関数 i=0,1,..,M-1
Polynomial <- function(x, i) {
  x ^ i
}

# 真の分布関数
TrueFunction <- function(x){
  sin(2 * pi * x)
}

# wの計算
GetW <- function(x, t) {
  # φの構築
  design.matrix <- matrix(ncol = M, nrow = N)
  for (i in 1:M) {
    design.matrix[, i] <- BaseFunction(x, i - 1)
  }
  
  if (lambda == 0) {
    # λ=0のとき精度エラーが出て計算できなかったのでRの組込関数を使う
    phi <- ginv(design.matrix) 
  } else {
    phi <- solve(diag(lambda, M) + crossprod(design.matrix)) %*%
      t(design.matrix)
  }
  phi %*% t
}

# y(x, w)の計算
GetY <- function(x, w) {
  y.value <- 0
  for (i in 1:M) {
    y.value <- y.value + w[i] * BaseFunction(x, i - 1)
  }
  y.value
}

