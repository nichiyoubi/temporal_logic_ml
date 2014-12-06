#####################################################
#
# Usage: R --vanilla --slave < linear_reqgression.r
# 	http://d.hatena.ne.jp/n_shuyo/20090709/predictive
#
#####################################################

M <- 9		# number of basis function
alpha <- 2	# hyper parameter α
beta <- 25	# hyper parameter β
Lattice <- 30	# number of graph's lattice
s<-0.1
u_i<-0.5

# training data
xlist <- seq(0, 1, length=25)
tlist <- sin(2*pi*xlist) + rnorm(length(xlist), sd=0.2)
D0 <- data.frame(x=xlist, t=tlist)

predictive <- function(D) {
    # design matrix
    gauss_base_func <- function(x)exp(-(x-u_i)^2/(2*s*s))
    phi <- function(x) sapply(x,function(x){exp(-(x-seq(0,1,length=9))^2/(2*s*s))})
    PHI <- t(phi(D$x))

    # convariance matrix & means
    S_N_inv <- alpha * diag(9) + beta * t(PHI) %*% PHI
    S_N <- solve(S_N_inv)
    m_N <- beta * S_N %*% t(PHI) %*% D$t

    # regression function
    y <- function(x)(t(phi(x)) %*% m_N)
    plot(y, xlim=c(0,1), ylim=c(-1.2, 1.2))
    par(new=T)
    plot(D, xlim=c(0,1), ylim=c(-1.2, 1.2), ylab="")

    # predictive distribution
    var_N <- function(x) {1/beta + (t(phi(x)) %*% S_N %*% phi(x))[1]}
    function(x,t) {
        mapply(function(x,t)dnorm(t,m=(t(m_N) %*% phi(x))[1], s=var_N(x), log=T), x, t)
    }
}
draw_dist <- function(p){
    x <- seq(0, 1, length=Lattice)
    t <- seq(-1.5, 1.5, length=Lattice*2)
    z <- outer(x, t, p)
    persp(x, t, z, theta=0, phi=60, shade=0.4)
}

p <- predictive(D0[sample(length(D0$x))[1:1],])
draw_dist(p);
p <- predictive(D0[sample(length(D0$x))[1:2],])
draw_dist(p);
p <- predictive(D0[sample(length(D0$x))[1:4],])
draw_dist(p);
p <- predictive(D0)
draw_dist(p);


# pdf("PRML_graph.pdf")
# dev.off()
