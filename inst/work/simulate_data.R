set.seed(12345)
n <- 1000
p <- 3
qo <- 2 # no of ordinal
qn <- 2 # no of normals
X <- matrix(rnorm(p * n), ncol = p)
colnames(X) <- paste0("X", 1:ncol(X))
beta <- c(2, 0, -2)

beta0n <- c(-1, 1)
theta1 <-  c(-1, 1)
theta2 <- c(-2, 2)


# S <- rWishart(1, df = 6,Sigma = diag(qo + qn))[,,1]

S <- matrix(c(1.0,  0.7,  0.8,  0.7,
              0.7,  1.0,  0.9,  0.8,
              0.8,  0.9,  1.0,  0.9,
              0.7,  0.8,  0.9,  1.0), nrow = 4)


err <- mvtnorm::rmvnorm(n, mean = rep(0, qo + qn), sigma = cov2cor(S))
y1tilde <- X %*% beta + err[,1]
y2tilde <- X %*% beta + err[,2]
y1 <- as.numeric(cut(y1tilde, c(-Inf, theta1, Inf)))
y2 <- as.numeric(cut(y2tilde, c(-Inf, theta2, Inf)))

sigma1 <- 1
sigma2 <- 2
z1 <- beta0n[1] + X %*% beta + sigma1 * err[,3]
z2 <- beta0n[2] + X %*% beta + sigma2 * err[,4]

data_toy <- cbind.data.frame(y1, y2, z1, z2, X)
save(data_toy, file = "data/data_toy.rda")



