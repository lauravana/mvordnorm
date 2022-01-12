## Libraries
library(pbivnorm)
library(mvtnorm)
library(optimx)


################
library("mvordnorm")
data("data_toy", package = "mvordnorm")
fit <- mvordnorm("y1 + y2 + z1 + z2 ~ 0 + X1 + X2 + X3", data = data_toy,
          response_types = c("ordinal", "ordinal","gaussian", "gaussian"),
          control = mvordnorm.control(se = TRUE, solver = "CG"))

print.mvordnorm(fit)
summary.mvordnorm(fit)

###### Compare estimates
beta <- c(2, 0, -2)

beta0n <- c(-1, 1)
theta1 <-  c(-1, 1)
theta2 <- c(-2, 2)
sigma1 <- 1
sigma2 <- 2
S <- matrix(c(1.0,  0.7,  0.8,  0.7,
              0.7,  1.0,  0.9,  0.8,
              0.8,  0.9,  1.0,  0.9,
              0.7,  0.8,  0.9,  1.0), nrow = 4)

pars_true <- c(theta1[1], log(diff(theta1)),
               theta2[1], log(diff(theta2)),
               beta0n,
               beta,
               log(sigma1), log(sigma2),
               backtransf_sigmas(cov2cor(S[1:4, 1:4])))


cbind(estimate = unlist(fit$res[seq_along(pars_true)]), true = pars_true)


### One response normal
library("mvordnorm")
data("data_toy", package = "mvordnorm")
fit <- mvordnorm("y1 + y2 + z1 ~ 0 + X1 + X2 + X3", data = data_toy,
                 response_types = c("ordinal", "ordinal","gaussian"),
                 control = mvordnorm.control(se = TRUE, solver = "CG"))

print.mvordnorm(fit)
summary.mvordnorm(fit)

### One response ordinal
library("mvordnorm")
data("data_toy", package = "mvordnorm")
fit <- mvordnorm("y1 + z1 + z2 ~ 0 + X1 + X2 + X3", data = data_toy,
                 response_types = c("ordinal","gaussian", "gaussian"),
                 control = mvordnorm.control(se = TRUE, solver = "CG"))

print.mvordnorm(fit)
summary.mvordnorm(fit)

### One response ordinal + one reponse normal
library("mvordnorm")
data("data_toy", package = "mvordnorm")
fit <- mvordnorm("y1 + z1 ~ 0 + X1 + X2 + X3", data = data_toy,
                 response_types = c("ordinal","gaussian"),
                 control = mvordnorm.control(se = TRUE, solver = "CG"))

print.mvordnorm(fit)
summary.mvordnorm(fit)

### binary instead of ordinal
library("mvordnorm")
data("data_toy", package = "mvordnorm")
data_toy$y1[data_toy$y1 == 3] <- 2
fit <- mvordnorm("y1 + y2 + z1 ~ 0 + X1 + X2 + X3", data = data_toy,
                 response_types = c("ordinal", "ordinal","gaussian"),
                 control = mvordnorm.control(se = TRUE, solver = "CG"))

print.mvordnorm(fit)
summary.mvordnorm(fit)

### With MISSINGS

library("mvordnorm")
data("data_toy", package = "mvordnorm")
data_toy$y1[sample(1:nrow(data_toy), 20)] <- NA
data_toy$y2[sample(1:nrow(data_toy), 20)] <- NA
fit <- mvordnorm("y1 + y2 + z1 + z2 ~ 0 + X1 + X2 + X3", data = data_toy,
                 na.action = na.pass,
                 response_types = c("ordinal", "ordinal","gaussian", "gaussian"),
                 control = mvordnorm.control(se = TRUE, solver = "CG"))

print.mvordnorm(fit)
summary.mvordnorm(fit)

