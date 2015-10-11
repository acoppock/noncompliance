library(noncompliance)
library(testthat)


load("~/Desktop/data_for_DSNC.rdata")

df$X <- runif(n = 10000)

test <- estimate_cace_DS(Y, D_1, D_2, Z_1, Z_2, covariates = ~ X, pr_Z2 = .5, data = df)
test

test_1 <- summary(test,bootstrap = TRUE, sims = 50)
print(test_1)


cace_test <- estimate_cace(Y=Y, D=D_1, Z=Z_1, data = df)
cace_test <- estimate_cace(Y=Y, D=D_1, Z=Z_1, covariates = ~X, data = df)
summary(cace_test)

debugonce(estimate_cace)

library(randomizr)
N <- 10000

type <- sample(c("NT", "C"), N, replace = TRUE, prob = c(0.3, 0.7))
Z <- complete_ra(N, condition_names = c(0, 1, 2))
Y <- rep(NA, N)
D <- rep(NA, N)

Y[Z %in% c(0, 2) | Z %in% c(1) & type=="NT"] <-
  rbinom(n = sum(Z %in% c(0, 2) | Z %in% c(1) & type=="NT"), size = 1, prob = 0.4)

Y[Z %in% c(1) & type=="C"] <-
  rbinom(n = sum(Z %in% c(1) & type=="C"), size = 1, prob = 0.5)

D[Z == 0 | Z %in% c(1, 2) & type=="NT"] <- 0
D[Z %in% c(1, 2) & type=="C"] <- 1

df <- data.frame(Y, D, Z)

estimate_3g(Y, D, Z, df)
