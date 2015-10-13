#' @export
#' @param Y The (unquoted) dependent variable
#' @param D_1 The (unquoted) variable indicating first-round treatment receipt. Must be numeric and take values between 0 and 1.
#' @param D_2 The (unquoted) variable indicating second-round treatment receipt. Must be numeric or missing. Units that are in the second-round assignment must have numeric values that are 0 or 1. All others will be considered as not part of the second-round experiment.
#' @param Z_1 The (unquoted) variable indicating first-round treatment assignment. Must be numeric and take values between 0 and 1.
#' @param Z_2 The (unquoted) variable indicating second-round treatment assignment. Must be numeric or missing. Units that are in the second-round assignment must have numeric values that are 0 or 1. All others will be considered as not part of the second-round experiment.
#' @param covariates A formula of the form ~ X1 + X2, where X1 and X2 are pre-treatment covariates
#' @param pr_Z2 a numeric scalar, indicating the probability of assignment in the second round. If not, provided, will be estimated from the data.
#' @param data a dataframe. If there is missingness in Y or covariates, listwise deletion will be employed.
estimate_cace_DS <-
  function(Y, D_1, D_2, Z_1, Z_2, covariates = NULL, pr_Z2 = NULL, data){

    Y <-  eval(substitute(Y), data)
    if(!is.numeric(Y)){stop("Y must be numeric.")}
    D_1 <- eval(substitute(D_1), data)
    if(!all(D_1 %in% c(0, 1))){stop("D_1 must be a numeric variable whose values are 0 or 1.")}
    D_2 <- eval(substitute(D_2), data)
    if(!is.numeric(D_2)){stop("D_2 must be numeric.")}
    Z_1 <- eval(substitute(Z_1), data)
    if(!all(Z_1 %in% c(0, 1))){stop("Z_1 must be a numeric variable whose values are 0 or 1.")}
    Z_2 <- eval(substitute(Z_2), data)
    if(!is.numeric(D_2)){stop("D_2 must be numeric.")}

    Z_2[!Z_2 %in% c(0, 1)] <- NA

    if(!is.null(covariates)){
      X <- model.frame(covariates, data, na.action = NULL)
      not_missing <- !is.na(Y) & complete.cases(X)
      df <- data.frame(Y, D_1, D_2, Z_1, Z_2, X)
      df <- subset(df, not_missing)
    }else{
      not_missing <- !is.na(Y)
      df <- data.frame(Y, D_1, D_2, Z_1, Z_2)
      df <- subset(df, not_missing)
    }

    df <- within(df,{
      D <- (D_1==1) + (D_2==1)
    })

    # Obtain formulae

    if(is.null(covariates)){
      pr_1c_formula <- formula(D_1 ~ Z_1)
      pr_2c_formula <- formula(D_2 ~ Z_2)
      ITT1_formula <- formula(Y ~ Z_1)
      ITT2_formula <- formula(Y ~ Z_2)
    }else{
      pr_1c_formula <- formula(paste0("D_1 ~ Z_1 +", as.character(covariates)[2]))
      pr_2c_formula <- formula(paste0("D_2 ~ Z_2 +", as.character(covariates)[2]))
      ITT1_formula <- formula(paste0("Y ~ Z_1 +", as.character(covariates)[2]))
      ITT2_formula <- formula(paste0("Y ~ Z_2 +", as.character(covariates)[2]))
    }

    estimates <- DS_estimator(pr_1c_formula = pr_1c_formula,
                               pr_2c_formula = pr_2c_formula,
                               ITT1_formula = ITT1_formula,
                               ITT2_formula = ITT2_formula,
                               pr_Z2 = pr_Z2, df=df)
    return_object <- list(estimates = estimates,
                       data = df,
                       pr_1c_formula = pr_1c_formula,
                       pr_2c_formula = pr_2c_formula,
                       ITT1_formula = ITT1_formula,
                       ITT2_formula = ITT2_formula,
                       pr_Z2 = pr_Z2)
    class(return_object) <- "dsnc"
    return(return_object)
  }

#' @export
print.dsnc <- function(object, ...){
  coef_matrix <- cbind(object$estimates)
  colnames(coef_matrix) <- c("Estimate")
  rownames(coef_matrix) <- c("1CACE", "2CACE", "CACE")
  print(coef_matrix)
}

#' @param bootstrap logical, FALSE by default. Should standard errors and 95% confidence intervals be obtained via the nonparametric bootstrap?
#' @param sims If bootstrap = TRUE, the number of bootstrap simulations calculated. By default, sims = 500.
#' @export
summary.dsnc <- function(object, bootstrap=FALSE, sims = 500,...){
  if(bootstrap){
    bootstrapped_estimates <- sapply(X = 1:sims,
                                     FUN = x <- function(i){
                                       DS_estimator(pr_1c_formula = object$pr_1c_formula,
                                                     pr_2c_formula = object$pr_2c_formula,
                                                     ITT1_formula = object$ITT1_formula,
                                                     ITT2_formula = object$ITT2_formula,
                                                     pr_Z2 = object$pr_Z2, df=object$data[sample(1:nrow(object$data),
                                                                                                 nrow(object$data),
                                                                                                 replace = TRUE),])
                                     })
    ses_vec <- apply(bootstrapped_estimates, 1, sd)
    cis <- apply(bootstrapped_estimates, 1, quantile, probs=c(0.025, 0.975))

    coef_matrix <- cbind(object$estimates, ses_vec, t(cis))

    colnames(coef_matrix) <- c("Estimate", "Std. Error", "CI Lower", "CI Upper")
    rownames(coef_matrix) <- c("1CACE", "2CACE", "CACE")
  }else{
    coef_matrix <- cbind(object$estimates)
    colnames(coef_matrix) <- c("Estimate")
    rownames(coef_matrix) <- c("1CACE", "2CACE", "CACE")
  }

  pr_1c_hat <-lm(object$pr_1c_formula, data = object$data)$coefficients[2]
  pr_firstroundnoncompliers <- with(object$data, mean(D_1[Z_1==1]==0))
  pr_2c_hat <- lm(object$pr_2c_formula, data = subset(object$data, !is.na(Z_2)))$coefficients[2] * pr_firstroundnoncompliers


  return_object <- list(coef_matrix = coef_matrix,
                        pr_1c_hat = pr_1c_hat,
                        pr_2c_hat = pr_2c_hat)
  class(return_object) <- "summary.dsnc"
  return(return_object)
}

#' @export
print.summary.dsnc <- function(object, ...){

  cat("\n Double Sampling for Noncompliance: \n")

  print(object$coef_matrix)

  cat(paste0("\n Estimate of proportion of first-round compliers: ", round(object$pr_1c_hat, 3)))
  cat(paste0("\n Estimate of proportion of second-round compliers: ", round(object$pr_2c_hat, 3)))

}

#' @export
DS_estimator <- function(pr_1c_formula, pr_2c_formula, ITT1_formula, ITT2_formula,pr_Z2, df){

  # Obtain building blocks
  pr_1c_hat <-lm(pr_1c_formula, data = df)$coefficients[2]
  pr_firstroundnoncompliers <- with(df, mean(D_1[Z_1==1]==0))
  pr_2c_hat <- lm(pr_2c_formula, data = subset(df, !is.na(Z_2)))$coefficients[2] * pr_firstroundnoncompliers

  ITT1_hat <- lm(ITT1_formula, data=df)$coefficients[2]
  ITT2_hat <- lm(ITT2_formula, data=subset(df,!is.na(Z_2)))$coefficients[2]

  # Obtain Estimates
  cace2_hat <- ITT2_hat/ (pr_2c_hat/pr_firstroundnoncompliers)
  cace1_hat <- (ITT1_hat - pr_2c_hat * pr_Z2 *cace2_hat)/pr_1c_hat
  CACE_hat <- (pr_1c_hat*cace1_hat + pr_2c_hat*cace2_hat)/(pr_1c_hat + pr_2c_hat)

  return(c(cace1_hat=cace1_hat,
           cace2_hat=cace2_hat,
           CACE_hat = CACE_hat))
}

# CACE Estimate
#' @export
estimate_cace <- function(Y, D, Z, covariates = NULL, data){

  require(AER)
  Y <-  eval(substitute(Y), data)
  if(!is.numeric(Y)){stop("Y must be numeric.")}
  D <- eval(substitute(D), data)
  if(!all(D %in% c(0, 1))){stop("D must be a numeric variable whose values are 0 or 1.")}
  Z <- eval(substitute(Z), data)
  if(!all(Z %in% c(0, 1))){stop("Z must be a numeric variable whose values are 0 or 1.")}

  if(!is.null(covariates)){
    X <- model.frame(covariates, data, na.action = NULL)
    not_missing <- !is.na(Y) & complete.cases(X)
    df <- data.frame(Y, D, Z, X)
    df <- subset(df, not_missing)
  }else{
    not_missing <- !is.na(Y)
    df <- data.frame(Y, D, Z)
    df <- subset(df, not_missing)
  }

  if(is.null(covariates)){
    CACE_formula <- "Y ~ D | Z"
  }else{
    CACE_formula <- paste0("Y ~ D +",
                           as.character(covariates)[2],
                           " | Z + ",
                           as.character(covariates)[2])
  }

  fit <- ivreg(CACE_formula, data = df)

  return_object <- list(fit = fit,
                        CACE_formula = CACE_formula,
                        df = df)
  class(return_object) <- "cace"
 return(return_object)
}

#' @export
print.cace <- function(object){
  print(object$fit)
}

#' @export
summary.cace <- function(object){
  print(summary(object$fit))
}



