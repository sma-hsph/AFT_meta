#' Fit naive (i.e., only with full observations) Gehan estimator to multipe studies with 
#' systematically missing variables.
#'
#' @param y length-n ceonsored log time-to-event variable
#' @param delta length-n censoring indicator 1=event 0=censored
#' @param matX n*p_X matrix of always observed covariates
#' @param matZ n*p_Z matrix of systematically missing covariates
#' @param study length-n study indicator
#' @param missing length-n TRUE/FALSE missingness indicator TRUE=Z is missing FALSE=not missing
#' @param beta.ini initial coefficient value for optimization; will estimate from lm if not
#' provided
#' @param B number of perturbations to perform
#' @param ncores number of parallel cores to run
#'
#' @return A list object with component coef for estimated coefficients and Sigma for estimated 
#' covariance matrix
#' @export
#'
#' @examples
gehan.obs <- function(y, delta, matX, matZ,
                      study, missing,
                      beta.ini = NULL,
                      B = 30,
                      ncores = 1, indirect = FALSE) {
  # dimensions must agree
  if(length(y) != length(delta) | length(y) != nrow(matX) | length(y) != nrow(matZ) |
     length(y) != length(study) | length(y) != length(missing))
    stop("Number of samples given by y, delta, matX, matZ, study, missing, and wt must agree!")
  if(!is.null(beta.ini) & length(beta.ini) != ncol(matX) + ncol(matZ))
    stop("Number of covariates given by beta.ini and matX + matZ must agree!")
  
  # missing must be per-study
  if(any(apply(table(study, missing) > 0, 1, sum) > 1))
    stop("Missingness must be systematic - either a study has missingness or it doesn't. ",
         "Check study and missing.")
  
  # naive estimator
  coef <- gehan.fit(y = y[!missing], delta = delta[!missing], 
                    matX = cbind(matX[!missing, , drop = F], matZ[!missing, , drop = F]),
                    study = study[!missing], beta.ini = beta.ini)
  
  # perturb to estimate Sigma
  n <- length(y)
  p <- ncol(matX) + ncol(matZ)
  matW <- matrix(rexp(n*B), nrow=n, ncol=B)
  matBetaPt <- perturbfn(f = gehan.fit, matW = matW[!missing, ], ncores = ncores,
                         y = y[!missing], delta = delta[!missing], 
                         matX = cbind(matX[!missing, , drop = F], matZ[!missing, , drop = F]),
                         study = study[!missing], beta.ini = beta.ini)
  Sigma <- cov(t(matBetaPt))

  fit <- list(coef = coef, Sigma = Sigma)

  if(indirect) {
    if(ncol(matZ) != 1)
      stop("Estimation of total indirect effect is not supported for",
           " more than one mediators!")
    matAlphaPt <- perturbfn(f = alpha.fit, matW = matW[!missing, ], 
                            ncores = ncores,
                            matX = matX[!missing, , drop = F], 
                            matZ = matZ[!missing, , drop = F],
                            study = study[!missing])
    
    coef_alpha <- apply(matAlphaPt, 1, mean)
    Sigma_alpha <- cov(t(matAlphaPt))
    
    coef_indirect <- coef_alpha * coef[(ncol(matX) + 1):p]
    Sigma_indirect <- cov(t(vapply(seq_len(B),
                                   function(b) 
                                     matAlphaPt[, b] * matBetaPt[(ncol(matX) + 1):p, b],
                                   rep(0.0, length(coef_alpha)))))
    fit <- 
      c(fit, 
        list(coef_alpha = coef_alpha, Sigma_alpha = Sigma_alpha,
             coef_indirect = coef_indirect, Sigma_indirect = Sigma_indirect))
  }
  
  return(fit)
}

#' Fit combined Gehan estimator to multipe studies with systematically missing variables.
#'
#' @param y length-n ceonsored log time-to-event variable
#' @param delta length-n censoring indicator 1=event 0=censored
#' @param matX n*p_X matrix of always observed covariates
#' @param matZ n*p_Z matrix of systematically missing covariates
#' @param study length-n study indicator
#' @param missing length-n TRUE/FALSE missingness indicator TRUE=Z is missing FALSE=not missing
#' @param beta.ini initial coefficient value for optimization; will estimate from lm if not
#' provided
#' @param B number of perturbations to perform
#' @param ncores number of parallel cores to run
#'
#' @return A list object with component coef for estimated coefficients and Sigma for estimated 
#' covariance matrix
#' @export
#'
#' @examples
gehan.full <- function(y, delta, matX, matZ,
                       study, missing,
                       beta.ini = NULL,
                       B = 30,
                       ncores = 1) {
  # dimensions must agree
  if(length(y) != length(delta) | length(y) != nrow(matX) | length(y) != nrow(matZ) |
     length(y) != length(study) | length(y) != length(missing))
    stop("Number of samples given by y, delta, matX, matZ, study, missing, and wt must agree!")
  if(!is.null(beta.ini) & length(beta.ini) != ncol(matX) + ncol(matZ))
    stop("Number of covariates given by beta.ini and matX + matZ must agree!")
  
  # missing must be per-study
  if(any(apply(table(study, missing) > 0, 1, sum) > 1))
    stop("Missingness must be systematic - either a study has missingness or it doesn't. ",
         "Check study and missing.")
  
  # combined estimator
  coef <- gehan.combined.fit(y = y, delta = delta, matX = matX, matZ = matZ,
                              study = study, missing = missing, beta.ini = beta.ini)
  
  
  # perturb to estimate Sigma
  n <- length(y)
  p <- ncol(matX) + ncol(matZ)
  matW <- matrix(rexp(n*B), nrow=n, ncol=B)
  matBetaPt <- perturbfn(f = gehan.combined.fit, matW = matW, ncores = ncores,
                          y = y, delta = delta, matX = matX, matZ = matZ,
                          study = study, missing = missing, beta.ini = beta.ini)
  Sigma <- cov(t(matBetaPt))
  
  return(list(coef = coef, Sigma = Sigma))
} 

#' Fit AFT version of Fibrinogen{19222087} estimator to multipe studies with 
#' systematically missing variables.
#'
#' @param y length-n ceonsored log time-to-event variable
#' @param delta length-n censoring indicator 1=event 0=censored
#' @param matX n*p_X matrix of always observed covariates
#' @param matZ n*p_Z matrix of systematically missing covariates
#' @param study length-n study indicator
#' @param missing length-n TRUE/FALSE missingness indicator TRUE=Z is missing FALSE=not missing
#' @param beta.ini initial full model coefficient value for optimization; will estimate 
#' from lm if not provided
#' @param gamma.ini initial marginal model coefficient value for optimization; will 
#' estimate from lm if not provided
#' @param B number of perturbations to perform
#' @param ncores number of parallel cores to run
#'
#' @return A list object with component coef for estimated coefficients and Sigma for estimated 
#' covariance matrix
#' @export
#'
#' @examples
gehan.fib <- function(y, delta, matX, matZ,
                      study, missing,
                      beta.ini = NULL, gamma.ini = NULL,
                      B = 30,
                      ncores = 1){
  # dimensions must agree
  if(length(y) != length(delta) | length(y) != nrow(matX) | length(y) != nrow(matZ) |
     length(y) != length(study) | length(y) != length(missing))
    stop("Number of samples given by y, delta, matX, matZ, study, missing, and wt must agree!")
  if(!is.null(beta.ini) & length(beta.ini) != ncol(matX) + ncol(matZ))
    stop("Number of covariates given by beta.ini and matX + matZ must agree!")
  if(!is.null(gamma.ini) & length(gamma.ini) != ncol(matX))
    stop("Number of covariates given by gamma.ini and matX must agree!")
  
  # missing must be per-study
  missingness <- table(study, missing) > 0
  if(any(apply(missingness, 1, sum) > 1))
    stop("Missingness must be systematic - either a study has missingness or it doesn't. ",
         "Check study and missing.")
  missingness <- c(FALSE, TRUE)[apply(missingness, 1, which)]
  
  
  # Fit full/marginal models per-study, record coefficient and covariance estimates
  studies <- unique(study)
  p_beta <- ncol(matX) + ncol(matZ)
  p_gamma <- ncol(matX)
  matBeta <- matrix(nrow = p_beta, ncol = 0)
  matGamma <- matrix(nrow = p_gamma, ncol = 0)
  ns <- c()
  lSigma <- list()
  lSigmaGamma <- list()
  # Estimations for studies with all covariates
  for (i.study in studies[!missingness]) {
    i.ind <- study == i.study
    i.n <- sum(i.ind)
    i.beta <- gehan.fit(y = y[i.ind], delta = delta[i.ind], 
                        matX = cbind(matX, matZ)[i.ind, , drop = FALSE],
                        study = study[i.ind], beta.ini = beta.ini)
    i.gamma <- gehan.fit(y = y[i.ind], delta = delta[i.ind], 
                         matX = matX[i.ind, , drop = FALSE],
                         study = study[i.ind], beta.ini = gamma.ini)
    
    # perturbing to estimate covariance matrix
    i.matW <- matrix(rexp(i.n*B), nrow=i.n, ncol=B)
    i.matBetaPt <- perturbfn(f = gehan.fit, matW = i.matW, ncores = ncores,
                             y = y[i.ind], delta = delta[i.ind], 
                             matX = cbind(matX, matZ)[i.ind, , drop = FALSE],
                             study = study[i.ind], beta.ini = beta.ini)
    i.matGammaPt <- perturbfn(f = gehan.fit, matW = i.matW, ncores = ncores,
                              y = y[i.ind], delta = delta[i.ind], 
                              matX = matX[i.ind, , drop = FALSE],
                              study = study[i.ind], beta.ini = gamma.ini)
    i.Sigma <- cov(t(rbind(i.matBetaPt, i.matGammaPt)))
    
    matBeta <- cbind(matBeta, i.beta)
    matGamma <- cbind(matGamma, i.gamma)
    ns <- c(ns, i.n)
    lSigma <- c(lSigma, list(i.Sigma))
  }
  # Estimations for studies where some covariates are systematically missing
  for (i.study in studies[missingness]) {
    i.ind <- study == i.study
    i.n <- sum(i.ind)
    i.gamma <- gehan.fit(y = y[i.ind], delta = delta[i.ind], 
                         matX = matX[i.ind, , drop = FALSE],
                         study = study[i.ind], beta.ini = gamma.ini)
    
    # perturbing to estimate covariance matrix
    i.matW <- matrix(rexp(i.n*B), nrow=i.n, ncol=B)
    i.matGammaPt <- perturbfn(f = gehan.fit, matW = i.matW, ncores = ncores,
                              y = y[i.ind], delta = delta[i.ind], 
                              matX = matX[i.ind, , drop = FALSE],
                              study = study[i.ind], beta.ini = gamma.ini)
    i.Sigma <- cov(t(i.matGammaPt))
    
    matGamma <- cbind(matGamma, i.gamma)
    ns <- c(ns, i.n)
    lSigmaGamma <- c(lSigmaGamma, list(i.Sigma))
  }
  
  # Estimate Sigma
  # Currently this is estimated as average of Sigmas across studies scaled by n
  ## FIXME
  Sigma <- 
    Reduce("+",
           lapply(1:length(lSigma), 
                  function(i) lSigma[[i]] * (ns[i] - 1))) / length(lSigma) / 
    (sum(ns[1:length(lSigma)]) - 1)
  Sigma[(p_beta+1):(p_beta + p_gamma), (p_beta+1):(p_beta + p_gamma)] <- 
    (Sigma[(p_beta+1):(p_beta + p_gamma), (p_beta+1):(p_beta + p_gamma)] * 
       (sum(ns[1:length(lSigma)]) - 1) +
       Reduce("+",
              lapply(1:length(lSigmaGamma), function(i) lSigmaGamma[[i]] * 
                       (ns[length(lSigma) + i] - 1))) / length(lSigmaGamma)) / 
    2 / (sum(ns) - 1)
  
  return(beta.mle(betas = matBeta, gammas = matGamma, ns = ns, Sigma = Sigma))
}

#' Fit AFT version of multiple imputation{23857554} estimator to multipe studies with 
#' systematically missing variables.
#'
#' @param y length-n ceonsored log time-to-event variable
#' @param delta length-n censoring indicator 1=event 0=censored
#' @param matX n*p_X matrix of always observed covariates
#' @param matZ n*p_Z matrix of systematically missing covariates
#' @param study length-n study indicator
#' @param missing length-n TRUE/FALSE missingness indicator TRUE=Z is missing FALSE=not missing
#' @param beta.ini initial coefficient value for optimization; will estimate from lm if not
#' provided
#' @param m number of multiple imputations
#'
#' @return A list object with component coef for estimated coefficients and Sigma for estimated 
#' covariance matrix. Currently covariance is not provided
#' @export
#'
#' @examples
gehan.mi <- function(y, delta, 
                     matX, matZ,
                     study, missing,
                     surv_est = "none",
                     surv_use = "H",
                     beta.ini = NULL,
                     m = 10, B = 30, ncores = 1){
  matZ_imp <- matZ
  matZ_imp[missing, ] <- NA
  pX <- ncol(matX)
  pZ <- ncol(matZ)
  if(pZ != 1)
    stop("Currently only one systematically missing covariate is supported!")
  
  # data frame for imputation
  df_imp <- data.frame(study, y, delta, matX, matZ_imp)
  
  if(surv_est != "none") {
    surv_val <- rep(NA_real_, length(y))
    if(surv_est == "marg") { 
      for(i.study in unique(study)) {
        i.ind <- study == i.study
        i.fit.surv <- survival::survfit(survival::Surv(exp(y[i.ind]),
                                                       delta[i.ind]) ~ 1)
        surv_val[i.ind] <- i.fit.surv$cumhaz[match(rank(y[i.ind]), 
                                                   rank(i.fit.surv$time))]
      }
    }
    if(surv_est == "cond") { 
      coef_marginal <- gehan.fit(y = y,
                                 delta = delta,
                                 matX = matX,
                                 study = study,
                                 beta.ini = beta.ini)
      epsilon <- as.vector(y - matX %*% coef_marginal)
      for(i.study in unique(study)) {
        i.ind <- study == i.study
        i.fit.surv <- survival::survfit(survival::Surv(exp(epsilon[i.ind]),
                                                     delta[i.ind]) ~ 1,
                                      type="kaplan-meier")
        surv_val[i.ind] <- i.fit.surv$cumhaz[match(rank(epsilon[i.ind]), 
                                                   rank(i.fit.surv$time))]
      }
    }
    
    if(surv_use == "S")
      surv_val <- exp(-surv_val)
    
    df_imp$surv_val <- surv_val
  }
  
  predictorMatrix <- matrix(NA, nrow = ncol(df_imp), ncol = ncol(df_imp))
  rownames(predictorMatrix) <- colnames(predictorMatrix) <- colnames(df_imp)
  method <- rep("", ncol(df_imp))
  names(method) <- colnames(df_imp)
  if(length(unique(study[!missing])) > 1) {
    # prediction matrix for multiple imputation
    predictorMatrix[,] <- 2
    diag(predictorMatrix) <- 0
    predictorMatrix[-1, 1] <- -2
    method[(4 + pX):(3 + pX + pZ)] <- "2l.2stage.norm"
  } else {
    predictorMatrix[,] <- 1
    diag(predictorMatrix) <- 0
    predictorMatrix[, 1] <- 0
    predictorMatrix[1, ] <- 0
    method[(4 + pX):(3 + pX + pZ)] <- "norm"
  }
  
  mi_fit <- mice::mice(data = df_imp, 
                       predictorMatrix = predictorMatrix,
                       method = method,
                       m = m,
                       printFlag = FALSE)
  
  # generate beta estimate from each imputation and then average
  matCoef <- sapply(1:m, function(j) {
    matZ_tmp <- matZ
    matZ_tmp[missing, ] <- sapply(mi_fit$imp[(4 + pX):(3 + pX + pZ)], 
                                  function(mat) mat[, j])
    gehan.fit(
      y = y,
      delta = delta,
      matX = cbind(matX, matZ_tmp),
      study = study,
      beta.ini = beta.ini
    )
  })
  coef <- apply(matCoef, 1, mean)
  
  return(list(coef = coef,
              Sigma = NULL ## FIXME It's not clear to me how to compute Sigma
              # Given that perturbation can't be incorporated into multiple 
              # imputation
              ))
}