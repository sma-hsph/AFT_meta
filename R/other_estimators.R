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
    stop("Number of covariates given by beta.ini and matX must agree!")
  
  # missing must be per-study
  missingness <- table(study, missing) > 0
  if(any(apply(missingness, 1, sum) > 1))
    stop("Missingness must be systematic - either a study has missingness or it doesn't. ",
         "Check study and missing.")
  missingness <- c(FALSE, TRUE)[apply(missingness, 1, which)]
  
  
  # Fit full/marginal models per-study, record coefficient and covariance estimates
  studies <- unique(study)
  p_beta <- ncol(matX)
  p_gamma <- ncol(matZ)
  matBeta <- matrix(nrow = p_beta, ncol = 0)
  matGamma <- matrix(nrow = p_gamma, ncol = 0)
  ns <- c()
  lSigma <- list()
  lSigmaGamma <- list()
  # Estimations for studies with all covariates
  for (i.study in studies[!missingness]) {
    i.ind <- studies == i.study
    i.n <- sum(i.ind)
    i.beta <- gehan.fit(y = y[i.ind], delta = delta[i.ind], 
                        matX = cbind(matX, matZ)[i.ind, , drop = FALSE],
                        study = study[i.ind])
    i.gamma <- gehan.fit(y = y[i.ind], delta = delta[i.ind], 
                         matX = matX[i.ind, , drop = FALSE],
                         study = study[i.ind])
    
    # perturbing to estimate covariance matrix
    i.matW <- matrix(rexp(i.n*B), nrow=n, ncol=B)
    i.matBetaPt <- perturbfn(f = gehan.fit, matW = i.matW, ncores = ncores,
                             y = y[i.ind], delta = delta[i.ind], 
                             matX = cbind(matX, matZ)[i.ind, , drop = FALSE],
                             study = study[i.ind])
    i.matGammaPt <- perturbfn(f = gehan.fit, matW = i.matW, ncores = ncores,
                              y = y[i.ind], delta = delta[i.ind], 
                              matX = matX[i.ind, , drop = FALSE],
                              study = study[i.ind])
    i.Sigma <- cov(t(rbind(i.matBetaPt, i.matGammaPt)))
    
    matBeta <- cbind(matBeta, i.beta)
    matGamma <- cbind(matGamma, i.gamma)
    ns <- c(ns, i.n)
    lSigma <- c(lSigma, i.Sigma)
  }
  # Estimations for studies where some covariates are systematically missing
  for (i.study in studies[missingness]) {
    i.ind <- studies == i.study
    i.n <- sum(i.ind)
    i.gamma <- gehan.fit(y = y[i.ind], delta = delta[i.ind], 
                         matX = matX[i.ind, , drop = FALSE],
                         study = study[i.ind])
    
    # perturbing to estimate covariance matrix
    i.matW <- matrix(rexp(i.n*B), nrow=n, ncol=B)
    i.matGammaPt <- perturbfn(f = gehan.fit, matW = i.matW, ncores = ncores,
                              y = y[i.ind], delta = delta[i.ind], 
                              matX = matX[i.ind, , drop = FALSE],
                              study = study[i.ind])
    i.Sigma <- cov(t(i.matGammaPt))
    
    matGamma <- cbind(matGamma, i.gamma)
    ns <- c(ns, i.n)
    lSigmaGamma <- c(lSigmaGamma, i.Sigma)
  }
  
  # Estimate Sigma
  # Currently this is estimated as average of Sigmas across studies scaled by root n
  ## FIXME??
  p_beta <- ncol(matX)
  p_gamma <- ncol(matZ)
  Sigma <- Reduce("+",
                  lapply(1:length(lSigma), function(i) lSigma[[i]] * sqrt(ns[i]))) /
    length(lSigma)
  Sigma[(p_beta+1):(p_beta + p_gamma), (p_beta+1):(p_beta + p_gamma)] <- 
    (Sigma[(p_beta+1):(p_beta + p_gamma), (p_beta+1):(p_beta + p_gamma)] * length(lSigma) +
       Reduce("+",
              lapply(1:length(lSigmaGamma), function(i) lSigmaGamma[[i]] * 
                       sqrt(ns[length(lSigma) + i])))) /
    length(ns)
  
  coef <- beta.mle(betas = matBeta, gammas = matGamma, ns = ns, Sigma = Sigma)
  return(list(coef = coef))
}

beta_mi <- function(y, delta, matX, matZ,
                    study, missing, m = 5,
                    beta.ini = NULL){
  matZ_imp <- matZ
  matZ_imp[missing, ] <- NA
  pX <- ncol(matX)
  pZ <- ncol(matZ)
  
  # data frame for imputation
  df_imp <- data.frame(y, delta, matX, matZ_imp)
  # prediction matrix for multiple imputation
  predictorMatrix <- rbind(matrix(0, nrow = 1, ncol = ncol(df_imp)),
                           matrix(0, nrow = 1, ncol = ncol(df_imp)),
                           matrix(0, nrow = pX, ncol = ncol(df_imp)),
                           cbind(rep(1, pZ),
                                 rep(1, pZ),
                                 matrix(1, nrow = pZ, ncol = pX),
                                 matrix(0, nrow = pZ, ncol = pZ))
  )
  mi_fit <- mice(data = df_imp,
                 m = m,
                 meth = 'norm',
                 predictorMatrix = predictorMatrix,
                 printFlag=F)
  
  # generate beta estimate from each dataset and then average
  sapply(1:m, function(j) {
    matZ_tmp <- matZ
    matZ_tmp[missing, ] <- mi_fit$imp[(3 + pX):(2 + pX + pZ)] %>%
      sapply(function(mat) mat[, j])
    gehan.fit(
      y = y,
      delta = delta,
      matX = matX,
      matZ = matZ_tmp,
      study = study,
      missing = rep(F, length(y)),
      beta.ini = beta.ini
    )
  }) %>% apply(1, mean)
}


# some numeric testing?
Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
betahat <- c(1)
gammahat <- c(1, 2)
likelihood <- function(mu) {
  (betahat - mu[1])^2 * Sigma[1, 1] + 
    2*(betahat - mu[1])*Sigma[1, 2]*(gammahat[1] - mu[2]) +
    (gammahat[1] - mu[2])^2 * Sigma[2, 2] +
    (gammahat[2] - mu[2])^2 * Sigma[2, 2] 
}
1 / (2*Sigma[1, 1] - 1*Sigma[1, 2]/Sigma[2, 2]*Sigma[1, 2]) * 
  Sigma[1, 2]*2*(mean(gammahat) - gammahat[1])