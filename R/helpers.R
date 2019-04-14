#' Gehan objective function
#'
#' @param beta length-p regression coefficient of AFT model
#' @param y length-n ceonsored log time-to-event variable
#' @param delta length-n censoring indicator 1=event 0=censored
#' @param matX n*p matrix of covariates
#' @param wt length-n weight for observations
#'
#' @return Gehan objective function value
#'
#' @examples
gehan.obj <- function(beta, y, delta, matX, wt = rep(1, length(y))) {
  # dimensions must agree
  if(length(y) != length(delta) | length(y) != nrow(matX) | length(y) != length(wt))
    stop("Number of samples given by y, delta, matX, and wt must agree!")
  if(length(beta) != ncol(matX))
    stop("Number of covariates given by beta and matX must agree!")
  residual <- y - as.vector(matX %*% beta)
  
  index.reorder <- order(residual)
  residual.order <- residual[index.reorder]
  delta.order <- delta[index.reorder]
  wt.order <- wt[index.reorder]
  
  s1 <- cumsum(delta.order*wt.order)
  s2 <- delta.order*rev(cumsum(rev(wt.order)))
  return(sum((s1-s2)*residual.order*wt.order))
}

#' lm estimator for multiple studies 
#'
#' @param matZ n*p_Z matrix of outcomes
#' @param matX n*p_X matrix of covariates
#' @param study length-n study indicator
#' @param wt length-n observation weight
#'
#' @return Estimated regression coefficient matrix
#'
#' @examples
alpha.fit <- function(matZ, matX, study, wt = rep(1, nrow(matZ))) {
  
  # dimensions must agree
  if(length(y) != length(delta) | length(y) != nrow(matX) | length(y) != nrow(matZ) |
     length(y) != length(study) | length(y) != length(missing))
    stop("Number of samples given by y, delta, matX, matZ, study, missing, and wt must agree!")
  if(!is.null(beta.ini) & length(beta.ini) != ncol(matX) + ncol(matZ))
    stop("Number of covariates given by beta.ini and matX must agree!")
  
  # missing must be per-study
  if(any(apply(table(study, missing) > 0, 1, sum) > 1))
    stop("Missingness must be systematic - either a study has missingness or it doesn't. ",
         "Check study and missing.")
  # dimensions must agree
  if(nrow(matZ) != nrow(matX) | nrow(matZ) != length(study) | nrow(matZ) != length(wt))
    stop("Number of samples given by matZ, matX, study, and wt must agree!")
  
  # centering per-study because intercept terms in Gehan obj get canceled out anyway
  studies <- unique(study)
  for(i.study in studies) {
    i.ind <- study == i.study
    matZ[i.ind, ] <- 
      t(t(matZ[i.ind, , drop = FALSE]) - apply(matZ[i.ind, , drop = FALSE], 2, mean))
    matX[ind, ] <- 
      t(t(matX[i.ind, , drop = FALSE]) - apply(matX[i.ind, , drop = FALSE], 2, mean))
  }
  
  # get alpha estimates (Z'WZ)^{-1}Z'WX
  return(solve(t(matX) %*% (wt * matX), t(matX) %*% (wt* matZ)))
}

#' Naive Gehan estimator for multiple studies
#'
#' @param y length-n ceonsored log time-to-event variable
#' @param delta length-n censoring indicator 1=event 0=censored
#' @param matX n*p matrix of covariates
#' @param study length-n study indicator
#' @param wt length-n weight for observations
#' @param beta.ini initial coefficient value for optimization; will estimate from lm if not
#' provided
#'
#' @return Estimated Gehan coefficients
#'
#' @examples
gehan.fit <- function(y, delta, matX,
                      study, wt = rep(1, length(y)),
                      beta.ini = NULL) {
  # dimensions must agree
  if(length(y) != length(delta) | length(y) != nrow(matX) | length(y) != length(wt))
    stop("Number of samples given by y, delta, matX, and wt must agree!")
  if(!is.null(beta.ini) & length(beta.ini) != ncol(matX))
    stop("Number of covariates given by beta.ini and matX must agree!")
  
  # if beta.ini is missing estimate from simple lm
  if(is.null(beta.ini))
    beta.ini <- lm(y ~ matX)$coef[-1]
  
  fn.obj <- function(beta) {
    obj.all <- sapply(unique(study), function(i.study) {
      i.ind <- study == i.study
      gehan.obj(beta = beta, y = y[i.ind], delta = delta[i.ind],
                matX = matX[i.ind, , drop = FALSE], wt = wt[i.ind])
    })
    sum(obj.all)
  }
  fit <- optim(beta.ini, fn.obj)
  fit$par
}


#' Combined (non-optimal) Gehan estimator for multiple studies
#'
#' @param y length-n ceonsored log time-to-event variable
#' @param delta length-n censoring indicator 1=event 0=censored
#' @param matX n*p_X matrix of always observed covariates
#' @param matZ n*p_Z matrix of systematically missing covariates
#' @param study length-n study indicator
#' @param missing length-n TRUE/FALSE missingness indicator TRUE=Z is missing FALSE=not missing
#' @param wt length-n weight for observations
#' @param beta.ini initial coefficient value for optimization; will estimate from lm if not
#' provided
#' 
#' @return Estimated Gehan coefficients
#'
#' @examples
gehan.combined.fit <- function(y, delta, matX, matZ,
                               study, missing, wt = rep(1, length(y)),
                               beta.ini = NULL) {
  # dimensions must agree
  if(length(y) != length(delta) | length(y) != nrow(matX) | length(y) != nrow(matZ) |
     length(y) != length(study) | length(y) != length(missing) | length(y) != length(wt))
    stop("Number of samples given by y, delta, matX, matZ, study, missing, and wt must agree!")
  if(!is.null(beta.ini) & length(beta.ini) != ncol(matX) + ncol(matZ))
    stop("Number of covariates given by beta.ini and matX must agree!")
  
  # missing must be per-study
  if(any(apply(table(study, missing) > 0, 1, sum) > 1))
    stop("Missingness must be systematic - either a study has missingness or it doesn't. ",
         "Check study and missing.")
  
  # fill in the systematicly missing covariates using lm coefficients estimated from full datasets
  mat.alpha <- alpha.fit(matX[!missing, , drop = FALSE],
                         matZ[!missing, , drop = FALSE],
                         study[!missing],
                         wt = wt[!missing])
  matZ[missing, ] <- matX[missing, , drop = FALSE] %*% mat.alpha
  
  # now it's just naive Gehan estimator
  gehan.fit(y = y, delta = delta, matX = cbind(matX, matZ),
            study = study, wt = wt, beta.ini = beta.ini)
}

#' Running a estimation function over B different weights and record coefficients
#'
#' @param f the estimator function
#' @param matW n*B matrix of weights
#' @param ncores number of parallel cores to run
#' @param ... additional parameters for f
#' 
#' @importFrom foreach %dopar%
#'
#' @return p*B matrix of recorded coefficients
#'
#' @examples
perturbfn <- function(f, matW, ncores = 1, ...) {
  l.beta <- list()
  if(ncores == 1) {
    for(i in 1:ncol(matW)) l.beta[[i]] <- f(..., wt = matW[, i])
    return(Reduce("cbind", l.beta))
  }
  if(ncores > 1) {
    doParallel::registerDoParallel(ncores)
    betas <- foreach::foreach(i = 1:ncol(matW),
                              .combine = cbind) %dopar% 
      {return(f(..., wt = matW[, i]))}
    doParallel::stopImplicitCluster()
    return(betas)
  }
}

#' Obtain bivariate MLE estimate for the full model coefficients in Fibrinogen estimator
#'
#' @param betas p_beta*k_avail matrix of per-study full model coefficient estimates
#' @param gammas p_gamma*k_total matrix of per-study marginal model coefficient estimates
#' @param ns length k_total per-study sample size
#' @param Sigma (p_beta + p_gamma)*(p_beta + p_gamma) matrix of combined covariance of
#' both beta and gamma
#'
#' @return Estimated full model coefficients additionally including studies with 
#' systematically missing covariates so fitting the full model not possible
#'
#' @examples
beta.mle <- function(betas, gammas, ns, Sigma) {
  # dimensions must agree
  if(ncol(gammas) != length(ns))
    stop("Number of samples given by gammas and ns must agree!")
  if(nrow(betas) + nrow(gammas) != nrow(Sigma))
    stop("Number of covariates given by betas + gammas and Sigma must agree!")
  
  p_beta <- nrow(betas)
  p_gamma <- nrow(gammas)
  
  # Inverse of covariance matrix and its blocks
  SigmaInv <- solve(Sigma)
  SigmaInv_beta <- SigmaInv[1:p_beta, 1:p_beta, drop = FALSE]
  SigmaInv_betagamma <- SigmaInv[1:p_beta, 
                                 (p_beta+1):(p_beta + p_gamma), 
                                 drop = FALSE]
  SigmaInv_gamma <- SigmaInv[(p_beta+1):(p_beta + p_gamma), 
                             (p_beta+1):(p_beta + p_gamma), 
                             drop = FALSE]
  
  # weighted mean of betas and gammas
  k_avail <- ncol(betas)
  k_total <- ncol(gammas)
  n_avail <- sum(ns[1:k_avail])
  n_total <- sum(ns)
  beta_avail <- apply(betas, 1, function(x) sum(x*ns[1:k_avail])/n_avail)
  gamma_avail <- apply(gammas[, 1:k_avail, drop = FALSE], 1, 
                       function(x) sum(x*ns[1:k_avail])/n_avail)
  gamma_total <- apply(gammas, 1, function(x) sum(x*ns)/n_total)
  
  beta_avail - as.vector(solve(n_total*SigmaInv_beta -
                                 n_avail*SigmaInv_betagamma %*%
                                 solve(SigmaInv_gamma, t(SigmaInv_betagamma)),
                               n_total*SigmaInv_betagamma*(gamma_total - gamma_avail)))
}


#' Bivariate normal likelihood function. This is to test that beta.mle is correct
#'
#' @param mu length-(p_beta + p_gamma) vector of combined true mean of beta and gamma
#' @param betas p_beta*k_avail matrix of per-study full model coefficient estimates
#' @param gammas p_gamma*k_total matrix of per-study marginal model coefficient estimates
#' @param ns length k_total per-study sample size
#' @param Sigma (p_beta + p_gamma)*(p_beta + p_gamma) matrix of combined covariance of
#' both beta and gamma
#'
#' @return Likelihood value (which can then be optimized)
#'
#' @examples
bivariate.likelihood <- function(mu, betas, gammas, ns, Sigma) {
  # dimensions must agree
  if(ncol(gammas) != length(ns))
    stop("Number of samples given by gammas and ns must agree!")
  if(nrow(betas) + nrow(gammas) != nrow(Sigma) | length(mu) != nrow(Sigma))
    stop("Number of covariates given by betas + gammas, mu, and Sigma must agree!")
  
  p_beta <- nrow(betas)
  p_gamma <- nrow(gammas)
  
  # Inverse of covariance matrix and its blocks
  SigmaInv <- solve(Sigma)
  SigmaInv_gamma <- SigmaInv[(p_beta+1):(p_beta + p_gamma), 
                             (p_beta+1):(p_beta + p_gamma), 
                             drop = FALSE]
  
  # root n centered and scaled version of betas and gammas for likelihood calculation
  k_avail <- ncol(betas)
  k_total <- ncol(gammas)
  betas.centered.scaled <- t(t(betas - mu[1:p_beta]) * sqrt(ns[1:k_avail]))
  gammas.centered.scaled <- t(t(gammas - mu[(p_beta + 1):(p_beta + p_gamma)]) * sqrt(ns))
  mus.centered.scaled <- rbind(betas.centered.scaled, 
                               gammas.centered.scaled[, 1:k_avail, drop = FALSE])
  
  sum(sapply(1:k_avail, function(k) {
    t(mus.centered.scaled[, k, drop = FALSE]) %*%
      SigmaInv %*%
      mus.centered.scaled[, k, drop = FALSE]
  })) + 
    sum(sapply((k_avail + 1):k_total, function(k) {
      t(gammas.centered.scaled[, k, drop = FALSE]) %*%
        SigmaInv_gamma %*%
        gammas.centered.scaled[, k, drop = FALSE]
    }))
}