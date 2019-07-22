#' Fit optimized Gehan estimator to multipe studies with systematically missing variables.
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
gehan.opt <- function(y, delta, matX, matZ,
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
  
  # naive estimator
  beta1 <- gehan.fit(y = y[!missing], delta = delta[!missing], 
                     matX = cbind(matX[!missing, , drop = F], matZ[!missing, , drop = F]),
                     study = study[!missing], beta.ini = beta.ini)
  # combined estimator
  beta2 <- gehan.combined.fit(y = y, delta = delta, matX = matX, matZ = matZ,
                              study = study, missing = missing, beta.ini = beta.ini)
  
  
  # perturb to estimate Sigma_{beta1, beta2}
  n <- length(y)
  p <- ncol(matX) + ncol(matZ)
  matW <- matrix(rexp(n*B), nrow=n, ncol=B)
  matBeta1Pt <- perturbfn(f = gehan.fit, matW = matW[!missing, ], ncores = ncores,
                          y = y[!missing], delta = delta[!missing], 
                          matX = cbind(matX[!missing, , drop = F], matZ[!missing, , drop = F]),
                          study = study[!missing], beta.ini = beta.ini)
  matBeta2Pt <- perturbfn(f = gehan.combined.fit, matW = matW, ncores = ncores,
                          y = y, delta = delta, matX = matX, matZ = matZ,
                          study = study, missing = missing, beta.ini = beta.ini)
  matBetaPt <- rbind(matBeta1Pt, matBeta2Pt)
  
  # Calculate the optimal weighting matrix
  matSigmaPt <- cov(t(matBetaPt))
  matSigmaA <- matSigmaPt[1:p, 1:p, drop = FALSE] +
    matSigmaPt[(p+1):(2*p), (p+1):(2*p), drop = FALSE] -
    matSigmaPt[1:p, (p+1):(2*p), drop = FALSE] -
    matSigmaPt[(p+1):(2*p), 1:p, drop = FALSE]
  matSigmaB <- matSigmaPt[(p+1):(2*p), (p+1):(2*p), drop = FALSE] - 
    matSigmaPt[(p+1):(2*p), 1:p, drop = FALSE]
  
  matOpt <- matSigmaB %*% solve(matSigmaA)
  
  # coefficient and covariance estimation
  coef <- as.vector((matOpt %*% beta1 + (diag(1, p) - matOpt) %*% beta2))
  Sigma <- matSigmaPt[(p+1):(2*p), (p+1):(2*p), drop = FALSE] - 
    matOpt %*% t(matSigmaB)
  names(coef) <- rownames(Sigma)
  
  return(list(coef = coef, Sigma = Sigma))
}
