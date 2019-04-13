#' Gehan objective function
#'
#' @param beta length-p regression coefficient of AFT model
#' @param y length-n ceonsored time-to-event variable
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
  if(nrow(matZ) != nrow(matX) | nrow(matZ) != length(study) | nrow(matZ) != length(wt))
    stop("Number of samples given by matZ, matX, study, and wt must agree!")
  
  # centering per-study because intercept terms in Gehan obj get canceled out anyway
  studies <- unique(study)
  for(i.study in studies) {
    ind <- study == i.study
    matZ[ind, ] <- 
      t(t(matZ[ind, , drop = FALSE]) - apply(matZ[ind, , drop = FALSE], 2, mean))
    matX[ind, ] <- 
      t(t(matX[ind, , drop = FALSE]) - apply(matX[ind, , drop = FALSE], 2, mean))
  }
  
  # get alpha estimates (Z'WZ)^{-1}Z'WX
  return(solve(t(matX) %*% (wt * matX), t(matX) %*% (wt* matZ)))
}

#' Naive Gehan estimator for multiple studies
#'
#' @param y length-n ceonsored time-to-event variable
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
      ind <- study == i.study
      gehan.obj(beta = beta, y = y[ind], delta = delta[ind],
                matX = matX[ind, , drop = FALSE], wt = wt[ind])
    })
    sum(obj.all)
  }
  fit <- optim(beta.ini, fn.obj)
  fit$par
}


#' Combined (non-optimal) Gehan estimator for multiple studies
#'
#' @param y length-n ceonsored time-to-event variable
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

beta_fib <- function(y, delta, matX, matZ,
                     study, missing, B = 500,
                     beta.ini = NULL, gamma.ini = NULL){
  beta1 <- sapply(unique(study[!missing]), function(iStudy) {
    gehan.fit(
      y = y[study == iStudy],
      delta = delta[study == iStudy],
      matX = matX[study == iStudy, , drop = F],
      matZ = matZ[study == iStudy, , drop = F],
      study = study[study == iStudy],
      missing = missing[study == iStudy],
      beta.ini = beta.ini
    ) * sum(study == iStudy)
  }) %>% apply(1, sum) %>% `/`(sum(!missing))
  gamma1 <- sapply(unique(study[!missing]), function(iStudy) {
    gehan.fit(
      y = y[study == iStudy],
      delta = delta[study == iStudy],
      matX = matX[study == iStudy, , drop = F],
      matZ = matrix(NA, sum(study == iStudy), 0),
      study = study[study == iStudy],
      missing = missing[study == iStudy],
      beta.ini = gamma.ini
    ) * sum(study == iStudy)
  }) %>% apply(1, sum) %>% `/`(sum(!missing))
  gamma2 <- sapply(unique(study), function(iStudy) {
    gehan.fit(
      y = y[study == iStudy],
      delta = delta[study == iStudy],
      matX = matX[study == iStudy, , drop = F],
      matZ = matrix(NA, sum(study == iStudy), 0),
      study = study[study == iStudy],
      missing = missing[study == iStudy],
      beta.ini = gamma.ini
    ) * sum(study == iStudy)
  }) %>% apply(1, sum) %>% `/`(length(y))
  
  n <- sum(!missing)
  pX <- ncol(matX)
  pZ <- ncol(matZ)
  p <- pX + pZ
  matWt <- matrix(rexp(n*B), nrow=n, ncol=B)
  matBeta1Pt <- perturbfn(
    y = y[!missing],
    delta = delta[!missing],
    matX = matX[!missing, , drop = F],
    matZ = matZ[!missing, , drop = F],
    study = study[!missing],
    missing = missing[!missing],
    matWt = matWt,
    B = B,
    beta.ini = beta.ini
  )
  matGammaPt <- perturbfn(
    y = y[!missing],
    delta = delta[!missing],
    matX = matX[!missing, , drop = F],
    matZ = matrix(NA, nrow = sum(!missing), ncol = 0),
    study = study[!missing],
    missing = missing[!missing],
    matWt = matWt,
    B = B,
    beta.ini = gamma.ini
  )
  matSigmaBetaGammaPt <- cov(cbind(matBeta1Pt, matGammaPt))
  matSigmaBetaGammaPtInv <- solve(matSigmaBetaGammaPt)
  
  (beta1 -
      solve(matSigmaBetaGammaPtInv[1:p, 1:p],
            matSigmaBetaGammaPtInv[1:p, (p + 1):(p + pX)]) %*%
      (gamma2 - gamma1)) %>% as.vector
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



