gehan.obj <- function(beta, y, delta, matX, wt) {
  if(is.null(wt)) wt <- rep(1, length(y))
  residual <- y - as.vector(matX %*% beta)

  index.reorder <- order(residual)
  residual.order <- residual[index.reorder]
  delta.order <- delta[index.reorder]
  wt.order <- wt[index.reorder]

  s1 <- cumsum(delta.order*wt.order)
  s2 <- delta.order*rev(cumsum(rev(wt.order)))
  return(sum((s1-s2)*residual.order*wt.order))
}

alpha.fit <- function(matX, matZ, study, wt) {
  mat.coef <- sapply(1:ncol(matX), function(i.col) {
    l.alpha <- lapply(unique(study), function(i.study) {
      ind <- study == i.study
      matrix(lm(matZ[i.study, i.col, drop = F] ~ matX[i.study, i.col, drop = F],
                weights = wt[ind])$coef[-1],
             ncol = 1)
    })
    apply(Reduce("cbind", l.alpha), 2, mean)
  })
  return(t(mat.coef))
}

gehan.fit <- function(y, delta, matX, matZ,
                      study, missing,
                      wt = rep(1, length(y)),
                      beta.ini = NULL) {
  if(any(missing) & ncol(matZ) != 0) {
    # estimating alphahat and imputing
    mat.alphahat <- alpha.fit(matX[!missing, , drop = FALSE],
                              matZ[!missing, , drop = FALSE],
                              study[!missing],
                              wt = wt[!missing])
    matZ[missing, , drop = FALSE] <- matX[missing, , drop = FALSE] %*% mat.alphahat
  }

  fn.obj <- function(beta) {
    obj.all <- sapply(unique(study), function(i.study) {
      ind <- study == i.study
      gehan.obj(beta = beta, y = y[ind], delta = delat[ind],
                matX = matX[ind, , drop = FALSE], wt = wt[ind])
    })
    sum(obj.all)
  }
  fit <- optim(beta.ini, fn.obj)
  fit$par
}

perturbfn <- function(y, delta, matX, matZ, study, missing,
                      matWt, B = 500, beta.ini=NULL)
{
  p <- ncol(matX) + ncol(matZ)

  #perturb and store new beta B times
  matBetaPt <- matrix(NA, B, p) #each column is a beta sample
  for(b in 1:B) matBetaPt[b ,] <- gehan.fit(
    y = y,
    delta = delta,
    matX = matX,
    matZ = matZ,
    study = study,
    missing = missing,
    wt = matWt[,b],
    beta.ini = beta.ini
  )

  return(matBetaPt)
}
