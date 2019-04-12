
simdata <- function(b, a, # regression parameters.
                    study, # vector specifying studies
                    dist_eT, sd_eT, # distribution of error terms
                    percCens, # distribution of censoring variable
                    missing
) {
  # establishing dimensions of data
  b <- matrix(b, ncol=1)
  if(is.vector(a)) a <- matrix(a, ncol=1)
  pall <- nrow(b)
  pZ <- ncol(a)
  pX <- pall - pZ
  if(pX != nrow(a)) stop('a should be pX x pZ matrix!')

  dfData <- lapply(1:length(unique(study)), function(i) {
    iStudy <- unique(study)[i]
    n_tmp <- sum(study %in% iStudy)
    matX <- mvrnorm(n_tmp, mu=rep(0, pX), Sigma=diag(rep(1, pX)))
    matZ <- matX %*% a + rnorm(n_tmp)
    eT <- switch(dist_eT[i],
                 normal = rnorm(n_tmp),
                 evd = rgev(n_tmp))
    eT <- (eT - mean(eT)) / sd(eT) * sd_eT[i] + mean(eT)
    logT <- (cbind(matX, matZ) %*% b + eT) %>% as.vector

    # generate censoring
    muCens <- 50
    sdCens <- 10
    obsCens <- 0
    while(obsCens < percCens[i]) {
      # Keep changing censoring dist'n until the censoring percentage gets close to the disired percentage.
      muCens <- muCens - 0.25
      logC <- rnorm(n_tmp, muCens, sdCens )
      delta <- 1 * (logT <= logC)
      obsCens <- 1 - sum(delta) / sum(n_tmp)
    }
    logY <- pmin(logT, logC)

    data.frame(matX, matZ,
               logY = logY, delta = delta, logT = logT,
               logC = logC, study = iStudy, missing = missing[i])
  }) %>% Reduce('rbind', .)

  colnames(dfData)[1:pX] <- paste0("X", 1:pX)
  colnames(dfData)[(pX + 1):(pX + pZ)] <- paste0("Z", 1:pZ)
  return(dfData)
}
