# Whole bunch of old functions used for fitting, simulation, and visualization
# matBeta1 <- sapply(1:200, function(i) {
#   df_sim <- simdata(b = b,
#                     a = a,
#                     study = study,
#                     dist_eT = dist_eT,
#                     sd_eT = sd_eT,
#                     percCens = percCens,
#                     missing = missing)
#   gehan.fit(y = df_sim$logY[!missing],
#             delta = df_sim$delta[!missing],
#             matX = (df_sim[, 1:2] %>% as.matrix)[!missing, ],
#             matZ = (df_sim[, 3:4] %>% as.matrix)[!missing, ],
#             study = df_sim$study[!missing],
#             missing = df_sim$missing[!missing],
#             beta.ini = b)
# })
#
# matBeta2 <- sapply(1:200, function(i) {
#   df_sim <- simdata(b = b,
#                     a = a,
#                     study = study,
#                     dist_eT = dist_eT,
#                     sd_eT = sd_eT,
#                     percCens = percCens,
#                     missing = missing)
#   gehan.fit(y = df_sim$logY,
#             delta = df_sim$delta,
#             matX = df_sim[, 1:2] %>% as.matrix,
#             matZ = df_sim[, 3:4] %>% as.matrix,
#             study = df_sim$study,
#             missing = df_sim$missing,
#             beta.ini = b)
# })
#
# matBetaOpt <- sapply(1:20, function(i) {
#   df_sim <- simdata(b = b,
#                     a = a,
#                     study = study,
#                     dist_eT = dist_eT,
#                     sd_eT = sd_eT,
#                     percCens = percCens,
#                     missing = missing)
#   beta_opt(
#     y = df_sim$logY,
#     delta = df_sim$delta,
#     matX = df_sim[, 1:2] %>% as.matrix,
#     matZ = df_sim[, 3:4] %>% as.matrix,
#     study = df_sim$study,
#     missing = df_sim$missing,
#     beta.ini = b,
#     B = 20
#   )
# })
#
# matBetaFib <- sapply(1:20, function(i) {
#   df_sim <- simdata(b = b,
#                     a = a,
#                     study = study,
#                     dist_eT = dist_eT,
#                     sd_eT = sd_eT,
#                     percCens = percCens,
#                     missing = missing)
#   beta_fib(
#     y = df_sim$logY,
#     delta = df_sim$delta,
#     matX = df_sim[, 1:2] %>% as.matrix,
#     matZ = df_sim[, 3:4] %>% as.matrix,
#     study = df_sim$study,
#     missing = df_sim$missing,
#     beta.ini = b,
#     gamma.ini = (b[1:2] + a %*% b[3:4]) %>% as.vector,
#     B = 20
#   )
# })
#
# matBetaMI <- sapply(1:20, function(i) {
#   df_sim <- simdata(b = b,
#                     a = a,
#                     study = study,
#                     dist_eT = dist_eT,
#                     sd_eT = sd_eT,
#                     percCens = percCens,
#                     missing = missing)
#   beta_mi(
#     y = df_sim$logY,
#     delta = df_sim$delta,
#     matX = df_sim[, 1:2] %>% as.matrix,
#     matZ = df_sim[, 3:4] %>% as.matrix,
#     study = df_sim$study,
#     missing = df_sim$missing,
#     beta.ini = b,
#     m = 5
#   )
# })
#
#
# df_sim <- simdata(b = b,
#                   a = a,
#                   study = study,
#                   dist_eT = dist_eT,
#                   sd_eT = sd_eT,
#                   percCens = percCens,
#                   missing = missing)
# gehan.fit(y = df_sim$logY,
#           delta = df_sim$delta,
#           matX = df_sim[, 1:2] %>% as.matrix,
#           matZ = df_sim[, 3:4] %>% as.matrix,
#           study = df_sim$study,
#           missing = df_sim$missing,
#           beta.ini = b)
#
# coefSummary <- function(est_name, l_results, trueValues, var_est=T) {
#   mat_coef <- l_results[[paste0('mat', est_name)]] # Discarding the first replicate
#                                                                # because this round has no sd estimates for beta_star
#   p_all <- length(trueValues)
#   p <- ncol( mat_coef )
#   # if ( p_all != p ) {
#   #   # stop( 'True value should be of the same length as the coefficients!' )
#   #   trueValues <- trueValues[1:p]
#   # }
#
#   matResults <- matrix(NA,nrow=p_all,ncol=3)
#   matResults[1:p,1] <- apply( mat_coef, 2, mean )
#   matResults[1:p,2] <- sqrt( apply( t( t( mat_coef )  - trueValues[1:p] )^2, 2, mean ) )
#   matResults[1:p,3] <- apply( mat_coef, 2, sd )
#
#   # if ( var_est & (paste0( 'lMatSigma', est_name, 'Pt' ) %in% names( l_results )) ) {
#   #   mat_var <- t( sapply( l_results[[paste0( 'lMatSigma', est_name, 'Pt' )]][-1], diag ) )
#   #   matResults[,4] <- sqrt( apply( mat_var, 2, mean ) )
#   #   matResults[,5] <- apply( t( mat_coef ) > trueValues - qnorm(0.975) * t( sqrt( mat_var ) ) &
#   #                              t( mat_coef ) < trueValues + qnorm(0.975) * t( sqrt( mat_var ) ), 1, mean )
#   # }
#   # colnames( matResults ) <- c("Mean", "sqrtMSE", "SDemp", "SDmean", "CovP")
#   colnames( matResults ) <- c("Mean", "sqrtMSE", "SDemp")
#   return( matResults )
# }
#
# runall <- function(b, # regression parameters.
#                    a, sd_eG2, dist_eG2,  # for specifying dist'n of mediation variables
#                    l_cov_eG2=NULL,
#                    dist_eT, sd_eT, mean_eT, # dist'n of error term in logT
#                    n1=150, n2=150,
#                    percCens=0.1,
#                    perturb=T,
#                    rep=500, B=500,
#                    quiet=F,
#                    track=T) {
#
#   n <- c(n1, n2)
#   b <- matrix(b, ncol=1)
#   if( is.vector( a ) ) a <- matrix( a, ncol=1 )
#   pall <- nrow( b ) # Number of dimensions.
#   pG <- ncol( a )
#   pZ <- nrow( a )
#   if( pG + pZ != pall ) stop( 'a should be pZ x pG, and b should be pZ + pG!' )
#
#   wt0 <- rep( 1, sum(n) )
#
#   lDfData <- list( ) # List of simulated datasets
#   lMatWt <- list( ) # List of perturbation weight matrices.
#
#   matBeta1 <- matrix( NA, rep, pall ) # Naive estimator
#   matBeta2 <- matrix( NA, rep, pall ) # Combined estimator
#   lMatSigmaPt <- list( ) # List of covariance matrices estimated from perturbation.
#   lMatSigma1Pt <- list( )
#   lMatSigma2Pt <- list( )
#
#   matBetaStar <- matrix( NA, rep, pall ) # Weighted estiator
#   lMatSigmaStarPt <- list( ) # List of covariance matrices of betaStar estimated from perturbation.
#
#   matBetaNan <- matrix( NA, rep, pall ) # Nan estimator, estimated IPW
#   lMatSigmaNanPt <- list( ) # List of covariance matrices of Nan estimated from perturbation.
#   # matBetaNan2 <- matrix( NA, rep, pall ) # Nan estimator, true IPW
#
#   matBetaMI <- matrix( NA, rep, pall ) # Multiple Imputation estimator
#   matBetaMI_tmp <- matrix( NA, 5, pall ) # for multiple imputation results
#
#   matBetaFib1 <- matrix( NA, rep, 1 ) # Fib estimator
#   lMatSigmaFib1Pt <- list( ) # List of covariance matrices of Fib estimated from perturbation.
#
#   matBetaFib2 <- matrix( NA, rep, pZ ) # Fib estimator
#   lMatSigmaFib2Pt <- list( ) # List of covariance matrices of Fib estimated from perturbation.
#
#   matBetaFib3 <- matrix( NA, rep, pall ) # Fib estimator
#   lMatSigmaFib3Pt <- list( ) # List of covariance matrices of Fib estimated from perturbation.
#
#   for( i in 1:rep ) {
#
#     if( !quiet & (i %% 20 == 0) ) cat( round( i / rep * 100, 0 ), '% ', sep='' )
#
#     dfData <- truedata(b=b, # regression parameters.
#                        a=a, sd_eG2=sd_eG2, dist_eG2=dist_eG2,  # for specifying dist'n of mediation variables
#                        l_cov_eG2=l_cov_eG2,
#                        dist_eT=dist_eT, sd_eT=sd_eT, mean_eT=mean_eT, # dist'n of error term in logT
#                        n1=n1, n2=n2,
#                        percCens=percCens)
#     # Sanity check to see if the data simulating function is correct.
#     # all( apply( dfData, 1, function(x){1*(x['logT'] == x['logX']) == x['delta']}) )
#     # all( apply( dfData, 1, function(x){ 1-is.na( x['G'] ) == x['R']}) )
#
#     matBeta1[i, ] <- gehan.fit1( y=dfData$logX, delta=dfData$delta, matX=as.matrix( dfData[, 1:pall] ),
#                                  wt=wt0, n1=n1, beta.ini=b[, 1] )
#     matBeta2[i, ] <- gehan.fit2( y=dfData$logX, delta=dfData$delta,
#                                  matZ=as.matrix( dfData[, 1:pZ] ), matG=as.matrix( dfData[, (pZ+1):(pZ+pG)] ),
#                                  wt=wt0, n1=n1, n2=n2, beta.ini=b[, 1] )
#
#     wtIp <- 1 / glm( dfData$R1 ~ as.matrix( cbind( dfData[, 1:pZ], dfData$logX, dfData$delta ) ), family=binomial() )$fitted.values
#     matBetaNan[i, ] <- gehan.fit1( y=dfData$logX, delta=dfData$delta,
#                                    matX=as.matrix( dfData[, 1:pall] ), wt=wtIp,
#                                    n1=n1, beta.ini=b[, 1] )
#     # matBetaNan2[i, ] <- gehan.fit1( y=dfData$logX, delta=dfData$delta, matX=as.matrix( dfData[, 1:pall] ), wt=1 / dfData$probMiss, n1=n1, beta.ini=b )
#
#
#
#
#     # non-weighted gammma estimates for calculating Fib estimator
#     gammaFib1 <- gehan.fit1( y=dfData$logX, delta=dfData$delta,
#                              matX=as.matrix( dfData[, 1:pZ] ),
#                              n1=n1, wt=wt0,
#                              beta.ini=(b[1:pZ, ] + a%*%b[(pZ+1):pall, ])[, 1] )
#     dfData2 <- dfData[(n1 + 1):sum(n), ]
#     gammaFib2 <- gehan.fit1( y=dfData2$logX, delta=dfData2$delta,
#                              matX=as.matrix( dfData2[, 1:pZ] ),
#                              n1=n2, wt=wt0,
#                              beta.ini=(b[1:pZ, ] + a%*%b[(pZ+1):pall, ])[, 1] )
#
#     if ( perturb ) {
#       # beta star, as well as its variance
#       if ( i >= 2 ) {
#         # Calculate the estimated variance matrix of betaStar from perturbation.
#         # The variance comes from current replicate and the weight comes from
#         # the previous run (matW hasn't been updated yet).
#         lMatSigmaStarPt[[i]] <- matW %*% matSigmaA %*% t( matW ) + matW %*% matSigmaB + t( matW %*% matSigmaB ) + matSigmaPt[(pall+1):(2*pall), (pall+1):(2*pall)]
#       } else lMatSigmaStarPt[[i]] <- NULL
#       lMatSigmaPt[[i]] <- matSigmaPt
#       lMatSigma1Pt[[i]] <- matSigmaPt[1:pall, 1:pall]
#       lMatSigma2Pt[[i]] <- matSigmaPt[(pall+1):(2*pall), (pall+1):(2*pall)]
#
#       # covariance for Nan estimator
#       matBetaNanPt <- matrix( NA, B, pall )
#       for( iPt in 1:B ) matBetaNanPt[iPt, ] <- gehan.fit1( y=dfData$logX, delta=dfData$delta,
#                                                            matX=as.matrix( dfData[, 1:pall] ),
#                                                            wt=matWt[, iPt] * 1 / glm( dfData$R1 ~ as.matrix( cbind( dfData[, 1:pZ], dfData$logX, dfData$delta ) ), family=binomial(), weights = matWt[, iPt] )$fitted.values,
#                                                            n1=n1, beta.ini=b[, 1] )
#       lMatSigmaNanPt[[i]] <- cov( matBetaNanPt )
#
#
#     }
#     if ( track ) {
#       lDfData[[i]] <- dfData
#       if ( perturb ) lMatWt[[i]] <- matWt
#     }
#   }
#
#   to.return <- list( matBeta1=matBeta1, matBeta2=matBeta2, matBetaStar=matBetaStar,
#                      matBetaNan=matBetaNan, matBetaMI=matBetaMI,
#                      matBetaFib=matBetaFib,
#                      lMatSigma1Pt=lMatSigma1Pt, lMatSigma2Pt=lMatSigma2Pt,
#                      lMatSigmaStarPt=lMatSigmaStarPt, lMatSigmaNanPt=lMatSigmaNanPt,
#                      lMatSigmaFibPt=lMatSigmaFibPt,
#                      setup=list( beta=b, # regression parameters.
#                                  alpha=a, sd_eG2, dist_eG2,  # for specifying dist'n of mediation variables
#                                  l_cov_eG2=NULL,
#                                  dist_eT, sd_eT, mean_eT, # dist'n of error term in logT
#                                  n1=150, n2=150,
#                                  percCens=0.1,
#                                  perturb=T,
#                                  rep=500, B=500,
#                                  quiet=F,
#                                  track=T ) )
#   if( track ) to.return <- c( to.return, list( lDfData=lDfData, lMatWt=lMatWt ) )
#   return( to.return )
# }
#
# runall_fib <- function( b, # regression parameters.
#                         a, sd_eG2, dist_eG2,  # for specifying dist'n of mediation variables
#                         l_cov_eG2=NULL,
#                         dist_eT, sd_eT, mean_eT, # dist'n of error term in logT
#                         n1=150, n2=150,
#                         percCens=0.1,
#                         perturb=T,
#                         rep=500, B=500,
#                         quiet=F,
#                         track=T ) {
#   n <- c(n1, n2)
#   b <- matrix( b, ncol=1 )
#   if( is.vector( a ) ) a <- matrix( a, ncol=1 )
#   pall <- nrow( b ) # Number of dimensions.
#   pG <- ncol( a )
#   pZ <- nrow( a )
#   if( pG + pZ != pall ) stop( 'a should be pZ x pG, and b should be pZ + pG!' )
#
#   wt0 <- rep( 1, sum(n) )
#
#   lDfData <- list( ) # List of simulated datasets
#   lMatWt <- list( ) # List of perturbation weight matrices.
#
#   matBeta1 <- matrix( NA, rep, pall ) # Naive estimator
#   matBeta2 <- matrix( NA, rep, pall ) # Combined estimator
#   # lMatSigmaPt <- list( ) # List of covariance matrices estimated from perturbation.
#   # lMatSigma1Pt <- list( )
#   # lMatSigma2Pt <- list( )
#
#   matBetaStar <- matrix( NA, rep, pall ) # Weighted estiator
#   # lMatSigmaStarPt <- list( ) # List of covariance matrices of betaStar estimated from perturbation.
#
#   # matBetaNan <- matrix( NA, rep, pall ) # Nan estimator, estimated IPW
#   # lMatSigmaNanPt <- list( ) # List of covariance matrices of Nan estimated from perturbation.
#   # # matBetaNan2 <- matrix( NA, rep, pall ) # Nan estimator, true IPW
#   #
#   # matBetaMI <- matrix( NA, rep, pall ) # Multiple Imputation estimator
#   # matBetaMI_tmp <- matrix( NA, 5, pall ) # for multiple imputation results
#
#   matBetaFib1 <- matrix( NA, rep, 1 ) # Fib estimator
#   # lMatSigmaFib1Pt <- list( ) # List of covariance matrices of Fib estimated from perturbation.
#
#   matBetaFib2 <- matrix( NA, rep, pZ ) # Fib estimator
#   # lMatSigmaFib2Pt <- list( ) # List of covariance matrices of Fib estimated from perturbation.
#
#   matBetaFib3 <- matrix( NA, rep, pZ ) # Fib estimator
#   # lMatSigmaFib3Pt <- list( ) # List of covariance matrices of Fib estimated from perturbation.
#
#   matBetaFib4 <- matrix( NA, rep, pall ) # Fib estimator
#   # lMatSigmaFib4Pt <- list( ) # List of covariance matrices of Fib estimated from perturbation.
#
#   matGamma1 <- matrix( NA, rep, pZ ) # Naive estimator
#   matGamma2 <- matrix( NA, rep, pZ ) # Combined estimator
#
#   for( i in 1:rep ) {
#
#     if( !quiet & (i %% 20 == 0) ) cat( round( i / rep * 100, 0 ), '% ', sep='' )
#
#     dfData <- truedata( b=b, # regression parameters.
#                         a=a, sd_eG2=sd_eG2, dist_eG2=dist_eG2,  # for specifying dist'n of mediation variables
#                         l_cov_eG2=l_cov_eG2,
#                         dist_eT=dist_eT, sd_eT=sd_eT, mean_eT=mean_eT, # dist'n of error term in logT
#                         n1=n1, n2=n2,
#                         percCens=percCens )
#     # Sanity check to see if the data simulating function is correct.
#     # all( apply( dfData, 1, function(x){1*(x['logT'] == x['logX']) == x['delta']}) )
#     # all( apply( dfData, 1, function(x){ 1-is.na( x['G'] ) == x['R']}) )
#
#     matBeta1[i, ] <- gehan.fit1( y=dfData$logX, delta=dfData$delta, matX=as.matrix( dfData[, 1:pall] ),
#                                  wt=wt0, n1=n1, beta.ini=b[, 1] )
#     matBeta2[i, ] <- gehan.fit2( y=dfData$logX, delta=dfData$delta,
#                                  matZ=as.matrix( dfData[, 1:pZ] ), matG=as.matrix( dfData[, (pZ+1):(pZ+pG)] ),
#                                  wt=wt0, n1=n1, n2=n2, beta.ini=b[, 1] )
#
#     # wtIp <- 1 / glm( dfData$R1 ~ as.matrix( cbind( dfData[, 1:pZ], dfData$logX, dfData$delta ) ), family=binomial() )$fitted.values
#     # matBetaNan[i, ] <- gehan.fit1( y=dfData$logX, delta=dfData$delta,
#     #                                matX=as.matrix( dfData[, 1:pall] ), wt=wtIp,
#     #                                n1=n1, beta.ini=b[, 1] )
#     # matBetaNan2[i, ] <- gehan.fit1( y=dfData$logX, delta=dfData$delta, matX=as.matrix( dfData[, 1:pall] ), wt=1 / dfData$probMiss, n1=n1, beta.ini=b )
#
#     # md_tmp <- mice( dfData, meth='norm',
#     #                 predictorMatrix=rbind( matrix( 0, nrow=pZ, ncol=ncol( dfData ) ),
#     #                                        cbind( matrix( 1, nrow=pG, ncol=pZ ),
#     #                                               1 - diag( 1, nrow=pG, ncol=pG ),
#     #                                               matrix( 0, nrow=pG, ncol=2 ),
#     #                                               matrix( 1, nrow=pG, ncol=1 ),
#     #                                               matrix( 0, nrow=pG, ncol=ncol( dfData ) - pZ - pG - 3 )
#     #                                        ),
#     #                                        matrix( 0, nrow=ncol(dfData) - pZ - pG, ncol=ncol( dfData ) ) ),
#     #                 printFlag=F )
#     # for ( j in 1:5 ) {
#     #   df_tmp <- dfData
#     #   df_tmp[(n1+1):(n1+n2), (pZ + 1):(pZ + pG)] <- Reduce( 'cbind', lapply( md_tmp$imp[paste0( 'G', 1:pG )], function(x)x[, j] ) )
#     #   matBetaMI_tmp[j, ] <- gehan.fit1( y=dfData$logX, delta=dfData$delta, matX=as.matrix( df_tmp[, 1:pall] ),
#     #                                     wt=wt0, n1=sum(n), beta.ini=b[, 1] )
#     # }
#     # matBetaMI[i, ] <- apply( matBetaMI_tmp, 2, mean )
#
#     # non-weighted gammma estimates for calculating Fib estimator
#     gammaFib1 <- gehan.fit1( y=dfData$logX, delta=dfData$delta,
#                              matX=as.matrix( dfData[, 1:pZ] ),
#                              n1=n1, wt=wt0,
#                              beta.ini=(b[1:pZ, ] + a%*%b[(pZ+1):pall, ])[, 1] )
#     dfData2 <- dfData[(n1 + 1):sum(n), ]
#     gammaFib2 <- gehan.fit1( y=dfData2$logX, delta=dfData2$delta,
#                              matX=as.matrix( dfData2[, 1:pZ] ),
#                              n1=n2, wt=wt0,
#                              beta.ini=(b[1:pZ, ] + a%*%b[(pZ+1):pall, ])[, 1] )
#     matGamma1[i, ] <- gammaFib1
#     matGamma2[i, ] <- gammaFib2
#
#     if ( perturb ) {
#       matWt <- matrix( rexp( sum(n)*B ), nrow=sum(n), ncol=B )
#       matBeta1Pt <- perturbfn1( y=dfData$logX, delta=dfData$delta, matX=as.matrix( dfData[, 1:pall] ),
#                                 n1=n1, matWt=matWt, B=B, beta.ini=b[, 1] )
#       matBeta2Pt <- perturbfn2( y=dfData$logX, delta=dfData$delta,
#                                 matZ=as.matrix( dfData[, 1:pZ] ), matG=as.matrix( dfData[, (pZ+1):(pZ+pG)] ),
#                                 n1=n1, n2=n2, matWt=matWt, B=B, beta.ini=b[, 1] )
#       matBetaPt <- cbind( matBeta1Pt, matBeta2Pt )
#       matSigmaPt <- cov( matBetaPt )
#       matSigmaA <- matSigmaPt[1:pall, 1:pall] + matSigmaPt[(pall+1):(2*pall), (pall+1):(2*pall)] - matSigmaPt[1:pall, (pall+1):(2*pall)] - matSigmaPt[(pall+1):(2*pall), 1:pall]
#       matSigmaB <- matSigmaPt[1:pall, (pall+1):(2*pall)] - matSigmaPt[(pall+1):(2*pall), (pall+1):(2*pall)]
#
#       # beta star, as well as its variance
#       # if ( i >= 2 ) {
#       #   # Calculate the estimated variance matrix of betaStar from perturbation.
#       #   # The variance comes from current replicate and the weight comes from
#       #   # the previous run (matW hasn't been updated yet).
#       #   lMatSigmaStarPt[[i]] <- matW %*% matSigmaA %*% t( matW ) + matW %*% matSigmaB + t( matW %*% matSigmaB ) + matSigmaPt[(pall+1):(2*pall), (pall+1):(2*pall)]
#       # } else lMatSigmaStarPt[[i]] <- NULL
#       matW <- - t(matSigmaB) %*% solve( matSigmaA )
#       matBetaStar[i, ] <- matW %*% matBeta1[i, ] + ( diag( 1, pall ) - matW ) %*% matBeta2[i, ]
#       # lMatSigmaPt[[i]] <- matSigmaPt
#       # lMatSigma1Pt[[i]] <- matSigmaPt[1:pall, 1:pall]
#       # lMatSigma2Pt[[i]] <- matSigmaPt[(pall+1):(2*pall), (pall+1):(2*pall)]
#       #
#       # covariance for Nan estimator
#       # matBetaNanPt <- matrix( NA, B, pall )
#       # for( iPt in 1:B ) matBetaNanPt[iPt, ] <- gehan.fit1( y=dfData$logX, delta=dfData$delta,
#       #                                                      matX=as.matrix( dfData[, 1:pall] ),
#       #                                                      wt=matWt[, iPt] * 1 / glm( dfData$R1 ~ as.matrix( cbind( dfData[, 1:pZ], dfData$logX, dfData$delta ) ), family=binomial(), weights = matWt[, iPt] )$fitted.values,
#       #                                                      n1=n1, beta.ini=b[, 1] )
#       # lMatSigmaNanPt[[i]] <- cov( matBetaNanPt )
#
#       # betaFib, as well as its variance
#       matGammaPt <- perturbfn1( y=dfData$logX, delta=dfData$delta, matX=as.matrix( dfData[, 1:pZ] ),
#                                 n1=n1, matWt=matWt, B=B, beta.ini=(b[1:pZ, ] + a%*%b[(pZ+1):pall, ])[, 1] )
#       matGamma2Pt <- perturbfn1( y=dfData2$logX, delta=dfData2$delta, matX=as.matrix( dfData2[, 1:pZ] ),
#                                  n1=n2, matWt=matWt[(n1 + 1):sum(n), ], B=B, beta.ini=(b[1:pZ, ] + a%*%b[(pZ+1):pall, ])[, 1] )
#       matSigmaBetaGammaPt <- cov( cbind( matBeta1Pt, matGammaPt ) )
#       matSigmaGamma2Pt <- cov( matGamma2Pt )
#
#       # Fib1, the original estimator from Fibrogen group
#       p_beta_tmp <- 1
#       p_gamma_tmp <- 1
#       matSigmaBetaGammaPt_tmp <- matSigmaBetaGammaPt[c( 1:p_beta_tmp, (pall + 1):(pall + p_gamma_tmp) ),
#                                                      c( 1:p_beta_tmp, (pall + 1):(pall + p_gamma_tmp) )]
#       matSigmaBetaGammaPtInv_tmp <- solve( matSigmaBetaGammaPt_tmp )
#       matSigmaGamma2Pt_tmp <- matrix( matSigmaGamma2Pt[1:p_gamma_tmp, 1:p_gamma_tmp], p_gamma_tmp, p_gamma_tmp )
#       matB_tmp <- solve( matSigmaBetaGammaPtInv_tmp[(p_beta_tmp + 1):(p_beta_tmp + p_gamma_tmp),
#                                                     (p_beta_tmp + 1):(p_beta_tmp + p_gamma_tmp)] +
#                            solve( matSigmaGamma2Pt_tmp ) )
#       matA_tmp <- solve( matSigmaBetaGammaPtInv_tmp[1:p_beta_tmp, 1:p_beta_tmp] -
#                            matSigmaBetaGammaPtInv_tmp[1:p_beta_tmp, (p_beta_tmp+1):(p_beta_tmp + p_gamma_tmp)] %*%
#                            matB_tmp %*% matSigmaBetaGammaPtInv_tmp[(p_beta_tmp+1):(p_beta_tmp + p_gamma_tmp), 1:p_beta_tmp],
#                          matSigmaBetaGammaPtInv_tmp[1:p_beta_tmp, (p_beta_tmp+1):(p_beta_tmp + p_gamma_tmp)] %*%
#                            matB_tmp %*% solve( matSigmaGamma2Pt_tmp ) )
#       matBetaFib1[i, ] <- matA_tmp%*%( gammaFib1[1:p_gamma_tmp] - gammaFib2[1:p_gamma_tmp] ) +
#         matBeta1[i, 1:p_beta_tmp]
#
#       # Fib2, estimator using all observed covariates and only the first covariate of gamma
#       p_beta_tmp <- pZ
#       p_gamma_tmp <- 1
#       matSigmaBetaGammaPt_tmp <- matSigmaBetaGammaPt[c( 1:p_beta_tmp, (pall + 1):(pall + p_gamma_tmp) ),
#                                                      c( 1:p_beta_tmp, (pall + 1):(pall + p_gamma_tmp) )]
#       matSigmaBetaGammaPtInv_tmp <- solve( matSigmaBetaGammaPt_tmp )
#       matSigmaGamma2Pt_tmp <- matrix( matSigmaGamma2Pt[1:p_gamma_tmp, 1:p_gamma_tmp], p_gamma_tmp, p_gamma_tmp )
#       matB_tmp <- solve( matSigmaBetaGammaPtInv_tmp[(p_beta_tmp + 1):(p_beta_tmp + p_gamma_tmp),
#                                                     (p_beta_tmp + 1):(p_beta_tmp + p_gamma_tmp)] +
#                            solve( matSigmaGamma2Pt_tmp ) )
#       matA_tmp <- solve( matSigmaBetaGammaPtInv_tmp[1:p_beta_tmp, 1:p_beta_tmp] -
#                            matSigmaBetaGammaPtInv_tmp[1:p_beta_tmp, (p_beta_tmp+1):(p_beta_tmp + p_gamma_tmp)] %*%
#                            matB_tmp %*% matSigmaBetaGammaPtInv_tmp[(p_beta_tmp+1):(p_beta_tmp + p_gamma_tmp), 1:p_beta_tmp],
#                          matSigmaBetaGammaPtInv_tmp[1:p_beta_tmp, (p_beta_tmp+1):(p_beta_tmp + p_gamma_tmp)] %*%
#                            matB_tmp %*% solve( matSigmaGamma2Pt_tmp ) )
#       matBetaFib2[i, ] <- matA_tmp%*%( gammaFib1[1:p_gamma_tmp] - gammaFib2[1:p_gamma_tmp] ) +
#         matBeta1[i, 1:p_beta_tmp]
#
#       # Fib3, estimator using all observed covariates and all of gamma
#       p_beta_tmp <- pZ
#       p_gamma_tmp <- pZ
#       matSigmaBetaGammaPt_tmp <- matSigmaBetaGammaPt[c( 1:p_beta_tmp, (pall + 1):(pall + p_gamma_tmp) ),
#                                                      c( 1:p_beta_tmp, (pall + 1):(pall + p_gamma_tmp) )]
#       matSigmaBetaGammaPtInv_tmp <- solve( matSigmaBetaGammaPt_tmp )
#       matSigmaGamma2Pt_tmp <- matrix( matSigmaGamma2Pt[1:p_gamma_tmp, 1:p_gamma_tmp], p_gamma_tmp, p_gamma_tmp )
#       matB_tmp <- solve( matSigmaBetaGammaPtInv_tmp[(p_beta_tmp + 1):(p_beta_tmp + p_gamma_tmp),
#                                                     (p_beta_tmp + 1):(p_beta_tmp + p_gamma_tmp)] +
#                            solve( matSigmaGamma2Pt_tmp ) )
#       matA_tmp <- solve( matSigmaBetaGammaPtInv_tmp[1:p_beta_tmp, 1:p_beta_tmp] -
#                            matSigmaBetaGammaPtInv_tmp[1:p_beta_tmp, (p_beta_tmp+1):(p_beta_tmp + p_gamma_tmp)] %*%
#                            matB_tmp %*% matSigmaBetaGammaPtInv_tmp[(p_beta_tmp+1):(p_beta_tmp + p_gamma_tmp), 1:p_beta_tmp],
#                          matSigmaBetaGammaPtInv_tmp[1:p_beta_tmp, (p_beta_tmp+1):(p_beta_tmp + p_gamma_tmp)] %*%
#                            matB_tmp %*% solve( matSigmaGamma2Pt_tmp ) )
#       matBetaFib3[i, ] <- matA_tmp%*%( gammaFib1[1:p_gamma_tmp] - gammaFib2[1:p_gamma_tmp] ) +
#         matBeta1[i, 1:p_beta_tmp]
#
#       # Fib4, full estimator
#       p_beta_tmp <- pZ + pG
#       p_gamma_tmp <- pZ
#       matSigmaBetaGammaPt_tmp <- matSigmaBetaGammaPt[c( 1:p_beta_tmp, (pall + 1):(pall + p_gamma_tmp) ),
#                                                      c( 1:p_beta_tmp, (pall + 1):(pall + p_gamma_tmp) )]
#       matSigmaBetaGammaPtInv_tmp <- solve( matSigmaBetaGammaPt_tmp )
#       matSigmaGamma2Pt_tmp <- matrix( matSigmaGamma2Pt[1:p_gamma_tmp, 1:p_gamma_tmp], p_gamma_tmp, p_gamma_tmp )
#       matB_tmp <- solve( matSigmaBetaGammaPtInv_tmp[(p_beta_tmp + 1):(p_beta_tmp + p_gamma_tmp),
#                                                     (p_beta_tmp + 1):(p_beta_tmp + p_gamma_tmp)] +
#                            solve( matSigmaGamma2Pt_tmp ) )
#       matA_tmp <- solve( matSigmaBetaGammaPtInv_tmp[1:p_beta_tmp, 1:p_beta_tmp] -
#                            matSigmaBetaGammaPtInv_tmp[1:p_beta_tmp, (p_beta_tmp+1):(p_beta_tmp + p_gamma_tmp)] %*%
#                            matB_tmp %*% matSigmaBetaGammaPtInv_tmp[(p_beta_tmp+1):(p_beta_tmp + p_gamma_tmp), 1:p_beta_tmp],
#                          matSigmaBetaGammaPtInv_tmp[1:p_beta_tmp, (p_beta_tmp+1):(p_beta_tmp + p_gamma_tmp)] %*%
#                            matB_tmp %*% solve( matSigmaGamma2Pt_tmp ) )
#       matBetaFib4[i, ] <- matA_tmp%*%( gammaFib1[1:p_gamma_tmp] - gammaFib2[1:p_gamma_tmp] ) +
#         matBeta1[i, 1:p_beta_tmp]
#
#       # lMatSigmaFibPt[[i]] <- matSigmaBetaGammaPt[1:pall, 1:pall] +
#       #   matA %*% matSigmaBetaGammaPt[(pall+1):(pall + pZ), (pall+1):(pall + pZ)] %*% t( matA ) +
#       #   matA %*% matSigmaBetaGammaPt[(pall+1):(pall + pZ), 1:pall] +
#       #   matSigmaBetaGammaPt[1:pall, (pall+1):(pall + pZ)] %*% t( matA ) +
#       #   matA %*% matSigmaGamma2Pt %*% t( matA )
#     }
#     if ( track ) {
#       lDfData[[i]] <- dfData
#       if ( perturb ) lMatWt[[i]] <- matWt
#     }
#   }
#
#   to.return <- list( matBeta1=matBeta1, matBeta2=matBeta2, matBetaStar=matBetaStar,
#                      # matBetaNan=matBetaNan, matBetaMI=matBetaMI,
#                      matBetaFib1=matBetaFib1, matBetaFib2=matBetaFib2,
#                      matBetaFib3=matBetaFib3, matBetaFib4=matBetaFib4,
#                      matGamma1=matGamma1, matGamma2=matGamma2,
#                      # lMatSigma1Pt=lMatSigma1Pt, lMatSigma2Pt=lMatSigma2Pt,
#                      # lMatSigmaStarPt=lMatSigmaStarPt,
#                      # # lMatSigmaNanPt=lMatSigmaNanPt,
#                      # lMatSigmaFibPt=lMatSigmaFibPt,
#                      setup=list( beta=b, # regression parameters.
#                                  alpha=a, sd_eG2, dist_eG2,  # for specifying dist'n of mediation variables
#                                  l_cov_eG2=NULL,
#                                  dist_eT, sd_eT, mean_eT, # dist'n of error term in logT
#                                  n1=150, n2=150,
#                                  percCens=0.1,
#                                  perturb=T,
#                                  rep=500, B=500,
#                                  quiet=F,
#                                  track=T ) )
#   if( track ) to.return <- c( to.return, list( lDfData=lDfData, lMatWt=lMatWt ) )
#   return( to.return )
# }
#

# runall.v2 <- function( n=c( 150, 150 ), a=c( 1, 1 ), b=c( 1, 1, 1 ),
#                     percCens=0.1, sdG=1, rep=500, B=500,
#                     quiet=F, perturb=T, sameErrDist=T, disteG='normal',
#                     track=T ) {
#   # Return: list of betas in all replicates?
#
#   if ( length( n ) > 1 ) {
#     n1 <- n[1]
#     n2 <- n[2]
#     n <- sum( n )
#   } else stop( 'n should be vector of length 2!' )
#
#   p <- length( b ) # Number of dimensions.
#   wt0 <- rep( 1, n )
#
#   lDfData <- list( ) # List of simulated datasets
#   lMatWt <- list( ) # List of perturbation weight matrices.
#
#   matBeta1 <- matrix( NA, rep, p ) # Naive estimator
#   matBeta2 <- matrix( NA, rep, p ) # Combined estimator
#   lMatSigmaPt <- list( ) # List of covariance matrices estimated from perturbation.
#   lMatSigma1Pt <- list( )
#   lMatSigma2Pt <- list( )
#
#   matBetaStar <- matrix( NA, rep, p ) # Weighted estiator
#   lMatSigmaStarPt <- list( ) # List of covariance matrices of betaStar estimated from perturbation.
#
#   matBetaNan <- matrix( NA, rep, p ) # Nan estimator, estimated IPW
#   lMatSigmaNanPt <- list( ) # List of covariance matrices of Nan estimated from perturbation.
#   # matBetaNan2 <- matrix( NA, rep, p ) # Nan estimator, true IPW
#
#
#   matBetaFib <- matrix( NA, rep, p ) # Fib estimator
#   lMatSigmaFibPt <- list( ) # List of covariance matrices of Fib estimated from perturbation.
#
#   for( i in 1:rep ) {
#
#     if( !quiet & (i %% 20 == 0) ) cat( round( i / rep * 100, 0 ), '% ', sep='' )
#
#     dfData <- truedata( a=a, b=b, n1=n1, n2=n2, percCens=percCens,
#                         sdG=sdG, sameErrDist=sameErrDist, disteG=disteG )
#     # Sanity check to see if the data simulating function is correct.
#     # all( apply( dfData, 1, function(x){1*(x['logT'] == x['logX']) == x['delta']}) )
#     # all( apply( dfData, 1, function(x){ 1-is.na( x['G'] ) == x['R']}) )
#
#     md_tmp <- mice( dfData, meth='norm', printFlag=F )
#     for ( j in 1:5 ) {
#       df_tmp <- dfData
#       df_tmp[(n1+1):(n1+n2), ]$G <- md_tmp$imp$G[, j]
#       matBeta_tmp[j, ] <- gehan.fit1( y=dfData$logX, delta=dfData$delta, matX=as.matrix( df_tmp[, 1:p] ),
#                                       wt=wt0, n1=n, beta.ini=b )
#     }
#     matBeta[i, ] <- apply( matBeta_tmp, 2, mean )
#   }
#
#
#   for( i in 1:rep ) {
#
#     if( !quiet & (i %% 20 == 0) ) cat( round( i / rep * 100, 0 ), '% ', sep='' )
#
#     dfData <- truedata( a=a, b=b, n1=n1, n2=n2, percCens=percCens,
#                         sdG=sdG, sameErrDist=sameErrDist, disteG=disteG )
#     # Sanity check to see if the data simulating function is correct.
#     # all( apply( dfData, 1, function(x){1*(x['logT'] == x['logX']) == x['delta']}) )
#     # all( apply( dfData, 1, function(x){ 1-is.na( x['G'] ) == x['R']}) )
#
#     matBeta1[i, ] <- gehan.fit1( y=dfData$logX, delta=dfData$delta, matX=as.matrix( dfData[, 1:p] ),
#                                  wt=wt0, n1=n1, beta.ini=b )
#     matBeta2[i, ] <- gehan.fit2( y=dfData$logX, delta=dfData$delta, matX=as.matrix( dfData[, 1:p] ),
#                                  wt=wt0, n1=n1, n2=n2, beta.ini=b )
#
#     wtIp <- 1 / glm( dfData$R ~ as.matrix( cbind( dfData[, 1:(p-1)], dfData$logX, dfData$delta ) ), family=binomial() )$fitted.values
#     matBetaNan[i, ] <- gehan.fit1( y=dfData$logX, delta=dfData$delta,
#                                    matX=as.matrix( dfData[, 1:p] ), wt=wtIp,
#                                    n1=n1, beta.ini=b )
#     # matBetaNan2[i, ] <- gehan.fit1( y=dfData$logX, delta=dfData$delta, matX=as.matrix( dfData[, 1:p] ), wt=1 / dfData$probMiss, n1=n1, beta.ini=b )
#
#     if ( perturb ) {
#       matWt <- matrix( rexp( n*B ), nrow=n, ncol=B )
#       matBeta1Pt <- perturbfn1( y=dfData$logX, delta=dfData$delta, matX=as.matrix( dfData[, 1:p] ),
#                                 n1=n1, matWt=matWt, B=B, beta.ini=b )
#       matBeta2Pt <- perturbfn2( y=dfData$logX, delta=dfData$delta, matX=as.matrix( dfData[, 1:p] ), n1=n1, n2=n2, matWt=matWt, B=B, beta.ini=b )
#       matBetaPt <- cbind( matBeta1Pt, matBeta2Pt )
#       matSigmaPt <- cov( matBetaPt )
#       matSigmaA <- matSigmaPt[1:p, 1:p] + matSigmaPt[(p+1):(2*p), (p+1):(2*p)] - matSigmaPt[1:p, (p+1):(2*p)] - matSigmaPt[(p+1):(2*p), 1:p]
#       matSigmaB <- matSigmaPt[1:p, (p+1):(2*p)] - matSigmaPt[(p+1):(2*p), (p+1):(2*p)]
#
#       # beta star, as well as its variance
#       if ( i >= 2 ) {
#         # Calculate the estimated variance matrix of betaStar from perturbation.
#         # The variance comes from current replicate and the weight comes from
#         # the previous run (matW hasn't been updated yet).
#         lMatSigmaStarPt[[i]] <- matW %*% matSigmaA %*% t( matW ) + matW %*% matSigmaB + t( matW %*% matSigmaB ) + matSigmaPt[(p+1):(2*p), (p+1):(2*p)]
#       } else lMatSigmaStarPt[[i]] <- NULL
#       matW <- - matSigmaB %*% solve( matSigmaA )
#       matBetaStar[i, ] <- matW %*% matBeta1[i, ] + ( diag( 1, p ) - matW ) %*% matBeta2[i, ]
#       lMatSigmaPt[[i]] <- matSigmaPt
#       lMatSigma1Pt[[i]] <- matSigmaPt[1:p, 1:p]
#       lMatSigma2Pt[[i]] <- matSigmaPt[(p+1):(2*p), (p+1):(2*p)]
#
#       # covariance for Nan estimator
#       matBetaNanPt <- matrix( NA, B, p )
#       for( iPt in 1:B ) matBetaNanPt[iPt, ] <- gehan.fit1( y=dfData$logX, delta=dfData$delta,
#                                                        matX=as.matrix( dfData[, 1:p] ),
#                                                     wt=matWt[, iPt] * 1 / glm( dfData$R ~ as.matrix( cbind( dfData[, 1:(p-1)], dfData$logX, dfData$delta ) ), family=binomial(), weights = matWt[, iPt] )$fitted.values,
#                                                     n1=n1, beta.ini=b )
#       lMatSigmaNanPt[[i]] <- cov( matBetaNanPt )
#
#       # betaFib, as well as its variance
#       matGammaPt <- perturbfn1( y=dfData$logX, delta=dfData$delta, matX=as.matrix( dfData[, 1:(p - 1)] ),
#                                 n1=n1, matWt=matWt, B=B, beta.ini=b[1:(p-1)] + b[p]*a )
#       dfData2 <- dfData[(n1 + 1):n, ]
#       matGamma2Pt <- perturbfn1( y=dfData2$logX, delta=dfData2$delta, matX=as.matrix( dfData2[, 1:(p - 1)] ),
#                                 n1=n2, matWt=matWt[(n1 + 1):n, ], B=B, beta.ini=b[1:(p-1)] + b[p]*a )
#       matSigmaBetaGammaPt <- cov( cbind( matBetaPt, matGammaPt ) )
#       matSigmaGamma2Pt <- cov( matGamma2Pt )
#       matSigmaBetaGammaPtInv <- solve( matSigmaBetaGammaPt )
#       matB <- solve( matSigmaBetaGammaPtInv[(p+1):(2*p - 1), (p+1):(2*p - 1)] +
#                        solve( matSigmaGamma2Pt ) )
#       matA <- solve( matSigmaBetaGammaPtInv[1:p, 1:p] -
#                        matSigmaBetaGammaPtInv[1:p, (p+1):(2*p - 1)] %*%
#                        matB %*% matSigmaBetaGammaPtInv[(p+1):(2*p - 1), 1:p],
#                      matSigmaBetaGammaPtInv[1:p, (p+1):(2*p - 1)] %*%
#                        matB %*% solve( matSigmaGamma2Pt ) )
#       matBetaFib[i, ] <- matA%*%( gehan.fit1( y=dfData$logX, delta=dfData$delta,
#                                               matX=as.matrix( dfData[, 1:(p - 1)] ),
#                                               n1=n1, wt=wt0,
#                                               beta.ini=b[1:(p-1)] + b[p]*a ) -
#                                     gehan.fit1( y=dfData2$logX, delta=dfData2$delta,
#                                                 matX=as.matrix( dfData2[, 1:(p - 1)] ),
#                                                 n1=n2, wt=wt0,
#                                                 beta.ini=b[1:(p-1)] + b[p]*a )  ) +
#                          matBeta1[i, ]
#       lMatSigmaFibPt[[i]] <- matSigmaBetaGammaPt[1:p, 1:p] +
#                              matA %*% matSigmaBetaGammaPt[(p+1):(2*p - 1), (p+1):(2*p - 1)] %*% t( matA ) +
#                              matA %*% matSigmaBetaGammaPt[(p+1):(2*p - 1), 1:p] +
#                              matSigmaBetaGammaPt[1:p, (p+1):(2*p - 1)] %*% t( matA ) +
#                              matA %*% matSigmaGamma2Pt %*% t( matA )
#     }
#     if ( track ) {
#       lDfData[[i]] <- dfData
#       if ( perturb ) lMatWt[[i]] <- matWt
#     }
#   }
#
#   to.return <- list( matBeta1=matBeta1, matBeta2=matBeta2, matBetaStar=matBetaStar,
#                      matBetaNan=matBetaNan, matBetaFib=matBetaFib,
#                      lMatSigma1Pt=lMatSigma1Pt, lMatSigma2Pt=lMatSigma2Pt,
#                      lMatSigmaStarPt=lMatSigmaStarPt, lMatSigmaNanPt=lMatSigmaNanPt,
#                      lMatSigmaFibPt=lMatSigmaFibPt )
#   if( track ) to.return <- c( to.return, list( lDfData=lDfData, lMatWt=lMatWt ) )
#   return( to.return )
# }
#
# runall.mi <- function( n=c( 150, 150 ), a=c( 1, 1 ), b=c( 1, 1, 1 ),
#                        percCens=0.1, sdG=1, rep=500, B=500,
#                        quiet=F, perturb=T, sameErrDist=T, disteG='normal',
#                        track=T ) {
#   # Return: list of betas in all replicates?
#
#   if ( length( n ) > 1 ) {
#     n1 <- n[1]
#     n2 <- n[2]
#     n <- sum( n )
#   } else stop( 'n should be vector of length 2!' )
#
#   p <- length( b ) # Number of dimensions.
#   wt0 <- rep( 1, n )
#
#   lDfData <- list( ) # List of simulated datasets
#
#   matBeta <- matrix( NA, rep, p ) # Naive estimator
#   matBeta_tmp <- matrix( NA, 5, p ) # for multiple imputation results
#   lMatSigmaPt <- list( ) # List of covariance matrices estimated from imputation.
#
#   for( i in 1:rep ) {
#
#     if( !quiet & (i %% 20 == 0) ) cat( round( i / rep * 100, 0 ), '% ', sep='' )
#
#     dfData <- truedata( a=a, b=b, n1=n1, n2=n2, percCens=percCens,
#                         sdG=sdG, sameErrDist=sameErrDist, disteG=disteG )
#     # Sanity check to see if the data simulating function is correct.
#     # all( apply( dfData, 1, function(x){1*(x['logT'] == x['logX']) == x['delta']}) )
#     # all( apply( dfData, 1, function(x){ 1-is.na( x['G'] ) == x['R']}) )
#
#     md_tmp <- mice( dfData, meth='norm', printFlag=F )
#     for ( j in 1:5 ) {
#       df_tmp <- dfData
#       df_tmp[(n1+1):(n1+n2), ]$G <- md_tmp$imp$G[, j]
#       matBeta_tmp[j, ] <- gehan.fit1( y=dfData$logX, delta=dfData$delta, matX=as.matrix( df_tmp[, 1:p] ),
#                                    wt=wt0, n1=n, beta.ini=b )
#     }
#     matBeta[i, ] <- apply( matBeta_tmp, 2, mean )
#   }
#
#   return( list( matBetaMI=matBeta ) )
# }
# runall <- function( n=300, a=c( 1, 1 ), b=c( 1, 1, 1 ),
#                     eta=c(1, 1, 0, 0, 0), percCens=0.1, sdG=1, rep=500, B=500,
#                     quiet=F, perturb=F, missing=F,
#                     sameErrDist=!missing, betaNan=missing, track=F ) {
#   # Return: list of betas in all replicates?
#   # n should be 1 number when missing=T, and a vector of length2 when missing=F
#
#   if ( length( n ) > 1 ) {
#     n1 <- n[1]
#     n2 <- n[2]
#     n <- sum( n )
#   }
#
#   p <- length( b ) # Number of dimensions.
#   wt0 <- rep( 1, n )
#
#   matBeta1 <- matrix( NA, rep, p ) # Naive estimator
#   matBeta2 <- matrix( NA, rep, p ) # Combined estimator
#   lMatSigmaPt <- list( ) # List of covariance matrices estimated from perturbation.
#
#   matBetaStar <- matrix( NA, rep, p ) # Weighted estiator
#   lMatSigmaStarPt <- list( ) # List of covariance matrices of betaStar estimated from perturbation.
#
#   matBetaNan <- matrix( NA, rep, p ) # Nan estimator, estimated IPW
#   lMatSigmaStarPt <- list( ) # List of covariance matrices of betaStar estimated from perturbation.
#
#   # matBetaNan2 <- matrix( NA, rep, p ) # Nan estimator, true IPW
#
#   lDfData <- list( ) # List of simulated datasets
#   lMatWt <- list( ) # List of perturbation weight matrices.
#
#   for( i in 1:rep ) {
#
#     if( !quiet & (i %% 20 == 0) ) cat( round( i / rep * 100, 0 ), '% ', sep='' )
#
#     if ( !missing ) dfData <- truedata( a=a, b=b, n1=n1, n2=n2, percCens=percCens,
#                                         sdG=sdG, sameErrDist=sameErrDist )
#     else {
#       dfData <- truedataMiss( a=a, b=b, eta=eta, n=n, percCens=percCens, sdG=sdG )
#       # Sanity check to see if the data simulating function is correct.
#       # all( apply( dfData, 1, function(x){1*(x['logT'] == x['logX']) == x['delta']}) )
#       # all( apply( dfData, 1, function(x){ 1-is.na( x['G'] ) == x['R']}) )
#     }
#
#     # Because for the missing data setting the missing number is unfixed, need to re-calculate this.
#     n2 <- sum( is.na( dfData$G ) )
#     n1 <- n - n2
#
#     matBeta1[i, ] <- gehan.fit1( y=dfData$logX, delta=dfData$delta, matX=as.matrix( dfData[, 1:p] ), wt=wt0, n1=n1, beta.ini=b )
#     matBeta2[i, ] <- gehan.fit2( y=dfData$logX, delta=dfData$delta, matX=as.matrix( dfData[, 1:p] ), wt=wt0, n1=n1, n2=n2, beta.ini=b )
#
#     if ( betaNan ) {
#       wtIp <- 1 / glm( dfData$R ~ as.matrix( cbind( dfData[, 1:(p-1)], dfData$logX, dfData$delta ) ), family=binomial() )$fitted.values
#       matBetaNan[i, ] <- gehan.fit1( y=dfData$logX, delta=dfData$delta,
#                                      matX=as.matrix( dfData[, 1:p] ), wt=wtIp,
#                                      n1=n1, beta.ini=b )
#       # matBetaNan2[i, ] <- gehan.fit1( y=dfData$logX, delta=dfData$delta, matX=as.matrix( dfData[, 1:p] ), wt=1 / dfData$probMiss, n1=n1, beta.ini=b )
#     }
#
#     if ( perturb ) {
#       matWt <- matrix( rexp( n*B ), nrow=n, ncol=B )
#       matBeta1Pt <- perturbfn1( y=dfData$logX, delta=dfData$delta, matX=as.matrix( dfData[, 1:p] ), n1=n1, matWt=matWt, B=B, beta.ini=b )
#       matBeta2Pt <- perturbfn2( y=dfData$logX, delta=dfData$delta, matX=as.matrix( dfData[, 1:p] ), n1=n1, n2=n2, matWt=matWt, B=B, beta.ini=b )
#       matBetaPt <- cbind( matBeta1Pt, matBeta2Pt )
#       matSigmaPt <- cov( matBetaPt )
#       matSigmaA <- matSigmaPt[1:p, 1:p] + matSigmaPt[(p+1):(2*p), (p+1):(2*p)] - matSigmaPt[1:p, (p+1):(2*p)] - matSigmaPt[(p+1):(2*p), 1:p]
#       matSigmaB <- matSigmaPt[1:p, (p+1):(2*p)] - matSigmaPt[(p+1):(2*p), (p+1):(2*p)]
#
#       if ( i >= 2 ) {
#         # Calculate the estimated variance matrix of betaStar from perturbation. The variance comes from current replicate and the weight comes from the previous run (matW hasn't been updated yet).
#         lMatSigmaStarPt[[i - 1]] <- matW %*% matSigmaA %*% t( matW ) + matW %*% matSigmaB + t( matW %*% matSigmaB ) + matSigmaPt[(p+1):(2*p), (p+1):(2*p)]
#       }
#
#       matW <- - matSigmaB %*% solve( matSigmaA )
#       matBetaStar[i, ] <- matW %*% matBeta1[i, ] + ( diag( 1, p ) - matW ) %*% matBeta2[i, ]
#       lMatSigmaPt[[i]] <- matSigmaPt
#     }
#
#
#     if ( track ) {
#       lDfData[[i]] <- dfData
#       if ( perturb ) lMatWt[[i]] <- matWt
#     }
#   }
#
#   to.return <- list( matBeta1=matBeta1, matBeta2=matBeta2, matBetaStar=matBetaStar, lMatSigmaPt=lMatSigmaPt, lMatSigmaStarPt=lMatSigmaStarPt, matBetaNan=matBetaNan, matBetaNan2=matBetaNan2 )
#   if( track ) to.return <- c( to.return, list( lDfData=lDfData, lMatWt=lMatWt ) )
#   return( to.return )
# }
#
# # "expression" type legend labels: c(expression(paste(hat(beta)^(1))),expression(paste(hat(beta)^(2))),expression(paste(hat(beta)^("*"))))
# boxplotBeta <- function( lMatBeta, lBetaTrue, nameVar=NULL, nameEstr=NULL, ylab='', legendLabel=nameEstr, colPalette=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), title='' ) {
#
#   if ( class( lMatBeta ) != 'list' ) lMatBeta <- list( lMatBeta )
#   if ( class( lBetaTrue ) != 'list' ) lBetaTrue <- lapply( lMatBeta, function( x ) lBetaTrue ) # Make sure that the true beta values is a list of equal length as the list of beta samples, even though different beta list components share the same true values.
#   p <- ncol( lMatBeta[[1]] )
#   nEstr <- length( lMatBeta )
#   if ( length( nameVar ) == 0 ) nameVar <- as.character( 1:p )
#   if ( length( nameEstr ) == 0 ) nameEstr <- as.character( 1:nEstr )
#   if ( length( legendLabel ) == 0 ) legendLabel <- nameEstr
#
#   lDfBetaLong <- lapply( 1:nEstr, function( i ) {
#     dfBeta <- data.frame( lMatBeta[[i]], check.names=F )
#     colnames( dfBeta ) <- nameVar
#     betaTrue <- lBetaTrue[[i]]
#     names( betaTrue ) <- nameVar
#     dfBeta$Estimator <- rep( nameEstr[i], length=nrow( lMatBeta[[i]] ) )
#     dfBetaLong <- gather_( data=dfBeta, key_col='Variable', value_col='Estimate', gather_cols=nameVar )
#     dfBetaLong$MSE <- coefSummary( lMatBeta[[i]], betaTrue )[dfBetaLong$Variable, 5]
#     dfBetaLong$betaTrue <- betaTrue[dfBetaLong$Variable]
#     return( dfBetaLong )
#   } )
#
#   dfToplot <- Reduce( rbind, lDfBetaLong )
#   dfToplot$Variable <- factor( dfToplot$Variable, levels=nameVar )
#   dfToplot$Estimator <- factor( dfToplot$Estimator, levels=nameEstr )
#   ggplot <- ggplot( aes( y=Estimate, x=Variable, fill=Estimator ), data=dfToplot ) + geom_boxplot( positon="dodge" ) + geom_errorbar( aes( y=betaTrue, ymax=betaTrue+MSE, ymin=betaTrue-MSE, color='red'), width=0.25, position=position_dodge(width=0.75) ) + geom_point( aes(y=betaTrue, color='red' ), shape=18, size=4, position=position_dodge(width=0.75) ) + scale_fill_manual( values=colPalette, name='Estimator', guide=guide_legend(override.aes = list(shape = NA)), labels=legendLabel ) + xlab('Coviariate') + ylab( ylab ) + scale_color_manual( values='red', name='Error Bar', labels=expression( paste( 'True Value' %+-% 'MSE' ) ) ) + theme(legend.text.align=0)
#   return( ggplot )
# }
#
# trendplotSigma <- function( lMatSigma, bG, nameVar=NULL, nameEstr=NULL, legendLabel=nameEstr, colPalette=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), title='' ){
#
#   if ( class( lMatSigma ) != 'list' ) lMatSigma <- list( lMatSigma )
#   p <- ncol( lMatSigma[[1]] )
#   nEstr <- length( lMatSigma )
#   if ( length( nameVar ) == 0 ) nameVar <- as.character( 1:p )
#   if ( length( nameEstr ) == 0 ) nameEstr <- as.character( 1:nEstr )
#   if ( length( legendLabel ) == 0 ) legendLabel <- nameEstr
#
#   lDfSigmaLong <- lapply( 1:nEstr, function( i ) {
#     dfSigma <- data.frame( lMatSigma[[i]], check.names=F )
#     colnames( dfSigma ) <- nameVar
#     dfSigma$Estimator <- rep( nameEstr[i], length=nrow( lMatSigma[[i]] ) )
#     dfSigma$bG <- bG
#     dfSigmaLong <- gather_( data=dfSigma, key_col='Variable', value_col='Estimate', gather_cols=nameVar )
#     return( dfSigmaLong )
#   } )
#
#   dfToplot <- Reduce( rbind, lDfSigmaLong )
#   dfToplot$Variable <- factor( dfToplot$Variable, levels=nameVar )
#   dfToplot$Estimator <- factor( dfToplot$Estimator, levels=nameEstr )
#   ggplot <- ggplot( aes( y=Estimate, x=bG, color=Estimator ), data=dfToplot ) + geom_line( ) + geom_point( size=3 ) + scale_color_manual( values=colPalette, name='Estimator', labels=legendLabel ) + facet_wrap( ~ Variable ) + ggtitle( title )
#   return( ggplot )
# }
#
#
#
# gehan.fit2=function(y, delta, matX, wt, n1, n2, beta.ini=NULL)
# {
#   # matX is the covariate matrix. The missing values in matX must only be contained in the last column.
#   p <- ncol( matX )
#   n <- n1 + n2
#
#   # Set initial values to be fitted from a lm if not given.
#   if( length(beta.ini)==0 ) {
#     beta1=lm( y ~ matX )$coef[-1]
#     beta2=lm( y ~ matX, subset=(delta==1) )$coef[-1]
#     beta.ini=(beta1+beta2)/2
#   }
#
#   alpha <- lm( matX[1:n1, p] ~ matX[1:n1, 1:(p-1)], weights=wt[1:n1] )$coef[-1] # Fitted alpha value from fitting Z and S on G in the first dataset.
#   matX[(n1+1):n, p] <- matX[(n1+1):n, 1:(p-1)] %*% alpha #impute G for D2
#
#   fit=optim(beta.ini, gehan.obj2, y=y, delta=delta, matX=matX, wt=wt, n1=n1, n2=n2)
#   beta=fit$par
#   return(beta)
# }
# truedata_for_MI <- function( a, b, n1, n2, percCens, sdG, sameErrDist=T, disteG='normal' )
# {
#   n <- n1 + n2
#   p <- length( b )
#   matZ <- mvrnorm( n, mu=rep( 0, length( a ) ) )
#   matG <- matZ %*%
#   matZ <- mvrnorm( n, mu=rep( 0, p-2 ), Sigma=diag( rep( 1, p-2 ) ) )
#   if ( sameErrDist ) {
#     if ( disteG == 'normal' ) G <- 3 + matZ %*% a[1:(p-2)] + S*a[p-1] +
#         c( rnorm(n1, 0, sdG[1]), rnorm( n2, 0, sdG[2] ) )
#     if ( disteG == 'exp' ) G <- 3 + matZ %*% a[1:(p-2)] + S*a[p-1] +
#         c( rexp(n1, rate=1/sdG[1]) - 1/sdG[1], rexp(n2, rate=1/sdG[2]) - 1/sdG[2] )
#   } else {
#     G <- 3 + matZ %*% a[1:(p-2)] + S*a[p-1] +
#       c( rnorm(n1, 0, sdG[1] ), rexp(n2, rate=1/sdG[2]) - 1/sdG[2] )
#   }
#   if ( disteG == 'normal' ) G <- 3 + matZ %*% a[1:(p-2)] + S*a[p-1] + rnorm(n, 0, sdG)
#   if ( disteG == 'exp' ) G <- 3 + matZ %*% a[1:(p-2)] + S*a[p-1] + rexp(n, rate=1/sdG) - 1/sdG
#   if ( sameErrDist ) logT <- matZ %*% b[1:(p-2)] + S*b[p-1] + G*b[p] + rgev(n)
#   else logT <- matZ %*% b[1:(p-2)] + S*b[p-1] + G*b[p] + c( rgev(n1, scale=1), rnorm(n2, mean=10, sd=1.5) )
#   muCens <- 10
#   sdCens <- 1
#   obsCens <- 0
#   while( obsCens < percCens ) {
#     # Keep changing censoring dist'n until the censoring percentage gets close to the disired percentage.
#     muCens <- muCens - 0.25
#     logC <- rnorm( n, muCens, sdCens )
#     delta <- 1 * ( logT <= logC )
#     obsCens <- 1 - sum( delta ) / n
#   }
#   logX <- pmin( logT, logC )
#
#   # Missingness
#   GAll <- G
#   G[(n1+1):n] <- NA
#   R <- 1 - is.na( G )
#   probMiss <- rep( 1, n ) # This isn't a missing data problem, but still include R and probMisss to fit into Nan et al.'s framework.
#   dfData <- data.frame( cbind( matZ, S, G, logX, delta, logT, logC, GAll, R, probMiss ) )
#   colnames( dfData ) <- c( paste0( 'Z', 1:(p-2) ), 'S', 'G', 'logX', 'delta', 'logT', 'logC', 'GAll', 'R', 'probMiss' )
#   return( dfData )
# }
#
# truedata2 <- function( a1, a2, b, n1, n2, percCens, sdG, sameErrDist=T, disteG='normal' )
# {
#   n <- n1 + n2
#   p <- length( b )
#   S <- rnorm( n, 0, 1)
#   matZ <- mvrnorm( n, mu=rep( 0, p-3 ), Sigma=diag( rep( 1, p-3 ) ) )
#   G1 <- 3 + matZ %*% a1[1:(p-3)] + S*a1[p-2] + rnorm(n, 0, sdG)
#   G2 <- 3 + matZ %*% a2[1:(p-3)] + S*a2[p-2] + rexp(n, rate=1/sdG) - 1/sdG
#   # if ( disteG == 'exp' ) G <- 3 + matZ %*% a[1:(p-2)] + S*a[p-1] + rexp(n, rate=1/sdG) - 1/sdG
#   if ( sameErrDist ) logT <- matZ %*% b[1:(p-3)] + S*b[p-2] + G1*b[p-1] + G2*b[p] + rgev(n)
#   else logT <- matZ %*% b[1:(p-3)] + S*b[p-1] + G*b[p] + c( rgev(n1, scale=1), rnorm(n2, mean=10, sd=1.5) )
#   muCens <- 10
#   sdCens <- 1
#   obsCens <- 0
#   while( obsCens < percCens ) {
#     # Keep changing censoring dist'n until the censoring percentage gets close to the disired percentage.
#     muCens <- muCens - 0.25
#     logC <- rnorm( n, muCens, sdCens )
#     delta <- 1 * ( logT <= logC )
#     obsCens <- 1 - sum( delta ) / n
#   }
#   logX <- pmin( logT, logC )
#
#   # Missingness
#   GAll <- cbind( G1, G2 )
#   G1[(n1+1):n] <- NA
#   G2[(n1+1):n] <- NA
#   R <- 1 - is.na( G1 )
#   probMiss <- rep( 1, n ) # This isn't a missing data problem, but still include R and probMisss to fit into Nan et al.'s framework.
#   dfData <- data.frame( cbind( matZ, S, G1, G2, logX, delta, logT, logC, GAll, R, probMiss ) )
#   colnames( dfData ) <- c( paste0( 'Z', 1:(p-3) ), 'S', 'G1', 'G2', 'logX', 'delta', 'logT', 'logC', 'GAll1', 'GAll2', 'R', 'probMiss' )
#   return( dfData )
# }

# truedataMiss <- function( a, b, eta, n, percCens, sdG )
# {
#   # Return a simulated data matrix, first part being the first dataset and second part being the second dataset. Note the order of the variables (columns) are Z - S - G - others.
#   p <- length( b )
#   S <- rnorm( n, 0, 1)
#   matZ <- mvrnorm( n, mu=rep( 0, p-2 ), Sigma=diag( rep( 1, p-2 ) ) )
#   G <- 3 + matZ %*% a[1:(p-2)] + S*a[p-1] + rnorm(n, 0, sdG)
#   logT <- matZ %*% b[1:(p-2)] + S*b[p-1] + G*b[p] + rgev(n)
#
#   # Ceonsoring
#   muCens <- 10
#   sdCens <- 1
#   obsCens <- 0
#   while( obsCens < percCens )
#   {
#     # Keep changing censoring dist'n until the censoring percentage gets close to the disired percentage.
#     muCens <- muCens - 0.25
#     logC <- rnorm( n, muCens, sdCens )
#     delta <- 1 * ( logT <= logC )
#     obsCens <- 1 - sum( delta ) / n
#   }
#   logX <- pmin( logT, logC )
#
#   # Missingness
#   probMiss <- exp( matZ %*% eta[1:(p-2)] + S*eta[p-1] + G*eta[p] + eta[p+1]*logX + eta[p+2]*delta ) / ( 1 + exp( matZ %*% eta[1:(p-2)] + S*eta[p-1] + G*eta[p] + eta[p+1]*logX + eta[p+2]*delta ) )
#   R <- rbinom( n=n, size=1, prob=probMiss )
#   GAll <- G
#   G[R==0] <- NA
#
#   dfData <- data.frame( cbind( matZ, S, G, logX, delta, logT, logC, GAll, R, probMiss ) )
#   colnames( dfData ) <- c( paste0( 'Z', 1:(p-2) ), 'S', 'G', 'logX', 'delta', 'logT', 'logC', 'GAll', 'R', 'probMiss' )
#   dfData <- rbind( subset( dfData, R==1 ), subset( dfData, R==0 ) )
#   return( dfData )
# }

# runall <- function( n1=150, n2=150, a=c(0.1, 0.1), b=c(1, 1, 1), perc.cens=0.1, sd.G=1, B=500, rep=500, quiet=F, track=F ) {
#   p <- length( b )
#   beta1s <- matrix( NA, rep, p )
#   beta2s <- matrix( NA, rep, p )
#   beta.stars <- matrix( NA, rep, p )
#   sigmas <- matrix( NA, rep, 2*p )
#   if( track ) {
#     datasets <- list( )
#     weights <- list( )
#   }
#   for( i in 1:rep ) {
#     if( !quiet & (i %% 20 == 0) ) cat( round( i / rep * 100, 0 ), '% ' )
#     n <- n1 + n2
#     data=truedata(a=a, b=b, n1=n1, n2=n2, perc.cens=perc.cens, sd.G=sd.G)
#     wt=matrix(rexp(n*B),nrow=n,ncol=B)
#     if ( track ) {
#       datasets[[i]] <- data
#       weights[[i]] <- wt
#     }
#     beta1=try( perturbfn1(y=data$logX, delta=data$delta, x=cbind(data$Z1,data$S,data$G), n1=n1, wt=wt, B=B, beta.ini=b ), silent=T )
#     beta2=try( perturbfn2(y=data$logX, delta=data$delta, x=cbind(data$Z1,data$S,data$G), n1=n1, n2=n2, wt=wt, B=B, beta.ini=b ), silent=T )
#     if( ( 'try-error' %in% class( beta1 ) ) | ( 'try-error' %in% class( beta2 ) ) ) return( list( data=data, wt=wt ) )
#     beta1.0 <- beta1[[1]]
#     beta2.0 <- beta2[[1]]
#     beta1s[i, ] <- beta1.0
#     beta2s[i, ] <- beta2.0
#     beta.perturb=rbind(beta1[[2]],beta2[[2]])
#     sigma=cov(t(beta.perturb))
#     sigmas[i, ] <- diag( sigma )
#     sigma.a <- sigma[1:p, 1:p] + sigma[(p+1):(2*p), (p+1):(2*p)] - sigma[1:p, (p+1):(2*p)] - sigma[(p+1):(2*p), 1:p]
#     sigma.b <- sigma[1:p, (p+1):(2*p)] - sigma[(p+1):(2*p), (p+1):(2*p)]
#     A <- - sigma.b %*% solve( sigma.a )
#     beta.stars[i, ] <- A %*% beta1.0 + ( diag( 1, p ) - A ) %*% beta2.0
#   }
#   to.return <- list( beta1=beta1s, beta2=beta2s, beta.star=beta.stars, sigma=sigmas )
#   if( track ) to.return <- c( to.return, list( data=datasets, wt=weights ) )
#   return( to.return )
# }

# truedata.mis=function(a,b,n1,n2,perc.cens,sd.G)
# {
#   while( T ) {
#     n = n1+n2
#     S = rbinom(n,1,.1)
#     Z1 = rnorm(n,0,1)
#     G = a[1]*Z1^2 + a[2]*S + rnorm(n,0,sd.G)
#     logT = b[1]*Z1 + b[2]*S + b[3]*G +rgev(n)
#     mu.cens=10
#     sd.cens=1
#     obs.cens=0
#     while(obs.cens<perc.cens)
#     {
#       mu.cens=mu.cens-0.25
#       logC=rnorm(n,mu.cens,sd.cens)
#       delta = 1*(logT<=logC)
#       obs.cens=1-sum(delta)/n
#     }
#     logX = pmin(logT,logC)
#     G[(n1+1):(n1+n2)]=NA
#     data=cbind(Z1,S,G,logX,delta)
#     if( sum( S[1:n1] * delta[1:n1] ) >= 1 )return(as.data.frame(data))
#   }
# }

  # runall.noperturb.mis <- function( n1=150, n2=150, a=c(0.1, 0.1), b=c(1, 1, 1), perc.cens=0.1, sd.G=1, rep=500, gamma=F, quiet=T, track=F ) {
  #   p <- length( b )
  #   wt <- rep(1, n1+n2)
  #   beta1s <- matrix( NA, rep, p )
  #   beta2s <- matrix( NA, rep, p )
  #   if ( gamma ) {
  #     gamma1.gehans <- matrix( NA, rep, p - 1 )
  #     gamma1.plugins <- matrix( NA, rep, p - 1 )
  #     gamma2.gehans <- matrix( NA, rep, p - 1 )
  #     gamma2.plugins <- matrix( NA, rep, p - 1 )
  #   }
  #   datasets <- list( )
  #   for( i in 1:rep ) {
  #     if( !quiet & (i %% 100 == 0) ) cat( round( i / rep * 100, 0 ), '% ' )
  #     n <- n1 + n2
  #     data=truedata.mis(a=a, b=b, n1=n1, n2=n2, perc.cens=perc.cens, sd.G=sd.G)
  #     beta1=gehan.fit1( y=data$logX,delta=data$delta,x=cbind(data$Z1,data$S,data$G),wt=wt,n1=n1, beta.ini=b )
  #     alpha=lm( G ~ Z1 + S, data=data[1:n1, ] )$coef[-1]
  #     beta2=gehan.fit2( y=data$logX, delta=data$delta, x=cbind(data$Z1,data$S,c( data$G[1:n1], cbind( data$Z1, data$S )[(n1+1):n, ] %*% alpha ) ), wt=wt, n1=n1, n2=n2, beta.ini=b )
  #     beta1s[i, ] <- beta1
  #     beta2s[i, ] <- beta2
  #     if( gamma ) {
  #       gamma1.gehan=gehan.fit1( y=data$logX,delta=data$delta,x=cbind(data$Z1,data$S),wt=wt,n1=n1, beta.ini=b[1:(p-1)]+a*b[p] )
  #       gamma1.plugin=beta1[1:(p-1)] + alpha*beta1[p]
  #       gamma2.gehan=gehan.fit1( y=data$logX,delta=data$delta,x=cbind(data$Z1,data$S),wt=wt, n1=n1+n2, beta.ini=b[1:(p-1)]+a*b[p] )
  #       gamma2.plugin=beta2[1:(p-1)] + alpha*beta2[p]
  #       gamma1.gehans[i, ] <- gamma1.gehan
  #       gamma1.plugins[i, ] <- gamma1.plugin
  #       gamma2.gehans[i, ] <- gamma2.gehan
  #       gamma2.plugins[i, ] <- gamma2.plugin
  #     }
  #     if ( track ) datasets[[i]] <- data
  #   }
  #   to.return <- list( beta1=beta1s, beta2=beta2s )
  #   if( gamma ) to.return <- c( to.return, list( gamma1.gehan=gamma1.gehans, gamma1.plugin=gamma1.plugins, gamma2.gehan=gamma2.gehans, gamma2.plugin=gamma2.plugins ) )
  #   if( track ) to.return <- c( to.return, datasets= )
  #   if( track ) return( list( beta1=beta1s, beta2=beta2s, gamma1.gehan=gamma1.gehans, gamma1.plugin=gamma1.plugins, gamma2.gehan=gamma2.gehans, gamma2.plugin=gamma2.plugins ), data=datasets )
  #   return( to.return )
  # }

#   # Function to get beta^1, beta^2, and four gamma estimates, without perturbation.
#   runall.noperturb <- function( n1=150, n2=150, a=c(0.1, 0.1), b=c(1, 1, 1), perc.cens=0.1, sd.G=1, rep=500, gamma=F, quiet=T, track=F ) {
#     p <- length( b )
#     wt <- rep(1, n1+n2)
#     beta1s <- matrix( NA, rep, p )
#     beta2s <- matrix( NA, rep, p )
#     if ( gamma ) {
#       gamma1.gehans <- matrix( NA, rep, p - 1 )
#       gamma1.plugins <- matrix( NA, rep, p - 1 )
#       gamma2.gehans <- matrix( NA, rep, p - 1 )
#       gamma2.plugins <- matrix( NA, rep, p - 1 )
#     }
#     datasets <- list( )
#     for( i in 1:rep ) {
#       if( !quiet & (i %% 100 == 0) ) cat( round( i / rep * 100, 0 ), '% ' )
#       n <- n1 + n2
#       data=truedata(a=a, b=b, n1=n1, n2=n2, perc.cens=perc.cens, sd.G=sd.G)
#       beta1=gehan.fit1( y=data$logX,delta=data$delta,x=cbind(data$Z1,data$S,data$G),wt=wt,n1=n1, beta.ini=b )
#       alpha=lm( G ~ Z1 + S, data=data[1:n1, ] )$coef[-1]
#       beta2=gehan.fit2( y=data$logX, delta=data$delta, x=cbind(data$Z1,data$S,c( data$G[1:n1], cbind( data$Z1, data$S )[(n1+1):n, ] %*% alpha ) ), wt=wt, n1=n1, n2=n2, beta.ini=b )
#       beta1s[i, ] <- beta1
#       beta2s[i, ] <- beta2
#       if( gamma ) {
#         gamma1.gehan=gehan.fit1( y=data$logX,delta=data$delta,x=cbind(data$Z1,data$S),wt=wt,n1=n1, beta.ini=b[1:(p-1)]+a*b[p] )
#         gamma1.plugin=beta1[1:(p-1)] + alpha*beta1[p]
#         gamma2.gehan=gehan.fit1( y=data$logX,delta=data$delta,x=cbind(data$Z1,data$S),wt=wt, n1=n1+n2, beta.ini=b[1:(p-1)]+a*b[p] )
#         gamma2.plugin=beta2[1:(p-1)] + alpha*beta2[p]
#         gamma1.gehans[i, ] <- gamma1.gehan
#         gamma1.plugins[i, ] <- gamma1.plugin
#         gamma2.gehans[i, ] <- gamma2.gehan
#         gamma2.plugins[i, ] <- gamma2.plugin
#       }
#       if ( track ) datasets[[i]] <- data
#     }
#     to.return <- list( beta1=beta1s, beta2=beta2s )
#     if( gamma ) to.return <- c( to.return, list( gamma1.gehan=gamma1.gehans, gamma1.plugin=gamma1.plugins, gamma2.gehan=gamma2.gehans, gamma2.plugin=gamma2.plugins ) )
#     if( track ) to.return <- c( to.return, datasets= )
#     if( track ) return( list( beta1=beta1s, beta2=beta2s, gamma1.gehan=gamma1.gehans, gamma1.plugin=gamma1.plugins, gamma2.gehan=gamma2.gehans, gamma2.plugin=gamma2.plugins ), data=datasets )
#     return( to.return )
#   }
#
# boxplot.sigma <- function( sigma.perturb, sigma.ep.beta1, sigma.ep.beta2, var.names=NULL, estimator.names=paste0( 'beta^(', 1:2, ')' ), ylab='', palette=cbbPalette ) {
#   # sigma.perturb should be a rep * 2p matrix, recording variance estimates from perturbations from all replicates and for beta_(1) and beta_(2)
#   # sigma.ep.beta1 and sigma.ep.beta2 should each be a p-dimensional vector, recording the estimated empirical variance acoss all replicates for beta_(1) and beta_(2)
#   p <- length( sigma.ep.beta1 )
#   rep <- nrow( sigma.perturb )
#   if ( length( var.names ) == 0 ) var.names=as.character( 1:p )
#   toplot <- data.frame( sigma.perturb=as.vector( sigma.perturb ), group.estimator=factor( rep( estimator.names, each=rep*p ), levels=estimator.names, ordered=T ), group.var=factor( rep( rep( var.names, each=rep ), times=ncol( sigma.perturb ) / p ), levels=var.names, ordered=T ), sigma.ep=c( rep( sigma.ep.beta1, each=rep ), rep( sigma.ep.beta2, each=rep ) ) )
#   plot <- ggplot(aes(y = sigma.perturb, x = group.var, fill = group.estimator), data = toplot) + geom_boxplot() + geom_point( aes(y=sigma.ep), shape=18, size=4, color='red', position=position_dodge(width=0.75) ) + scale_fill_manual(values=palette, name='Estimator') + xlab('Covariate') + ylab( ylab )
#   print( plot )
# }
# trendplot.sigma <- function( lMatSigma, varValues, bG, var.varying.name=NULL, colPalette=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), xlab=F, title ){
#   # sigmas should be a length( k ) list m * r column matrices.
#   # r is the number of estimators.
#   # k is the number of different varying variable values.
#   # m is the number of differnet bGs used.
#   # var.varying.values is length k.
#   # bGs is length m.
#   r <- ncol( sigmas[[1]] )
#   m <- nrow( sigmas[[1]] )
#   toplot <- data.frame( )
#   for( i in 1:length( sigmas ) ) {
#     toplot.tmp <- data.frame( sigma=as.vector( sigmas[[i]] ), var.varying.values=factor( rep( var.varying.values[i], r*m ), levels=var.varying.values, ordered=T ), bGs=rep( bGs, r ), group.estimator=factor( rep( estimator.names, each=m ), levels=estimator.names, ordered=T ) )
#     toplot <- rbind( toplot, toplot.tmp )
#   }
#   plot <- ggplot(aes(y=sigma, x=bGs, color=group.estimator ), data = toplot) + geom_line( ) + geom_point( size=3 ) + scale_color_manual(values=palette, name='Estimator', labels=labels ) + facet_wrap( ~ var.varying.values ) + ylab( ylab ) + ggtitle( title )
#   if ( xlab ) { plot <- plot + xlab( expression( paste( 'True ', beta[G] ) ) ) } else {
#     plot <- plot + theme( axis.title.x=element_blank() )
#   }
#   print( plot )
# }
}# gehan.obj1<-function( beta, y, delta, matX, wt, n1)
# {
#   # just use D1
#   # This also works for IPW setting where the wt is the given IPW.
#   wt <- wt[1:n1]
#   delta <- delta[1:n1]
#   y <- y[1:n1]
#   matX <- matX[1:n1,]
#
#   residual <- y - as.vector( matX %*% beta )
#   index <- order( residual )
#   residual.order <- residual[index]
#   delta.order <- delta[index]
#   wt.order <- wt[index]
#
#   s1 <- cumsum(delta.order*wt.order)
#   s2 <- delta.order*rev(cumsum(rev(wt.order)))
#   return(sum((s1-s2)*residual.order*wt.order))
# }
#
# gehan.obj2<-function( beta, y, delta, matX, wt, n1, n2)
# {
#   #D1
#   wt1 <- wt[1:n1]
#   delta1 <- delta[1:n1]
#   y1 <- y[1:n1]
#   matX1 <- matX[1:n1,]
#
#   residual <- y1 - as.vector( matX1 %*% beta )
#   index <- order(residual)
#   residual.order <- residual[index]
#   delta.order <- delta1[index]
#   wt.order <- wt1[index]
#
#   s1 <- cumsum(delta.order*wt.order)
#   s2 <- delta.order*rev(cumsum(rev(wt.order)))
#   firstpart <- sum((s1-s2)*residual.order*wt.order)
#
#   #D2
#   wt2 <- wt[(n1+1):(n1+n2)]
#   delta2 <- delta[(n1+1):(n1+n2)]
#   y2 <- y[(n1+1):(n1+n2)]
#   matX2 <- matX[(n1+1):(n1+n2),]
#
#   residual <- y2-as.vector( matX2 %*% beta )
#   index <- order(residual)
#   residual.order <- residual[index]
#   delta.order <- delta2[index]
#   wt.order <- wt2[index]
#
#   s1 <- cumsum(delta.order*wt.order)
#   s2 <- delta.order*rev(cumsum(rev(wt.order)))
#   secondpart <- sum((s1-s2)*residual.order*wt.order)
#
#   #add together
#   return(firstpart+secondpart)
# }