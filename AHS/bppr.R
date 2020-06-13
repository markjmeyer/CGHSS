###### Bayesian Penalized Probit Regression ######
library(truncnorm)
library(mvtnorm)
library(coda)
library(numbers)

### bppr ###
bppr  <- function(fixed, random, data, B = 10000, burnin = B, 
                  al = 1, bl = 1, up = 100, dots = up/10,
                  cut.start = c(0, 0.5, 1), nu = 0){
  ### data ###
  mf      <- model.frame(formula = fixed, data = data)
  Y       <- model.response(mf)
  X       <- model.matrix(attr(mf,"terms"), data = mf)
  randf   <- as.formula(paste(deparse(random),"- 1"))
  L       <- model.matrix(randf, data = data)

  Vb      <- solve(t(X)%*%X)
  LtL     <- t(L)%*%L

  nL      <- ncol(L)
  nX      <- ncol(X)
  n       <- nrow(X)
  J       <- length(levels(factor(Y)))
  Yl      <- sort(unique(Y))

  ### build matrices ###
  betam   <- matrix(0, nrow = B + burnin, ncol = nX)
  alpham  <- matrix(0, nrow = B + burnin, ncol = nL)
  lambda  <- rep(1, B + burnin)
  eta     <- rep(1, n)
  gammas  <- matrix(cut.start, nrow = B + burnin, ncol = J-1, byrow = TRUE)
  Z       <- rep(0, n)

  ### sampler ###
  for(b in 2:(B + burnin)){
  
    ### sample latent Zs ###
    low   <- c(-Inf, gammas[b-1,])
    high  <- c(gammas[b-1,], Inf)
    muz   <- X%*%betam[b-1,] + L%*%alpham[b-1,]
    for(i in 1:n){
      lev   <- Y[i]
      Z[i]  <- rtruncnorm(1, a = low[lev], b = high[lev], mean = muz[i], sd = 1/eta[i])
    }
  
    ### sample betas ###
    mub       <- Vb%*%t(X)%*%(Z - L%*%alpham[b-1,])
    betam[b,] <- c(rmvnorm(1, mub, Vb))
  
    ### sample alphas ###
    Pa          <- LtL + lambda[b-1]*diag(nL)
    Va          <- solve(Pa)
    mua         <- Va%*%t(L)%*%(Z - X%*%betam[b-1,])
    alpham[b,]  <- c(rmvnorm(1, mua, Va))
  
    ### sample lambda ###
    ata         <- t(alpham[b-1,])%*%alpham[b-1,]
    lambda[b]   <- rgamma(1, nL/2 + al, ata/2 + bl)
  
    ### sample gammas ###
    cuts        <- c(gammas[b-1,], Inf)
    for(j in 2:length(gammas[b-1,])){
      au  <- max(max(Z[Y == Yl[j]]), cuts[j-1])
      bu  <- min(min(Z[Y == Yl[j+1]]), cuts[j+1])
      gammas[b,j] <- runif(1, au, bu)
    }
    
    ### sample etas ###
    if(nu > 0){
      Zhat  <- Z - X%*%betam[b,]
      ae    <- (nu+1)/2
      be    <- 2/(nu + Zhat^2)
      eta   <- rgamma(n, ae, be)
    }

    if(mod(b, dots) == 0){
      cat('.')
    }
    
    if(mod(b, up) == 0){
      cat(paste("\n",b,"samples completed\n"))
    }
    
  }

  colnames(betam)   <- colnames(X)
  colnames(alpham)  <- colnames(L)
  colnames(gammas)  <- paste('g', 1:(J-1), sep = '')
  lambv             <- matrix(lambda, ncol = 1)
  colnames(lambv)   <- 'lambda'
  geweke  <- c(geweke.diag(betam[-c(1:burnin),])$z, geweke.diag(alpham[-c(1:burnin),])$z,
               geweke.diag(gammas[-c(1:burnin),2:(J-1)])$z, geweke.diag(lambv[-c(1:burnin)])$z)
  
  out <- list(betas = betam[-c(1:burnin),], alphas = alpham[-c(1:burnin),],
              lambda = lambda[-c(1:burnin)], gammas = gammas[-c(1:burnin),],
              Z = Z, eta = eta, geweke = geweke)
  
  return(out)
}


### marginal effects ###
me <- function(eff, model, data){
  PPC     <- cbind(model$betas, model$alphas)
  FX      <- model.matrix(~ pop0 + gdp0 + factor(year) + factor(region) + lag0 + lag1 + lag2, 
                          data = data)
  FX0     <- FX1  <- FX2  <- FX
  FX0[,which(colnames(FX) == 'lag0')] <- eff
  FX1[,which(colnames(FX) == 'lag1')] <- eff
  FX2[,which(colnames(FX) == 'lag2')] <- eff
  
  pred0   <- pnorm(FX0%*%t(PPC))
  pred1   <- pnorm(FX1%*%t(PPC))
  pred2   <- pnorm(FX2%*%t(PPC))

  pred0m  <- apply(pred0, 1, median)
  pred1m  <- apply(pred1, 1, median)
  pred2m  <- apply(pred2, 1, median)
  
  mme     <- c(mean(pred0m), mean(pred1m), mean(pred2m))

  pred0l  <- apply(pred0, 1, quantile, probs = 0.025)
  pred1l  <- apply(pred1, 1, quantile, probs = 0.025)
  pred2l  <- apply(pred2, 1, quantile, probs = 0.025)
  
  lme     <- c(mean(pred0l), mean(pred1l), mean(pred2l))
  
  pred0h  <- apply(pred0, 1, quantile, probs = 0.975)
  pred1h  <- apply(pred1, 1, quantile, probs = 0.975)
  pred2h  <- apply(pred2, 1, quantile, probs = 0.975)

  hme     <- c(mean(pred0h), mean(pred1h), mean(pred2h))

  out     <- list(mme = mme, lme = lme, hme = hme)
  
  return(out)
}



