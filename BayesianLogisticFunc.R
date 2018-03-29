
sigmaPrior = 0.1
logPrior = function(beta, uPrior) {
  logp = 0

    for(k in 1 : ncol(uPrior)) {
      logp = logp + log(dnorm(beta[k], uPrior[k], sigmaPrior))
    }
  
  return(logp)
}

logLikelihood = function(beta, x, y) {
  likelihood = 0.0
 
    for (h in 1 : nrow(x)) {
      bb = beta
      zz = sum(bb * x[h,])
      pp = exp(zz) / (1.0 + exp(zz))
        yy = y[h]

      likelihoodSample = log(pp) * yy + log(1.0 - pp) * (1.0 - yy)
      likelihood = likelihood + likelihoodSample
    }
 
  return(likelihood)
}

getLikelihoodSample = function(nSamples, y, x) {
  nSamples = 1000 # the length of the MCMC samples
  sigmaProposal = 0.05#0.05;
  yCols = ncol(y)
  if(is.null(yCols)){
    yCols = 1
  }
    
  xCols = ncol(x)
  betaCurr = array(0, dim=c(yCols, xCols)) 
  betaNew = array(0, dim=c(yCols, xCols)) 
  
  betaSamples = array(0, dim=c(nSamples, yCols, xCols)) 

  for (i in 1 : nSamples) {
    if (i %% 50 ==0) {
      print(i)
    }
  
    for (j in 1 : yCols) {
      for (k in 1 : xCols) {
        #sample liklihood
        betaNew[j, k] = betaCurr[j, k] + sigmaProposal * rnorm(1, mean = 0, sd = 1)
        logLikelihoodNew = logLikelihood(betaNew, x, y)
        logLikelihoodCurr = logLikelihood(betaCurr, x, y)
        logAlpha = logLikelihoodNew - logLikelihoodCurr
        alpha = exp(logAlpha)
        #print(alpha)
        u = runif(1, min = 0, max = 1);
        if (u < min(1.0, alpha)) {
          betaCurr[j, k] = betaNew[j, k];
        }
       
      }
    }
    betaSamples[i,,] = betaCurr

  }
  return(betaSamples) 
}

getPosteriorSample = function(nSamples, y, x) {
  nSamples = 1000 # the length of the MCMC samples
  sigmaProposal = 0.05#0.05;
  yCols = ncol(y)
  xCols = ncol(x)
  betaCurr = array(0, dim=c(yCols, xCols)) 
  betaNew = array(0, dim=c(yCols, xCols)) 
  
  betaSamples = array(0, dim=c(nSamples, yCols, xCols)) 
  
  for (i in 1 : nSamples) {
    if (i %% 50 ==0) {
      print(i)
    }
    
    for (j in 1 : yCols) {
      for (k in 1 : xCols) {
        #sample liklihood
        betaNew[j, k] = betaCurr[j, k] + sigmaProposal * rnorm(1, mean = 0, sd = 1)
        priorMean = array(0, dim=c(yCols, xCols)) 
        logPosteriorNew = logLikelihood(betaNew, x, y) + logPrior(betaNew, priorMean)
        logPosteriorCurr = logLikelihood(betaCurr, x, y) + logPrior(betaCurr, priorMean)
        logAlphaPosterior = logPosteriorNew - logPosteriorCurr
        alpha = exp(alphaPosterior)
        #print(alpha)
        u = runif(1, min = 0, max = 1);
        if (u < min(1.0, alpha)) {
          betaCurr[j, k] = betaNew[j, k];
        }
        
      }
    }
    betaSamples[i,,] = betaCurr
    
  }
  return(betaSamples) 
}

estimateModel = function(samples, yCols, xCols) {
  betaEst = rep(0, xCols) 


    for (k in 1 : xCols) {
      beta = samples[,  k]
      hi = hist(beta, breaks = 15, freq = FALSE, main = NULL)
      lines(density(beta))
      idx = which.is.max(hi$density)
      betaEst = hi$mids[idx]
      betaEst[k] = as.double(betaEst)
    }
  
  return(betaEst)
}

predict = function(beta, y, x) {

  N = nrow(x)
  cnt = 0
  for (h in 1 : N) {

      z = exp(sum(beta * x[h,]))
      p = z/(1 + z)
      if(y[h]==1 && p>0.5) {
        cnt = cnt + 1
      }

  }
  accuracy = cnt / N
  return(accuracy)
}

