# MCMC fit for a multi-class logistic regression model
# author: Ying Wang
# version: 1.0

ls()
rm(list=ls())
gc()
library(statip)
library(nnet)
library(sm)
#======   define functions   ======
sigmaPrior = 0.1
logPrior = function(beta, uPrior) {
  logp = 0
  for(j in 1 : nrow(uPrior)) {
    for(k in 1 : ncol(uPrior)) {
      logp = logp + log(dnorm(beta[j, k], uPrior[j,k], sigmaPrior))
    }
  }
  return(logp)
}

logLikelihood = function(beta, x, y) {
  likelihood = 0.0
 
    for (h in 1 : nrow(x)) {
      z1 = sum(beta[1,] * x[h,])
      z2 = sum(beta[2,] * x[h,])
      z3 = sum(beta[3,] * x[h,])
      z4 = sum(beta[4,] * x[h,])
      ezsum = exp(z2) + exp(z1) + exp(z3)+exp(z4)
      p1 = exp(z1) / ezsum
      p2 = exp(z2) / ezsum
      p3 = exp(z3) / ezsum
      if(y[h]==1){
        likelihoodSample = log(p1)
      } else if(y[h]==2){
        likelihoodSample = log(p2)
      }else if(y[h]==3){
        likelihoodSample = log(p3)
      }else {
        likelihoodSample = log(p4)
      }

      likelihood = likelihood + likelihoodSample
    }
 
  return(likelihood)
}

#======   prepare data   ======

  xCols = 4
  yCols = 4
  betaTrue = array(0, dim=c(yCols, xCols))
  uPrior = array(0, dim=c(yCols, xCols)) 
  uPrior[1,1] = 1.5
  uPrior[1,2] = 2.5
  uPrior[1,3] = 3.5
  uPrior[1,4] = 4.5
  uPrior[2,1] = 0.5
  uPrior[2,2] = 0.5
  uPrior[2,3] = 0.5
  uPrior[2,4] = 0.5
  uPrior[3,1] = 1.5
  uPrior[3,2] = 2.5
  uPrior[3,3] = 3.5
  uPrior[3,4] = 4.5
  uPrior[4,1] = 4.5
  uPrior[4,2] = 1.5
  uPrior[4,3] = 3.5
  uPrior[4,4] = 4.5
  for (j in 1 : yCols) {
    for (k in 1 : xCols) {
      betaTrue[j, k] = uPrior[j,k]#rnorm(1, uPrior[j,k], sigmaPrior)
    }
  }

  nRecord = 1000
  x = array(1, dim=c(nRecord, xCols)) 
  y = array(0, dim=c(nRecord, 1)) 
  for (k in 2 : xCols) {
    x[,k] = runif(nRecord, -1, 1)
  }
  for (h in 1 : nRecord) {
    
    # bb = betaTrue[1, ] - betaTrue[2, ]
    # zz = sum(bb * x[h,])
    # p2 = exp(zz)/(1+exp(zz))
    # y[h] =rbern(1,p2)
    
    z1 = sum(betaTrue[1,] * x[h,])
    z2 = sum(betaTrue[2,] * x[h,])
    z3 = sum(betaTrue[3,] * x[h,])
    z4 = sum(betaTrue[4,] * x[h,])
    ezsum = exp(z1) + exp(z2) +exp(z3) + exp(z4)
    p1 = exp(z1) / ezsum
    p2 = exp(z2) / ezsum
    p3 = exp(z3) / ezsum
    p4 = exp(z4) / ezsum

    uu = runif(1, min = 0, max = 1)
    if(uu < p1){
      y[h] = 1
    } else if(uu >=p1&&uu<(p1+p2)) {
      y[h] = 2
    } else if(uu >=(p1+p2)&&uu<(p1+p2+p3)) {
    y[h] = 3
  } else {
      y[h] = 4
    }
    
    
  }


#======   sampling beta   ======
nSamples = 5000 # the length of the MCMC samples
sigmaProposal = 0.05#0.05;
sigmaProposalPost = 0.05
betaCurr = array(0, dim=c(yCols, xCols)) 
betaNew = array(0, dim=c(yCols, xCols)) 
betaPosCurr = array(0, dim=c(yCols, xCols)) 
betaPosNew = array(0, dim=c(yCols, xCols))
init = seq(0.5, xCols)/1.0
for (j in 1 : yCols) {
  betaCurr[j,] = init
  betaNew[j,] = init
  betaPosCurr[j,] = init
  betaPosNew[j,] = init
}

betaSamplesLikelihood = array(0, dim=c(nSamples, yCols, xCols)) 
betaSamplesPosterior = array(0, dim=c(nSamples, yCols, xCols)) 

for (i in 1 : nSamples) {
  if (i %% 50 ==0) {
    print(i)
  }

 
    for (k in 1 : xCols) {
      #sample liklihood
      dd1 = sigmaProposal * rnorm(1, mean = 0, sd = 1)
      betaNew[1, k] = betaCurr[1, k] + dd1
      betaNew[2, k] = betaCurr[2, k] #+ sigmaProposal * rnorm(1, mean = 0, sd = 1)
      dd3 = sigmaProposal * rnorm(1, mean = 0, sd = 1)
      betaNew[3, k] = betaCurr[3, k] + dd3
      dd4 = sigmaProposal * rnorm(1, mean = 0, sd = 1)
      betaNew[4, k] = betaCurr[4, k] + dd3
      logLikelihoodNew = logLikelihood(betaNew, x, y)
      logLikelihoodCurr = logLikelihood(betaCurr, x, y)
      logAlpha = logLikelihoodNew - logLikelihoodCurr
      alpha = exp(logAlpha)
      #print(alpha)
      u = runif(1, min = 0, max = 1);
      if (u < min(1.0, alpha)) {
        betaCurr[1, k] = betaNew[1, k];
        betaCurr[2, k] = betaNew[2, k];
        betaCurr[3, k] = betaNew[3, k];
        betaCurr[4, k] = betaNew[4, k];
      }

    }
  
 # bbSamplesLikelihood[i,] = bbc
  betaSamplesLikelihood[i,,] = betaCurr
  #betaSamplesPosterior[i,,] = betaPosCurr
}



burnin = 500
dirName = "output2/"
for (j in 1 : yCols) {
  for (k in 1 : xCols) {
    beta = betaSamplesLikelihood[, j, k]
    traceTitle = sprintf("Trace of beta%i%i for likelihood", j, k)
    fName = paste(dirName, traceTitle, ".png", sep="")
    png(filename = fName)
    plot(beta, type = "l")
    title(traceTitle)
    dev.off()

  }
}
