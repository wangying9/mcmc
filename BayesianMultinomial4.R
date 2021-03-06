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
sigmaPrior = 1
logPrior = function(beta, uPrior) {
  logpropAccum = 0
  for(j in 1 : nrow(uPrior)) {
    for(k in 1 : ncol(uPrior)) {
      prop = dnorm(beta[j, k], uPrior[j,k], sigmaPrior)
      logprop = log(prop)
      logpropAccum = logpropAccum + logprop
    }
  }
  return(logpropAccum)
}

logLikelihood = function(beta, x, y) {
  likelihood = 0.0
  L = nrow(beta)
  for (h in 1 : nrow(x)) {
    sumez = 0.0
    for (j in 1 : L) {
      z = sum(beta[j,] * x[h,])
      ez = exp(z)
      sumez = sumez + ez
    }
      z = sum(beta[y[h],] * x[h,])
      ez = exp(z)
      sp = ez / sumez

    likelihood = likelihood + log(sp)
  }
  return(likelihood)
}

#======   prepare data   ======
loadData =  "MyData" 
if (loadData == "MockData") {
  dataDF = read.csv(file="mockData.txt",head=FALSE,sep=",")
  xCols = 4
  yCols = 2
  yLevel = 2 ^ yCols
  data = data.matrix(dataDF)
  nRecord = nrow(data)
  x = array(0, dim = c(nRecord, xCols)) 
  y = data[, 1]
  for (k in 1 : xCols) {
    x[, k] = data[, k + 1]
  }
  uDF = read.csv(file="mockDataPrior.txt",head=FALSE,sep=",")
  uPrior = data.matrix(uDF)
} else if (loadData == "MyData") {
  dataDF = read.csv(file="Moz_T.csv",head=TRUE,sep=",")
  yCols = 3
  yLevel = 2 ^ yCols
  xCols = 6
  nRecord = nrow(dataDF)
  trainPercent = 0.8
  nRecordTrain = as.integer(nRecord*trainPercent)
  nRecordTest = nRecord - nRecordTrain
  allIdx = 1:nRecord
  trainIdx = sample(allIdx, nRecordTrain)
  testIdx = trainIdx * (-1)
  x = array(1, dim=c(nRecordTrain, xCols)) 
  y = array(0, dim=c(nRecordTrain, 1)) 
  y1 = dataDF$Blood.test[trainIdx]
  y2 = dataDF$Heard[trainIdx]
  y3 = dataDF$H_Comp[trainIdx]
  x[,2] = dataDF$Wealth[trainIdx]
  x[,3] = dataDF$Religion[trainIdx]
  x[,4] = dataDF$Marital.St[trainIdx]
  x[,5] = dataDF$Sup.Soc[trainIdx]
  x[,6] = dataDF$Thin.Risk[trainIdx]
  y = y1*4 + y2*2 + y3 + 1
  
  xTest = array(1, dim=c(nRecordTest, xCols)) 
  yTest = array(0, dim=c(nRecordTest, 1)) 
  y1Test = dataDF$Blood.test[testIdx]
  y2Test = dataDF$Heard[testIdx]
  y3Test = dataDF$H_Comp[testIdx]
  xTest[,2] = dataDF$Wealth[testIdx]
  xTest[,3] = dataDF$Religion[testIdx]
  xTest[,4] = dataDF$Marital.St[testIdx]
  xTest[,5] = dataDF$Sup.Soc[testIdx]
  xTest[,6] = dataDF$Thin.Risk[testIdx]
  yTest = y1Test*4 + y2Test*2 + y3Test + 1
} else { #generate new mock data
  xCols = 4
  yCols = 2
  yLevel = 2 ^ yCols
  betaTrue = array(0, dim=c(yLevel, xCols))
  uPrior = array(0, dim=c(yLevel, xCols)) 
  uPrior[1,1] = 1.0
  uPrior[1,2] = 2.0
  uPrior[1,3] = 3.0
  uPrior[1,4] = 4.0
  uPrior[2,1] = 1.2
  uPrior[2,2] = 2.2
  uPrior[2,3] = 3.2
  uPrior[2,4] = 4.2
  uPrior[3,1] = 1.4
  uPrior[3,2] = 2.4
  uPrior[3,3] = 3.4
  uPrior[3,4] = 4.4
  uPrior[4,1] = 0
  uPrior[4,2] = 0
  uPrior[4,3] = 0
  uPrior[4,4] = 0
  for (j in 1 : yLevel) {
    for (k in 1 : xCols) {
      betaTrue[j, k] = rnorm(1, uPrior[j,k], sigmaPrior)
    }
  }
  nRecord = 1000
  ng = nRecord / yLevel
  x = array(1, dim=c(nRecord, xCols)) 
  for (k in 2 : xCols) {
    x[,k] = runif(nRecord, -1, 1)
  }
  y = rep(0, nRecord) 
  z = rep(0, yLevel)
  p = rep(0, yLevel)
  for (h in 1 : nRecord) {
    ezsum = 0
    for(j in 1 : yLevel)
    {
      z[j] = exp(sum(betaTrue[j,] * x[h,]))
      ezsum = ezsum + z[j]
    }
    for(j in 1 : yLevel)
    {
      p[j] = z[j]/ezsum
    }
    
    uu = runif(1, min = 0, max = 1)
    pacum = 0
    for(j in 1 : yLevel)
    {
      pacum = p[j]+pacum
      if(uu<=pacum)
      {
        y[h] = j
        break
      }
    }

  } 
  data = cbind(y, x)
  write.table(data, file = "mockData.txt", sep = ",", row.names = FALSE, col.names = FALSE)
  #write.table(uPrior, file = "mockDataPrior.txt", sep = ",", row.names = FALSE, col.names = FALSE)
  #write.table(betaTrue, file = "mockDataPriorTrue.txt", sep = ",", row.names = FALSE, col.names = FALSE)
}

#======   sampling beta   ======
nSamples = 5000 # the length of the MCMC samples
sigmaProposal = 0.05;
betaCurr = array(0, dim=c(yLevel, xCols)) 
betaNew = array(0, dim=c(yLevel, xCols)) 
betaPosCurr = array(0, dim=c(yLevel, xCols)) 
betaPosNew = array(0, dim=c(yLevel, xCols))

betaSamplesLikelihood = array(0, dim=c(nSamples, yLevel, xCols)) 
betaSamplesPosterior = array(0, dim=c(nSamples, yLevel, xCols)) 

for (i in 1 : nSamples) {
  if (i %% 50 ==0) {
    print(i)
  }
  for (k in 1 : xCols) {
    for (j in 1 : yLevel-1) {
      betaNew[j, k] = rnorm(1, betaCurr[j, k], sigmaProposal)
    }
    betaNew[yLevel, k] = betaCurr[yLevel, k]
    logLikelihoodNew = logLikelihood(betaNew, x, y)
    logLikelihoodCurr = logLikelihood(betaCurr, x, y)
    logAlpha = logLikelihoodNew - logLikelihoodCurr
    alpha = exp(logAlpha)
    #print(alpha)
    u = runif(1, min = 0, max = 1)
    if (u < min(1.0, alpha)) {
      betaCurr[,k] = betaNew[,k]
    }
    #sample posterior
    for (j in 1 : yLevel-1) {
      betaPosNew[j, k] = rnorm(1, betaPosCurr[j, k], sigmaProposal)
    }
    betaPosNew[yLevel, k] = betaPosCurr[yLevel, k]
    if (loadData == "MyData") {
      priorMean = array(0, dim=c(yCols, xCols)) 
    } else {
      priorMean = uPrior
    }
    logPosteriorNew = logLikelihood(betaPosNew, x, y) + logPrior(betaPosNew, priorMean)
    logPosteriorCurr = logLikelihood(betaPosCurr, x, y) + logPrior(betaPosCurr, priorMean)
    logAlphaPosterior = logPosteriorNew - logPosteriorCurr
    alphaPosterior = exp(logAlphaPosterior)
    u = runif(1, min = 0, max = 1);
    if (u < min(1.0, alphaPosterior)) {
      betaPosCurr[, k] = betaPosNew[, k];
    }
  }
  
  betaSamplesLikelihood[i,,] = betaCurr
  betaSamplesPosterior[i,,] = betaPosCurr
}

#======   display beta samples   ======
saveRDS(betaSamplesLikelihood, file="betaSamplesLikelihood.Rds")
saveRDS(betaSamplesPosterior, file="betaSamplesPosterior.Rds")
#betaSamplesLikelihood = readRDS("betaSamplesLikelihood.Rds")
#betaSamplesPosterior = readRDS("betaSamplesPosterior.Rds")

burnin = 500
for (j in 1 : yLevel) {
  for (k in 1 : xCols) {
    #display likelihood estimation
    # trace
    beta = betaSamplesLikelihood[, j, k]
    traceTitle = sprintf("Trace of beta%i%i for likelihood", j, k)
    fName = paste("output/", traceTitle, ".png", sep="")
    png(filename = fName)
    plot(beta, type = "l")
    title(traceTitle)
    dev.off()
    # histogram
    beta = beta[burnin : nSamples]
    histTitle = sprintf("Histogram of beta%i%i for likelihood", j, k)
    fName = paste("output/", histTitle, ".png", sep="")
    png(filename = fName)
    hi = hist(beta, breaks = 15, freq = FALSE, main = NULL)
    lines(density(beta))
    title(histTitle)
    dev.off()
    idx = which.is.max(hi$density)
    bMLK = hi$mids[idx]
    estimate = sprintf("MLK estimation of beta%i%i is %f;", j, k, bMLK)
    print(estimate)
    #display posterior estimation
    # trace
    beta = betaSamplesPosterior[, j, k]
    traceTitle = sprintf("Trace of beta%i%i for posterior", j, k)
    fName = paste("output/", traceTitle, ".png", sep="")
    png(filename = fName)
    plot(beta, type = "l")
    title(traceTitle)
    dev.off()
    # histogram
    beta = beta[burnin : nSamples]
    histTitle = sprintf("Histogram of beta%i%i for posterior", j, k)
    fName = paste("output/", histTitle, ".png", sep="")
    png(filename = fName)
    hi = hist(beta, breaks = 15, freq = FALSE, main = NULL)
    lines(density(beta))
    title(histTitle)
    dev.off()
    idx = which.is.max(hi$density)
    bMAP = hi$mids[idx]
    estimate = sprintf("MAP estimation of beta%i%i is %f;", j, k, bMAP)
    print(estimate)
  }
}

#overlay prior, likelihood and posterior
m = nSamples - burnin + 1
n0 = rnorm(m, 0, sigmaPrior)
if (loadData != "MyData") {
  for (j in 1 : yLevel) {
    for (k in 1 : xCols) {
      compareTitle = sprintf("Compare of beta%i%i for prior, likelihood and posterior", j, k)
      fName = paste("output/", compareTitle, ".png", sep="")
      png(filename = fName)
      betaLk = betaSamplesLikelihood[, j, k]
      betaLk = betaLk[burnin : nSamples]

      betaPost = betaSamplesPosterior[, j, k]
      betaPost = betaPost[burnin : nSamples]
      
      betaPrior = n0 + uPrior[j, k]# + 0.1
      cond = factor( rep(c("Prior", "Likelihood", "Posterior"), each = m) )
      
      dataDensity = data.frame(cond, betac = c(betaPrior, betaLk, betaPost))
      #sm.density.compare(dataDensity$betac, dataDensity$cond)
      dpr = density(betaPrior)
      dlk = density(betaLk)
      dpo = density(betaPost)
      yAll = c(dpr$y, dlk$y, dpo$y)
      yMax = max(yAll)
      plot(dpr$x,dpr$y,type = "l",col="red", xlab = "beta", ylab="density", ylim=c(0, yMax))
      #par(new=TRUE)
      lines(dlk$x,dlk$y,col="green")
      #par(new=TRUE)
      lines(dpo$x,dpo$y, col="blue",add=TRUE)
      legend("topright", c("Prior", "Likelihood", "Posterior"), fill=2+(0:nlevels(dataDensity$cond)))
      title(compareTitle)
      dev.off()

    }
  }
}
