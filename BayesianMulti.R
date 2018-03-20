# MCMC fit for a multi-class logistic regression model
# author: Ying Wang
# version: 1.0

ls()
rm(list=ls())
gc()

#======   define functions   ======
sigmaPrior = 0.01
logPrior = function(beta, uPrior) {
  p = 1
  for(j in 1 : nrow(uPrior)) {
    for(k in 1 : ncol(uPrior)) {
      p = p * dnorm(beta[j, k], uPrior[j,k], sigmaPrior)
    }
  }
  return(p)
}

logLikelihood = function(beta, x, y) {
  likelihood = 0.0
  for (j in 1 : ncol(y)) {
    for (h in 1 : nrow(y)) {
      z = sum(beta[j,] * x[h,])
      p = exp(z) / (1.0 + exp(z))
      likelihoodSample = log(p) * y[h, j] + log(1.0 - p) * (1.0 - y[h, j])
      likelihood = likelihood + likelihoodSample
    }
  }
  return(likelihood)
}

#======   prepare data   ======
loadData = "MockData" # or "MyData" 
if (loadData == "MockData") {
  dataDF = read.csv(file="mockData.txt",head=TRUE,sep=",")
  colNames = colnames(dataDF)
  idx = 0
  for (i in 1 : length(colNames)) {
    char = substring(colNames[i], 1, 1)
    if (char=="X") {
      idx = i
      break
    }
  }
  yCols = idx - 1
  xCols = length(colNames) - yCols
  data = data.matrix(dataDF)
  nRecord = nrow(data)
  x = array(0, dim=c(nRecord, xCols)) 
  y = array(0, dim=c(nRecord, yCols)) 
  for (j in 1 : yCols) {
    y[,j] = data[, j]
  }
  for (k in 1 : xCols) {
    x[,k] = data[, k+yCols]
  }
  uDF = read.csv(file="mockDataPrior.txt",head=TRUE,sep=",")
  uPrior = data.matrix(uDF)
} else if (loadData == "MyData") {
  dataDF = read.csv(file="Moz_T.csv",head=TRUE,sep=",")
  yCols = 3
  xCols = 6
  nRecord = nrow(dataDF)
  x = array(1, dim=c(nRecord, xCols)) 
  y = array(0, dim=c(nRecord, yCols)) 
  y[,1] = dataDF$Blood.test
  y[,2] = dataDF$Heard
  y[,3] = dataDF$H_Comp
  x[,2] = dataDF$Wealth
  x[,3] = dataDF$Religion
  x[,4] = dataDF$Marital.St
  x[,5] = dataDF$Sup.Soc
  x[,6] = dataDF$Thin.Risk
} else { #generate new mock data
  xCols = 4
  yCols = 2
  betaTrue = array(0, dim=c(yCols, xCols))
  uPrior = array(0, dim=c(yCols, xCols)) 
  uPrior[1,1] = 1
  uPrior[1,2] = 2
  uPrior[1,3] = 3
  uPrior[1,4] = 4
  uPrior[2,1] = 6
  uPrior[2,2] = 7
  uPrior[2,3] = 5
  uPrior[2,4] = 3
  for (j in 1 : yCols) {
    for (k in 1 : xCols) {
      betaTrue[j, k] = rnorm(1, uPrior[j,k], sigmaPrior)
    }
  }
  nRecord = 1000
  x = array(1, dim=c(nRecord, xCols)) 
  y = array(0, dim=c(nRecord, yCols)) 
  for (k in 2 : xCols) {
    x[,k] = runif(nRecord, -1, 1)
  }
  for (h in 1 : nRecord) {
    for (j in 1 : yCols) {
      z = sum(betaTrue[j,] * x[h,])
      p = exp(z) / (1 + exp(z))
      y[h,j] = rbern(1,p)
    }
  }
  colNames = vector()
  for (j in 1 : yCols) {
    colName = paste('Y', j, sep = "")
    colNames = c(colNames, colName)
  }
  for (k in 1 : xCols) {
    colName = paste('X', k, sep = "")
    colNames = c(colNames, colName)
  }
  data = cbind(y, x)
  write.table(data, file = "mockData.txt", sep = ",", row.names = FALSE, col.names = colNames)
  write.table(u, file = "mockDataPrior.txt", sep = ",", row.names = FALSE, col.names = FALSE)
}

#======   sampling beta   ======
nSamples = 5000 # the length of the MCMC samples
sigmaProposal = 0.05;
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

  for (j in 1 : yCols) {
    for (k in 1 : xCols) {
      #sample liklihood
      betaNew[j, k] = betaCurr[j, k] + sigmaProposal * rnorm(1, mean = 0, sd = 1)
      logLikelihoodNew = logLikelihood(betaNew, x, y)
      logLikelihoodCurr = logLikelihood(betaCurr, x, y)
      logAlpha = logLikelihoodNew - logLikelihoodCurr
      alpha = exp(logAlpha)
      u = runif(1, min = 0, max = 1);
      if ( u < min(1.0, alpha)) {
        betaCurr[j, k]=betaNew[j, k];
      }
      #sample posterior
      betaPosNew[j, k] = betaPosCurr[j, k] + sigmaProposal * rnorm(1, mean = 0, sd = 1)
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
      if ( u < min(1.0, alphaPosterior)) {
        betaPosCurr[j, k]=betaPosNew[j, k];
      }
    }
  }
  betaSamplesLikelihood[i,,] = betaCurr
  betaSamplesPosterior[i,,] = betaPosCurr
}

#======   display beta samples   ======
#saveRDS(betaSamplesLikelihood, file="betaSamplesLikelihood.Rds")
#saveRDS(betaSamplesPosterior, file="betaSamplesPosterior.Rds")
#betaSamplesLikelihood = readRDS("betaSamplesLikelihood.Rds")
#betaSamplesPosterior = readRDS("betaSamplesPosterior.Rds")

burnin = 500
for (j in 1 : yCols) {
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
    h<-hist(beta, breaks = 15, freq = FALSE, main = NULL)
    lines(density(beta))
    title(histTitle)
    dev.off()
    idx = which.is.max(h$density)
    bMLK = h$mids[idx]
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
    h = hist(beta, breaks = 15, freq = FALSE, main = NULL)
    lines(density(beta))
    title(histTitle)
    dev.off()
    idx = which.is.max(h$density)
    bMAP = h$mids[idx]
    estimate = sprintf("MAP estimation of beta%i%i is %f;", j, k, bMAP)
    print(estimate)
  }
}

#overlay prior, likelihood and posterior
if (loadData != "MockData") {
  for (j in 1 : yCols) {
    for (k in 1 : xCols) {
      compareTitle = sprintf("Compare of beta%i%i for prior, likelihood and posterior", j, k)
      fName = paste("output/", compareTitle, ".png", sep="")
      png(filename = fName)
      betaLk = betaSamplesLikelihood[, j, k]
      betaLk = betaLk[burnin : nSamples]

      betaPost = betaSamplesPosterior[, j, k]
      betaPost = betaPost[burnin : nSamples]

      lines(density(betaLk), col="red")
      lines(density(betaPost), col="blue")
      title(compareTitle)
      dev.off()

    }
  }
}

betaLk = betaSamplesLikelihood[, 1, 1]
betaLk = betaLk[burnin : nSamples]

betaPost = betaSamplesPosterior[, 1, 1]
betaPost = betaPost[burnin : nSamples]
m = nSamples - burnin + 1
cond = factor( rep(c("A", "B"), each = m) )

dataDensity = data.frame(cond, beta = c(betaLk,betaPost))
sm.density.compare(dataDensity$beta, dataDensity$cond)
dev.off()



 
 