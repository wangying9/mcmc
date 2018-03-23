# MCMC fit for a multi-class logistic regression model
# author: Ying Wang
# version: 1.0

ls()
rm(list=ls())
gc()

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
 
    for (h in 1 : nrow(y)) {
      z1 = sum(beta[1,] * x[h,])
      z2 = sum(beta[2,] * x[h,])
      ezsum = exp(z2) + exp(z1)
      p1 = exp(z1) / ezsum
      p2 = exp(z2) / ezsum
      if(y==1){
        likelihoodSample = log(p1)
      } else {
        likelihoodSample = log(p2)
      }
      bb = beta[1,] - beta[2,] 
      zz = sum(bb * x[h,])
      pp = exp(zz) / (1.0 + exp(zz))
      if(y[h]==1){
        yy = 1
      }else {
        yy = 0
      }
      likelihoodSample1 = log(pp) * yy + log(1.0 - pp) * (1.0 - yy)
      sd = likelihoodSample -likelihoodSample1
      print(sd)
      message(paste('ok made it this far with x=',sd))
      show(sd)
      likelihood = likelihood + likelihoodSample1
    }
 
  return(likelihood)
}

#======   prepare data   ======
loadData = "MockData1" # or "MyData" 
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
}  else { #generate new mock data
  xCols = 4
  yCols = 2
  betaTrue = array(0, dim=c(yCols, xCols))
  uPrior = array(0, dim=c(yCols, xCols)) 
  uPrior[1,1] = 1.0
  uPrior[1,2] = 2.0
  uPrior[1,3] = 3.0
  uPrior[1,4] = 4.0
  uPrior[2,1] = 0.5
  uPrior[2,2] = 1.5
  uPrior[2,3] = 2.5
  uPrior[2,4] = 3.5
  for (j in 1 : yCols) {
    for (k in 1 : xCols) {
      betaTrue[j, k] = rnorm(1, uPrior[j,k], sigmaPrior)
    }
  }
  nRecord = 100
  x = array(1, dim=c(nRecord, xCols)) 
  y = array(0, dim=c(nRecord, 1)) 
  for (k in 2 : xCols) {
    x[,k] = runif(nRecord, -1, 1)
  }
  for (h in 1 : nRecord) {
    
      z1 = sum(betaTrue[1,] * x[h,])
      z2 = sum(betaTrue[2,] * x[h,])
      ezsum = exp(z1) + exp(z2) 
      p1 = exp(z1) / ezsum
# bb = betaTrue[1,] - betaTrue[2,] 
# zz = sum(bb * x[h,])
# p2 = exp(zz)/(1+exp(zz))
# p3 = (exp(z1)/exp(z2)) /(ezsum/exp(z2))
# p4 = (exp(z1-z2)) /(1+exp(z1-z2))
# pd = p1 -p2
# print(pd)
      uu = runif(1, min = 0, max = 1)
      if(uu < p1){
        y[h] = 1
      } else {
        y[h] = 2
      }
    
  }
  data = cbind(y, x)
  write.table(data, file = "mockData.txt", sep = ",", row.names = FALSE, col.names = FALSE)
  write.table(uPrior, file = "mockDataPrior.txt", sep = ",", row.names = FALSE, col.names = FALSE)
}

#======   sampling beta   ======
nSamples = 1000 # the length of the MCMC samples
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
      #sample posterior
      betaPosNew[j, k] = betaPosCurr[j, k] + sigmaProposalPost * rnorm(1, mean = 0, sd = 1)
      if (loadData == "MyData") {
        priorMean = array(0, dim=c(yCols, xCols)) 
      } else {
        priorMean = uPrior
      }
      logPosteriorNew = logLikelihood(betaPosNew, x, y) + logPrior(betaPosNew, priorMean)
      logPosteriorCurr = logLikelihood(betaPosCurr, x, y) + logPrior(betaPosCurr, priorMean)
      logAlphaPosterior = logPosteriorNew - logPosteriorCurr
      alphaPosterior = exp(logAlphaPosterior)
      #print(alphaPosterior)
      u = runif(1, min = 0, max = 1);
      if (u < min(1.0, alphaPosterior)) {
        betaPosCurr[j, k] = betaPosNew[j, k];
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
m = nSamples - burnin + 1
n0 = rnorm(m, 0, sigmaPrior)
if (loadData != "MyData") {
  for (j in 1 : yCols) {
    for (k in 1 : xCols) {
      compareTitle = sprintf("Compare of beta%i%i for prior, likelihood and posterior", j, k)
      fName = paste("output/", compareTitle, ".png", sep="")
      png(filename = fName)
      betaLk = betaSamplesLikelihood[, j, k]
      betaLk = betaLk[burnin : nSamples]

      betaPost = betaSamplesPosterior[, j, k]
      betaPost = betaPost[burnin : nSamples]
      betaPrior = n0 + uPrior[j, k]

      cond = factor( rep(c("Prior", "Likelihood", "Posterior"), each = m) )
      
      dataDensity = data.frame(cond, betac = c(betaPrior, betaLk, betaPost))
      
      dpr = density(betaPrior)
      dlk = density(betaLk)
      dpo = density(betaPost)
      yAll = c(dpr$y, dlk$y, dpo$y)
      yMax = max(yAll)
      plot(dpr$x,dpr$y,type = "l",col="red", xlab = "beta", ylab="density", ylim=c(0, yMax))
      lines(dlk$x,dlk$y,col="green")
      lines(dpo$x,dpo$y, col="blue")
      legend("topright", c("Prior", "Likelihood", "Posterior"), fill=2+(0:nlevels(dataDensity$cond)))
      title(compareTitle)
      dev.off()

    }
  }
}
