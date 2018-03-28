# MCMC fit for a multi-class logistic regression model
# author: Ying Wang
# version: 1.0

ls()
rm(list=ls())
gc()
library(MCMCpack)
library(nnet)
library(sm)

nSample = 10000
trainingDF = read.csv(file="trainingData.csv",head=TRUE,sep=",")
testingDF = read.csv(file="testingData.csv",head=TRUE,sep=",")


post2<- MCMCmnl(Y ~
                  X1 + X2 + X3 + X4 + X5,
                baseline=1, mcmc.method="IndMH", B0=0,
                verbose=500, mcmc=nSample, thin=1, tune=0.5,
                data=trainingDF)
#plot(post2)
#summary(post2)
betaSamples = data.matrix(post2)
burnin = 500
yLevel = 7
xCols = 6 
betaEstLk = array(0, dim=c(yLevel, xCols)) 
for (k in 1 : xCols) {
  for (j in 1 : yLevel) {
  
    i = j + (k-1)*yLevel
    beta = betaSamples[, i]
    traceTitle = sprintf("Trace of beta%i%i for likelihood", j, k)
    fName = paste("output/", traceTitle, ".png", sep="")
    png(filename = fName)
    plot(beta, type = "l")
    title(traceTitle)
    dev.off()
    # histogram
    beta = beta[burnin : nSample]
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
    betaEstLk[j,k] = as.double(bMLK)

  }
}


testingData = data.matrix(testingDF)
y= testingData[, 1]
  N = nrow(testingData)
  one = array(1, dim=c(N, 1)) 
  x = cbind(one, testingData[,2:6])
  L = yLevel
  cnt = 0
  for (h in 1 : N) {
    sumez = 0.0
    z = rep(0, L)
    p = rep(0, L+1)
    for (j in 1 : L) {
      z[j] = exp(sum(betaEstLk[j,] * x[h,]))
      sumez = sumez + z[j]
    }
    sump = 0.0
    for (j in 1 : L) {
      p[j+1] = z[j]/(1+sumez)
      sump = sump + p[j]
    }
    p[1] = 1 - sump
    yPred = which.is.max(p)
    if (yPred == y[h]){
      cnt = cnt + 1
    }
  }
  acc = cnt / N
  accuracy = sprintf("MLK accuracy is %f;", acc)
  print(accuracy)
  