

# MCMC fit for a logistic regression model
ls()
rm(list=ls())
gc()

#======   define functions   ======
logprior = function(beta) {
  p = 1
  for(j in 1 : length(beta)) {
    p = p * dnorm(beta[j],0, 0.00001)
  }
  #p = dmvnorm(beta, u, sigma)
  return(p)
}

logLikelihood = function(beta, x, y) {
  H = nrow(y)
  J = ncol(y)
  likelihood = 0.0
  for (j in 1 : J) {
    for (h in 1 : H) {
      z = sum(beta[j,] * x[h,])
      p = exp(z) / (1.0 + exp(z))
      likelihoodSample = log(p) * y[h, j] + log(1.0 - p) * (1.0 - y[h, j])
      likelihood = likelihood + likelihoodSample
    }
  }
  return(likelihood)
}

#======   prepare data   ======
loadData = "MockData"
idx = 0
if (loadData == "MockData") {
  dataDF = read.csv(file="mockData.txt",head=TRUE,sep=",")
  colNames = colnames(dataDF)
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
} else if (loadData == "MyData") {
  
} else {
  xCols = 4
  yCols = 2
  betaTrue = array(0, dim=c(yCols, xCols)) 
  betaTrue[1,1] = 1
  betaTrue[1,2] = 2
  betaTrue[1,3] = 3
  betaTrue[1,4] = 4
  betaTrue[2,1] = 6
  betaTrue[2,2] = 7
  betaTrue[2,3] = 5
  betaTrue[2,4] = 3
  nRecord = 1000
  x = array(0, dim=c(nRecord, xCols)) 
  y = array(0, dim=c(nRecord, yCols)) 
  for (k in 2 : xCols) {
    x[,k] = runif(nRecord, -1, 1)
  }
  x[,1] = rep(1, nRecord)
  for (h in 1 : nRecord) {
    for (j in 1 : yCols) {
      z = sum(betaTrue[j,] * x[h,])
      #print(z)
      p = exp(z) / (1 + exp(z))
      #print(p)
      #cat(sprintf("z = %f and p= %f\n", z, p))
      y[h,j] = rbern(1,p)
    }
  }
  colNames = vector()
  for (j in 1 : yCols) {
    colName = paste('Y', j, sep = "")
    colNames = c(colNames, colName)
  }
  for (kk in 1 : xCols) {
    colName = paste('X', kk, sep = "")
    colNames = c(colNames, colName)
  }
  data = cbind(y, x)
  write.table(data, file = "mockData.txt", sep = ",", row.names = FALSE, col.names = colNames)
}

#======   sampling beta   ======
nSamples = 500 # the length of the MCMC samples
sigmaProposal = 0.05;
init = seq(0.5, xCols)/1.0
betaCurr = array(0, dim=c(yCols, xCols)) 
betaNew = array(0, dim=c(yCols, xCols)) 
betaPosteriorCurr = array(0, dim=c(yCols, xCols)) 
betaPosteriorNew = array(0, dim=c(yCols, xCols)) 
for (j in 1 : yCols) {
  betaCurr[j,] = init
  betaNew[j,] = init
  betaPosteriorCurr[j,] = init
  betaPosteriorNew[j,] = init
}
betaSamplesLikelihood = array(0, dim=c(nSamples, yCols, xCols)) 
betaSamplesPosterior = array(0, dim=c(nSamples, yCols, xCols)) 

for (i in 1 : nSamples)
{
  if (i %% 50 ==0){
    print(i)
  }

  for (j in 1 : yCols) {
    for (k in 1 : xCols)
    {
      #sample liklihood
      betaNew[j, k] = betaCurr[j, k] + sigmaProposal * rnorm(1, mean = 0, sd = 1)
      logLikelihoodNew = logLikelihood(betaNew, x, y)
      logLikelihoodCurr = logLikelihood(betaCurr, x, y)
      logAlpha = logLikelihoodNew - logLikelihoodCurr
      alpha = exp(logAlpha)
      u = runif(1, min = 0, max = 1);
      if ( u < min(1.0, alpha))
      {
        betaCurr[j, k]=betaNew[j, k];
      }
      #sample posterior
      betaPosteriorNew[j, k] = betaPosteriorCurr[j, k] + sigmaProposal * rnorm(1, mean = 0, sd = 1)
      logPosteriorNew = logLikelihood(betaNew, x, y) + logprior(betaNew)
      logPosteriorCurr = logLikelihood(betaCurr, x, y) + logprior(betaCurr)
      logAlphaPosterior = logPosteriorNew - logPosteriorCurr
      alphaPosterior = exp(logAlphaPosterior)
      u = runif(1, min = 0, max = 1);
      if ( u < min(1.0, alphaPosterior))
      {
        betaPosteriorCurr[j, k]=betaPosteriorNew[j, k];
      }
    }
  }
  betaSamplesLikelihood[i,,] = betaCurr
  betaSamplesPosterior[i,,] = betaCurr
}

#======   display beta samples   ======
b0a = betaSamplesLikelihood[, 1, 1]
plot(b0a, type = "l")
burnin = 100
b0=b0a[burnin : nSamples]
h<-hist(b0,breaks=15, freq=FALSE)
lines(density(b0))

idx = which.is.max(h$density)
b0map = h$mids[idx]





 
 