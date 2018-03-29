# MCMC fit for a multi-class logistic regression model
# author: Ying Wang
# version: 1.0


library(MCMCpack)
library(nnet)
library(sm)

nSample = 1000
dataDF = read.csv(file="Moz_T.csv",head=TRUE,sep=",")

xCols = 5
nRecord = nrow(dataDF)
trainPercent = 0.1
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
religion = dataDF$Religion
religion[which(religion==96)] = 0
x[,2] = dataDF$Wealth[trainIdx]
x[,3] = religion[trainIdx]
x[,4] = dataDF$Marital.St[trainIdx]
x[,5] = dataDF$Sup.Soc[trainIdx]
x[,1] = dataDF$Thin.Risk[trainIdx]


samplesLk = getLikelihoodSample(nSamples, y1, x) 
burnin = 500
s = samplesLk[burnin:nSamples,,]
betaLk = estimateModel(s, 1, xCols)

accuracy = predict(betaLk, y, x)

  print(accuracy)
  