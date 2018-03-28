# MCMC fit for a multi-class logistic regression model
# author: Ying Wang
# version: 1.0

ls()
rm(list=ls())
gc()


dataDF = read.csv(file="Moz_T.csv",head=TRUE,sep=",")
yCols = 3
yLevel = 2 ^ yCols
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
y = y1*4 + y2*2 + y3 + 1

xTest = array(1, dim=c(nRecordTest, xCols)) 
yTest = array(0, dim=c(nRecordTest, 1)) 
y1Test = dataDF$Blood.test[testIdx]
y2Test = dataDF$Heard[testIdx]
y3Test = dataDF$H_Comp[testIdx]
xTest[,2] = dataDF$Wealth[testIdx]
xTest[,3] = religion[testIdx]
xTest[,4] = dataDF$Marital.St[testIdx]
xTest[,5] = dataDF$Sup.Soc[testIdx]
xTest[,1] = dataDF$Thin.Risk[testIdx]
yTest = y1Test*4 + y2Test*2 + y3Test + 1

colNames = vector()
colNames = c(colNames, "Y")

for (k in 1 : xCols) {
  colName = paste('X', k, sep = "")
  colNames = c(colNames, colName)
}
dataTraining = cbind(y, x)
write.table(dataTraining, file = "trainingData.csv", sep = ",", row.names = TRUE, col.names = colNames)
dataTesting = cbind(yTest, xTest)
write.table(dataTesting, file = "testingData.csv", sep = ",", row.names = TRUE, col.names = colNames)
