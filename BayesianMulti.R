

# MCMC fit for a logistic regression model
ls()
rm(list=ls())
gc()

logprior <- function(beta) {
  p = 1
  for(j in 1 : length(beta)) {
    p = p * dnorm(beta[j],0, 0.00001)
  }
  
  #p = dmvnorm(beta, u, sigma)
  return(p)
}

#load data
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


cols=4;# one intercept + 3 attributes
rows=10000;
sigma2=0.05;

datax <- read.csv(file="R-simulated-logist-data-x.csv",head=FALSE,sep=",")
x1=datax$V1
x2=datax$V2
x3=datax$V3
beta_curr=seq(0.5,cols)/1.0
beta_new=seq(0.5,cols)/1.0
datay <- read.csv(file="R-simulated-logist-data-y.csv",head=FALSE,sep=",")
y=datay$V1
d1 <- data.frame(1.5, t(beta_curr))                
write.table(d1, "R-mcmc_logistic_version_2.txt", row.names = FALSE)
m=5000 # the length of the MCMC samples
betaSamples =array(0, dim=c(m,cols)) 
for (i in 1:m)
{
  if (i%%500 ==0){
    print(i)
  }
  for (k in 1:cols)
  {
    temp1=beta_curr[k] +sigma2*rnorm(1, mean = 0, sd = 1)
    beta_new[k]=temp1;
    likelihood=0.0;
    for (h in 1:rows)
    {
      temp_new=beta_new[1] +beta_new[2]*x1[h] +beta_new[3]*x2[h]+beta_new[4]*x3[h];
      temp_new=exp(temp_new);
      temp_curr=beta_curr[1] +beta_curr[2]*x1[h] +beta_curr[3]*x2[h]+beta_curr[4]*x3[h];
      temp_curr=exp(temp_curr);
      likelihood=likelihood+ log( temp_new /(1.0 +temp_new))*y[h];
      likelihood=likelihood+ log(1.0- temp_new /(1.0 +temp_new))*(1.0-y[h]);
      likelihood=likelihood- log(temp_curr /(1.0 +temp_curr))*y[h];
      likelihood=likelihood- log(1.0- temp_curr /(1.0 +temp_curr))*(1.0-y[h]);
    }
    
    likelihood=likelihood + logprior(temp_new)
    likelihood=likelihood - logprior(temp_curr)
    u=runif(1, min = 0, max = 1);
    if ( u<min(1.0, exp(likelihood)))
    {
      beta_curr[k]=beta_new[k];
    }
  }
  betaSamples[i,] = beta_curr
  d1 <- data.frame(likelihood, t(beta_curr))
  write.table(d1, "R-mcmc_logistic_version_2.txt", row.names = FALSE, col.names = FALSE, append = TRUE)
}
b0a = betaSamples[,1]
plot(b0a, type = "l")
burnin = 500
b0=b0a[burnin:m]
h<-hist(b0,breaks=15, freq=FALSE)
lines(density(b0))

idx = which.is.max(h$density)
b0map = h$mids[idx]





 
 