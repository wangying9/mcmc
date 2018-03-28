ls()
rm(list=ls())
gc()

data(Nethvote)
## just a choice-specific X var
data0 = cbind(Nethvote$vote, Nethvote$relig, Nethvote$class, Nethvote$income, Nethvote$urban, Nethvote$age)
df = data.frame(data0)

post2<- MCMCmnl(X1 ~
                   X2 + X3 + X4 + X5 + X6,
                baseline="1", mcmc.method="IndMH", B0=0,
                verbose=500, mcmc=1000, thin=10, tune=0.5,
                data=df)
plot(post2)
summary(post2)
