data(Nethvote)
## just a choice-specific X var

post2<- MCMCmnl(vote ~
                  relig + class + income + educ + age + urban,
                baseline="D66", mcmc.method="IndMH", B0=0,
                verbose=500, mcmc=10000, thin=10, tune=0.5,
                data=Nethvote)
plot(post2)
summary(post2)


X     = model.matrix(Species ~ ., data=iris);
y.all = model.matrix(~ Species - 1, data=iris);
J = nlevels(iris$Species)
y     = y.all[,-J];
