
fName = "output1/a.png"
png(filename = fName)

dpr = density(betaPrior)
dlk = density(betaLk)
dpo = density(betaPost)
x1=dpr$x
x2= dlk$x
x3= dpo$x
y3=dpo$y
y1=dpr$y
y2= dlk$y
plot(x1,y1,ylim=range(c(y1,y2,y3)),xlim=range(c(x1,x2,x3)), type="l",col="red")
lines(x2,y2,col="green")
lines(x3,y3,col="blue")

legend("topright", levels(dataDensity$cond), fill=2+(0:nlevels(dataDensity$cond)))

dev.off()

fName = "output1/apo.png"
png(filename = fName)

hist(betaPost, breaks = 15, freq = FALSE, main = NULL)
dpo = density(betaPost)

lines(dpo,col="blue")

legend("topright", levels(dataDensity$cond), fill=2+(0:nlevels(dataDensity$cond)))

dev.off()
