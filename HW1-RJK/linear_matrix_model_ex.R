A <- matrix(c(0, 3.75, 0.5, 0.5), nrow=2, ncol=2, byrow=T)
myeigs <- eigen(A)
myeigs$vectors[,1]
myeigs

initcond <- c(1,1)
A %*% initcond

maxtime <- 50

out <- matrix(rep(NA,2*(maxtime+1)), nrow = (maxtime+1))
out[1,] <- initcond
for(i in 1:maxtime){
  out[i+1,] <- A %*% out[i,]
}

out <- as.data.frame(out)
names(out) <- c("juv", "adult")
out$time <- 0:maxtime
out$N <- out$juv + out$adult

## plot total population
plot(N~time, out, type="b", col="blue")
points(juv~time, out, type="b", col="red")
points(adult~time, out, type="b", col="purple")

legend("topleft", c("N", "juv", "adults"), col=c("blue", "red", "black"), lty=1)

## plot relative population growth rate N(t+1)/N(t)
out$growthrate <- c(NA,out$N[2:(maxtime+1)]/out$N[1:(maxtime)])
out$growthrate
plot(growthrate~time, out)

## plot proportion of juveniles
plot((juv/N)~time, out)

propjuvtheory <- myeigs$vectors[,1][1]/sum(myeigs$vectors[,1])
abline(h=propjuvtheory, col="red", lty=2)