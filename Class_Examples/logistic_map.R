#### Numerically solve the logistic map equation: X(t+1) = r*X(t)*(1-X(t))
### Then make bifurcation diagram

## first start with a single value for r
r <- 3.9
maxtime <- 1000
times <- 1:maxtime
x <- rep(NA, maxtime+1)
x[1] <- 0.1
### now loop through times and save value of x
for(i in 1:maxtime){
  x[i+1] <- r*x[i]*(1-x[i])
}

## now plot it!
out <- data.frame(x=x, time=0:maxtime)
plot(x~time, out,  type="b")
abline(h=0)


###### Now let's make the bifurcation diagram!!
#### Remember what we need: for each r value, we need the steady state value of x
## We want to loop through values of r and save the last handful of values of x after it settles out

rs <- seq(0,4,0.01)

out <- NULL
for(r in rs){
  ### solve the model with this r
  for(i in 1:maxtime){
    x[i+1] <- r*x[i]*(1-x[i])
  }
  ### Now save the last 100 of values of x and the value r
  outies <- x
  outies <- data.frame(tail(x, n=100))
  outies$r <- r
  
  out <- rbind(out, outies)
  
}
names(out) <- c("x", "r")

plot(x~r, out, pch=19, cex=0.05)
