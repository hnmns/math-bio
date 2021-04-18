#PROBLEM 3

data <- c(0, 0.0043, 0.1132, 0, 0.9775, 0.9111, 0, 0, 0, 0.0736, 0.9534, 0, 0, 0, 0.0452, 0.9804)
proj_matr <- matrix(data, nrow = 4, ncol = 4, byrow = TRUE)
print(proj_matr)
eigs <- eigen(proj_matr)
#eigen returns named list with eigenvalues and eigenvectors, labeled vectors and values
eigs
values <- eigs$values
vectors <- eigs$vectors
#a) Return dominant eigenvalue and stable-stage distribution
d_eigval <- max(values) #dominant eigval is [1] by default anyway
print(d_eigval) #dominant eigenvalue
d_evec <- eigs$vectors[,1]
w <- d_evec/(sum(d_evec))
#make components sum to 1 for stable-stage distribution, w
print(w)

init_cond <- c(10,60,110,70)

d_eigval <- max(values) 
d_eigval #dominant eigenvalue

maxtime <- 50

out <- matrix(rep(NA,4*(maxtime+1)), nrow = (maxtime+1))
out[1,] <- init_cond
for(i in 1:maxtime){
	out[i+1,] <- proj_matr %*% out[i,]
}

#3.b) projecting population dynamics, storing in 'out'
out <- as.data.frame(out)
names(out) <- c("year","juv","mature","postrep")
out$time <- 0:maxtime
out$N <- out$year + out$juv + out$mature + out$postrep

#3.c) plotting
#graph total population, N
plot(N~time, out, type="b", xlab = "Time (years)", ylab = "Total Population", 
		 main = "Total Population", col="blue")
#~ creates dep-indep relation

#annual population growth rate, N(t+1)/N(t)
annual_growth_rates <- c(NA, out$N[2:(maxtime+1)]/out$N[1:maxtime])
plot(annual_growth_rates~time, out, type="b", xlab = "Time (years)",
		 ylab = "Annual Population Growth Rate", main = "Annual Population Growth Rate",
		 col="red")

#proportion of individuals at each stage
plot(year~time, out, type="b", main = "Population Proportions", ylab = "Population", 
		 xlab = "Time (years)", ylim = c(0,350), col="green") 
#'points' puts given data on previous plot
points(juv~time, out, type="b", col="black") 
points(mature~time, out, type="b", col="blue")
points(postrep~time, out, type="b", col="red")
legend("topleft", c("yearling", "juvenile","mature","postreproductive"), 
			 col=c("green", "black", "blue", "red"), lty=1)

#PROBLEM 4
init_cond_w <- 250*w #initial condition with stable-stage dist and 250 total
print(init_cond_w)

#Juvenile harvests
#make matrix 'harvests' to store 10 population projects with 1:10 juv harvests
juv_harvests <- matrix(rep(NA, 10*(maxtime+1)), nrow = (maxtime+1))
harvest_max <- 10
for(j in 1:harvest_max){
	out2 <- matrix(rep(NA,4*(maxtime+1)), nrow = (maxtime+1))
	out2[1,] <- init_cond_w
	for(i in 1:maxtime){
		out2[i+1,] <- ((proj_matr %*% out2[i,]) - c(0,j,0,0))
		#subtracting constant juv harvest j from each generation after breeding season
	}
	out2 <- as.data.frame(out2)
	names(out2) <- c("year","juv","mature","postrep")
	out2$N <- out2$year + out2$juv + out2$mature + out2$postrep
	juv_harvests[,j] <- out2$N
}
juv_harvests <- as.data.frame(juv_harvests)
names(juv_harvests) <- c("h1","h2","h3","h4","h5","h6","h7","h8","h9","h10")
juv_harvests$time <- 0:maxtime
plot(h1~time, juv_harvests, ylab = "Population", xlab="Time (Years)", 
		 main="Total Population Projections for Various Junior Harvests",type="l", col=10)
points(h2~time, juv_harvests, type="p", col=1)
points(h3~time, juv_harvests, type="l", col=2)
points(h4~time, juv_harvests, type="p", col=3)
points(h5~time, juv_harvests, type="l", col=4)
points(h6~time, juv_harvests, type="p", col=5)
points(h7~time, juv_harvests, type="l", col=6)
points(h8~time, juv_harvests, type="p", col=7)
points(h9~time, juv_harvests, type="l", col=8)
points(h10~time, juv_harvests, type="p", col=9)

legend("topleft", c("h=1","h=2","h=3","h=4","h=5","h=6","h=7","h=8","h=9","h=10"),
			 col=c(10,1:9), lty=rep(c(1,2),times=5))

#Now, for harvest of mature/reproductive adults
mature_harvests <- matrix(rep(NA, 10*(maxtime+1)), nrow = (maxtime+1))
harvest_max <- 10
for(j in 1:harvest_max){
	out2 <- matrix(rep(NA,4*(maxtime+1)), nrow = (maxtime+1))
	out2[1,] <- init_cond_w
	for(i in 1:maxtime){
		out2[i+1,] <- ((proj_matr %*% out2[i,]) - c(0,0,j,0))
	}
	out2 <- as.data.frame(out2)
	names(out2) <- c("year","juv","mature","postrep")
	out2$N <- out2$year + out2$juv + out2$mature + out2$postrep
	mature_harvests[,j] <- out2$N
}
mature_harvests <- as.data.frame(mature_harvests)
names(mature_harvests) <- c("h1","h2","h3","h4","h5","h6","h7","h8","h9","h10")
mature_harvests$time <- 0:maxtime
plot(h1~time, mature_harvests, ylab = "Population", xlab="Time (Years)", 
		 main="Total Population Projections for Various Mature Harvests",type="l", col=10)
points(h2~time, mature_harvests, type="p", col=1)
points(h3~time, mature_harvests, type="l", col=2)
points(h4~time, mature_harvests, type="p", col=3)
points(h5~time, mature_harvests, type="l", col=4)
points(h6~time, mature_harvests, type="p", col=5)
points(h7~time, mature_harvests, type="l", col=6)
points(h8~time, mature_harvests, type="p", col=7)
points(h9~time, mature_harvests, type="l", col=8)
points(h10~time, mature_harvests, type="p", col=9)

legend("topleft", c("h=1","h=2","h=3","h=4","h=5","h=6","h=7","h=8","h=9","h=10"),
			 col=c(10,1:9), lty=rep(c(1,2),times=5))