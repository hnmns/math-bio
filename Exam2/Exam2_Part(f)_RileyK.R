library(deSolve)

###Exam 2, Part f
###parameters
N = 10000
bet1 = 2.5e-05
bet2 = 2.75e-05
gam1 = 0.08
gam2 = 0.08
delt = 0.02
endemic_eq_1 = N - (gam1+delt)/bet1 #from part (c), shown in (e)

prams <- c(N=N,bet1=bet1,bet2=bet2,gam1=gam1,gam2=gam2,delt=delt)

ic <- c(I_1 = 11, I_2 = 11) #first init cond., to be changed later
days <- seq(1,365*5,1)

sii_1 <- function(t, ic, prams){
	with(as.list(c(ic, prams)),{
		dI_1 <- I_1*(bet1*(N-I_1-I_2)-gam1-delt)
		dI_2 <- I_2*(bet2*(N-I_1-I_2)-gam2-delt)

		list(c(dI_1,dI_2))
	})
}

out <- as.data.frame(ode(y=ic, times=days, func=sii_1,parms=prams)) #solving

plot(I_1~time, out, col="orange3", type="l", main = "Both Strains Initially Rare",
		 xlab="Days",ylab="Number of Infected", ylim=c(0,N))
lines(I_2~time, out, col="blue2")
lines((I_1+I_2)~time, out, col="black")
legend("topright", c("I_1","I_2","Total Infected"), 
			 col=c("orange3","blue2","black"), lty=1)

plot((I_1/N)~time, out, col="purple2", type="l", 
		 main = "Both Strains Initially Rare",
		 xlab="Days", ylab="Fraction of population infected", ylim=c(0,1))
lines((I_2/N)~time, out, col="green3")
legend("topright", c("I_1/N", "I_2/N"), col=c("purple2","green3"), lty=1)


###Now, plot again for second initial condition
ic <- c(I_1 = endemic_eq_1-30, I_2 = 11) #somewhere near endemic eq for strain 1
sii_2 <- function(t, ic, prams){
	with(as.list(c(ic, prams)),{
		dI_1 <- I_1*(bet1*(N-I_1-I_2)-gam1-delt)
		dI_2 <- I_2*(bet2*(N-I_1-I_2)-gam2-delt)
		
		list(c(dI_1,dI_2))
	})
}

out <- as.data.frame(ode(y=ic, times=days, func=sii_2, parms=prams))

plot(I_1~time, out, col="orange3", type="l", main = "Start Near Strain 1 Endemic Equilibrium",
		 xlab="Days",ylab="Number of Infected", ylim=c(0,N))
lines(I_2~time, out, col="blue2")
lines((I_1+I_2)~time, out, col="black")
legend("topright", c("I_1","I_2","Total Infected"), 
			 col=c("orange3","blue2","black"), lty=1)

plot((I_1/N)~time, out, col="purple2", type="l", 
		 main = "Start Near Strain 1 Endemic Equilibrium",
		 xlab="Days", ylab="Fraction of population infected", ylim=c(0,1))
lines((I_2/N)~time, out, col="green3")
legend("topright", c("I_1/N", "I_2/N"), col=c("purple2","green3"), lty=1)