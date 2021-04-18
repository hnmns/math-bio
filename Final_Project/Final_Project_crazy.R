library(deSolve)
library(matlib)
N = 10000

mu = 0.02
beta_u = 0.1
beta_s = 0.1
eta = 0.8
k = 0.3 #esc. from U_l to U_h
##Make rehab rate a function of policing, fear/spite toward law
sig = 0.8 #uptake rate into rehab, sigma
theta_1 = 0.9
theta_2 = 0.4 ##from U_l,U_r into Q
###1/theta_i = average duration of use, fraction of 12 months
r_1 = 0.01
r_2 = 0.07
r_3 = 0.04 #relapse rates
rho = 0.2 #amelioration U_h to U_l
alpha_1 = 0.6
phi = 1.4 #phi>1
alpha_2 = phi*alpha_1 #p2s contact rates, U_l and U_h
nu = 0.3 #1/nu = avg duration of drugs in circulation / supply chain
beta_hats <- beta_s*alpha_1/nu #transmission pram, S2P2S transmission cycle
epsilon = 0.001 #timescale parameter
cond_pol = 1 #conditional policing, if many heavy users, bigger crackdown policy 

beta_hats = beta_s*alpha_1/nu #NOT in prams yet
Q_1 = mu + k + theta_1
Q_2 = mu + sig + rho
Q_3 = mu + r_1 + theta_2
Q_4 = mu + r_2 + r_3

#My originally determined, 100% accurate expression for R_0, don't worryabout it 
#R_0 prime, not the original one from Nyabadza et al.
#R_0p <- beta_u*((Q_2*Q_3-r_1*sig+eta*Q_3*k)/(Q_3*(Q_1*Q_2-k*rho)-Q_1*r_1*sig))
	

time <- seq(0, 250, 1)
state <- c(s=0.95,v=0.04,w=0.005,x=0.004,y=0.001,z=0.05)
prams <- c(mu=mu, beta_u=beta_u, beta_s=beta_s, eta=eta, k=k, sig=sig, 
					 theta_1=theta_1, theta_2=theta_2, r_1=r_1, r_2=r_2, r_3=r_3,
					 rho=rho, alpha_1=alpha_1, alpha_2=alpha_2, nu=nu, cond_pol=cond_pol)

tik_d <- function(time, state, prams){
	with(as.list(c(state,prams)),{
		lam_u <- bet_u*((U_l+eta*U_h)/N) #forces of infection
		lam_s <- bet_s*D/N
		
		dS <- mu*N - (lam_u+lam_s)*S - mu*S
		dU_l <- (lam_u+lam_s)*S + rho*U_h + r_3*Q - (mu+k+theta_1)*U_l
		dU_h <- k*U_l + r_1*U_r + r_2*Q - (mu+sig+rho)*U_h
		dU_r <- sig*U_h - (mu+r_1+theta_2)*U_r
		dQ <- theta_1*U_l + theta_2*U_r - (mu+r_2+r_3)*Q
		dD <- alpha_1*(U_l + phi*U_h) - nu*D
		
		return(list(c(dS,dU_l,dU_h,dU_r,dQ,dD)))
	})
}

#nondim model with further conditional policing based on w threshold
tik_nd <- function(time, state, prams){
	with(as.list(c(state,prams)),{
		if (w<0.02){
			ds <- mu - beta_u*(v+eta*w)*s - beta_hats*z*s - mu*s
			dv <- beta_u*(v+eta*w)*s + beta_hats*z*s + rho*w + r_3*y - (mu+k+theta_1)*v
			dw <- k*v + r_1*x + r_2*y - (mu + sig* + rho)*w
			dx <- sig*w - (mu+r_1+theta_2)*x
			dy <- theta_1*v + theta_2*x - (mu + r_2 + r_3)*y
			dz <- (nu)*(v + phi*w - z)
		}else{
			ds <- mu - beta_u*(v+eta*w)*s - beta_hats*z*s - mu*s
			dv <- beta_u*(v+eta*w)*s + beta_hats*z*s + rho*w + r_3*y - (mu+k+theta_1)*v
			dw <- k*v + r_1*x + r_2*y - (mu + sig + rho)*w
			dx <- sig*w - (mu+r_1+theta_2)*x
			dy <- theta_1*v + theta_2*x - (mu + r_2 + r_3)*y
			dz <- (nu)*(v + phi*w - z*cond_pol)
		}
		return(list(c(ds,dv,dw,dx,dy,dz)))
	})
}

#The 5x5 F from Nyabadza et al., excludes the S equation/state var
F1 <- matrix(data = c(beta_u, eta*beta_u, 0, 0, beta_hats,
												 0,0,0,0,0,
												 0,0,0,0,0,
												 0,0,0,0,0,
												 0,0,0,0,0), nrow=5,ncol=5,byrow=TRUE)
#V also from Nyabadza et al., has quitter terms bc not treated
#as new infections, despite van den Driessche/Watmough ruling on the matter
V1 <- matrix(data = c(Q_1, -rho, 0, -r_3, 0,
										 -k, Q_2, -r_1, -r_2, 0,
										 0, -sig, Q_3, 0, 0,
										 -theta_1, 0, -theta_2, Q_4, 0,
										 -nu, -phi*nu, 0, 0, nu), nrow=5,ncol=5,byrow=TRUE)

#My F and V for next gen matrix. Does not produce meaningful R_0 results
#with the parameters determined by Nyabadza et al., but they also 
#determined their parameter ranges based on their own expression for R_0
F2 <- matrix(data = c(beta_u, eta*beta_u, 0,
											0,0,0,
											0,0,0), nrow=3, byrow=TRUE)
V2 <- matrix(data = c(Q_1, -rho, 0,
											-k, Q_2, -r_1,
											0, -sig, Q_3), nrow=3, byrow=TRUE)

#find FV^-1
n_gen1 <- F1 %*% inv(V1)
n_gen2 <- F2 %*% inv(V2)
#find eigs and R_0
eigs1 <- eigen(n_gen1)
eigs2 <- eigen(n_gen2)
R_0 <- eigs1$values[1]

slv <- as.data.frame(ode(y=state, times=time, 
												 func=tik_nd,parms=prams))

plot(s~time, slv, main = "Tik_nd",ylim = c(0,1.0), 
		 col = 'blue', xlab="Time (months)",ylab="Fraction of Pop. Density", type="l")
lines(v~time, slv, col="red", type="l")
lines(w~time, slv, col="green2", type="l")
lines(x~time, slv, col="purple3", type="l")
lines(y~time, slv, col="orange2", type="l")
lines(z~time, slv, col="yellow3", type="l")
legend("topright", c("s", "v", "w", "x", "y", "z"), 
			 col=c("blue", "red", "green2", "purple3","orange2","yellow3"), lty=1)