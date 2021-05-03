library(deSolve)
library(matlib)
library(pracma)
N = 5000000

mu = 0.02
beta_u = 0.01
beta_s = 0.15
eta = 0.8
k = 0.3 #esc. from U_l to U_h
##Make rehab rate a function of policing, fear/spite toward law
sig = 0.8 #uptake rate into rehab, sigma
theta_1 = 0.1 #about 6 or 7 months of use for 1.8
theta_2 = 0.1 ##from U_l,U_r into Q
###1/theta_i = average duration of use, fraction of 12 months
r_1 = 0.01
r_2 = 0.001
r_3 = 0.001 #relapse rates
rho = 0.2 #amelioration U_h to U_l
alpha_1 = 0.6
phi = 1.4 #phi>1
alpha_2 = phi*alpha_1 #p2s contact rates, U_l and U_h
nu = 0.75 #1/nu = avg duration of drugs in circulation / supply chain
beta_hats <- beta_s*alpha_1/nu #transmission pram, S2P2S transmission cycle
epsilon = 0.001 #timescale parameter

cond_pol = 1.5 #conditional policing, if many heavy users, bigger crackdown policy 
z_max = 0.35
k_m = 0.18

beta_hats = beta_s*alpha_1/nu #NOT in prams yet
Q_1 = mu + k + theta_1
Q_2 = mu + sig + rho
Q_3 = mu + r_1 + theta_2
Q_4 = mu + r_2 + r_3
PHI_1 = (sig*r_1)/(Q_2*Q_3)
PHI_2 = (k*rho)/(Q_1*Q_2)
PHI_3 = (sig*r_2*theta_2)/(Q_2*Q_3*Q_4)
PHI_4 = (rho*r_2*theta_1)/(Q_1*Q_2*Q_4)
PHI_5 = (r_3*theta_1)/(Q_1*Q_4)
PHI_6 = (k*sig*r_3*theta_2)/(Q_1*Q_2*Q_3*Q_4)
PI = 1/(1-(PHI_1+PHI_2+PHI_3+PHI_4+PHI_5*(1-PHI_1)+PHI_6))
R_1 = (beta_u/Q_1)*((1-PHI_1-PHI_3)/(1-(PHI_1+PHI_2+PHI_3+PHI_4)))
R_2 = (eta*beta_u/Q_2)*((k/Q_1)+(r_2*theta_1)/(Q_1*Q_4))*PI
R_3 = ((phi*beta_hats/Q_2)*((k/Q_1) + (r_2*theta_1/(Q_1*Q_4))) +
         (beta_hats/Q_1)*(1-PHI_1-PHI_3))*PI
R = R_1 + R_2 + R_3
xi = (k*Q_4+r_2*theta_1)/(Q_2*Q_4*(1-(PHI_1 + PHI_3)))

s_star = Q_1*((1-(PHI_1+PHI_2+PHI_3+PHI_4+PHI_5*(1-PHI_1)+PHI_6))/
                (1-(PHI_1+PHI_3))*(beta_u*(1+eta*xi)+beta_hats*(1+phi*xi)))


time <- seq(0, 250, 1)
state <- c(s=0.951,v=0.04,w=0.005,x=0.004,y=0.000,z=0.05)
prams <- c(mu=mu, beta_u=beta_u, beta_s=beta_s, eta=eta, k=k, sig=sig, 
					 theta_1=theta_1, theta_2=theta_2, r_1=r_1, r_2=r_2, r_3=r_3,
					 rho=rho, phi=phi, alpha_1=alpha_1, alpha_2=alpha_2, nu=nu, 
					 beta_hats=beta_hats, cond_pol=cond_pol, z_max=z_max)

# tik_nd <- function(time, state, prams){
#   with(as.list(c(state,prams)),{
#       ds <- mu - beta_u*(v+eta*w)*s - beta_hats*z*s - mu*s
#       dv <- beta_u*(v+eta*w)*s + beta_hats*z*s + rho*w + r_3*y - (mu+k+theta_1)*v
#       dw <- k*v + r_1*x + r_2*y - (mu + sig + rho)*w
#       dx <- sig*w - (mu+r_1+theta_2)*x
#       dy <- theta_1*v + theta_2*x - (mu + r_2 + r_3)*y
#       dz <- nu*(v + phi*w) - nu*z #try K_m=0.04 for bifurcation w/ z
#     return(list(c(ds,dv,dw,dx,dy,dz)))
#   })
# }

#nondim model with further conditional policing based on w threshold
tik_nd <- function(time, state, prams){
	with(as.list(c(state,prams)),{
			ds <- mu - beta_u*(v+eta*w)*s - beta_hats*z*s - mu*s
			dv <- beta_u*(v+eta*w)*s + beta_hats*z*s + rho*w + r_3*y - (mu+k+theta_1)*v
			dw <- k*v + r_1*x + r_2*y - (mu + sig + rho)*w
			dx <- sig*w - (mu+r_1+theta_2)*x
			dy <- theta_1*v + theta_2*x - (mu + r_2 + r_3)*y
			dz <- nu*(v + phi*w) - (z_max*z)/(k_m+z) #try K_m=0.04 for bifurcation w/ z
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
########## For NGM with new saturating policing
F1_sat <- matrix(data = c(beta_u, eta*beta_u, 0, 0, beta_hats,
                      0,0,0,0,0,
                      0,0,0,0,0,
                      0,0,0,0,0,
                      0,0,0,0,0), nrow=5,ncol=5,byrow=TRUE)

V1_sat <- matrix(data = c(Q_1, -rho, 0, -r_3, 0,
                      -k, Q_2, -r_1, -r_2, 0,
                      0, -sig, Q_3, 0, 0,
                      -theta_1, 0, -theta_2, Q_4, 0,
                      -nu, -phi*nu, 0, 0, z_max/k_m), nrow=5,ncol=5,byrow=TRUE)

#My F and V for next gen matrix. Does not produce meaningful R_0 results
#with the parameters determined by Nyabadza et al.
F2 <- matrix(data = c(beta_u, eta*beta_u, 0, beta_hats,
											0,0,0,0,
											0,0,0,0,
											nu, nu*phi,0,0), nrow=4, byrow=TRUE)
V2 <- matrix(data = c(Q_1, -rho, 0, 0,
											-k, Q_2, -r_1, 0,
											0, -sig, Q_3, 0,
											0, 0, 0, nu), nrow=4, byrow=TRUE)

#find FV^-1
n_gen1 <- F1 %*% inv(V1)
n_gen2 <- F2 %*% inv(V2)
#find eigs and R_0
eigs1 <- eigen(n_gen1)
eigs2 <- eigen(n_gen2)
R_0 <- eigs1$values[1]

slv <- as.data.frame(ode(y=state, times=time, 
												 func=tik_nd,parms=prams))

###Stability of DFE
sys <- function(x){ #vector value function of nondim system
  return(c(mu - beta_u*(x[2]+eta*x[3])*x[1] - beta_hats*x[6]*x[1] - mu*x[1],
           beta_u*(x[2]+eta*x[3])*x[1] + beta_hats*x[6]*x[1] + rho*x[3] + r_3*x[5] - (mu+k+theta_1)*x[2],
           k*x[2] + r_1*x[4] + r_2*x[5] - (mu + sig + rho)*x[3],
           sig*x[3] - (mu+r_1+theta_2)*x[4],
           theta_1*x[2] + theta_2*x[4] - (mu + r_2 + r_3)*x[5],
           nu*(x[2] + phi*x[3] - x[6]))
         )
}
jac <- jacobian(f=sys, x0=c(1,0,0,0,0,0)) #jacobian eval at disease-free equilibrium
jac_eigs <- eigen(jac)

####Plotting
plot(s~time, slv, main = "Tik_nd", ylim = c(0,1), 
		 col = 'blue', xlab="Time",ylab="Fraction of Pop. Density", type="l")
lines(v~time, slv, col="red", type="l")
lines(w~time, slv, col="green2", type="l")
lines(x~time, slv, col="purple3", type="l")
lines(y~time, slv, col="orange2", type="l")
lines(z~time, slv, col="yellow3", type="l")
legend("topright", c("s", "v", "w", "x", "y", "z"), 
			 col=c("blue", "red", "green2", "purple3","orange2","yellow3"), lty=1)

### quasi-steady state z eqn for endemic equilibrium
z <- function(t){
	return((state[6]-sapply(slv,tail,1)[3]*(1+phi*xi))*exp(-nu*t)+sapply(slv,tail,1)[3]*(1+phi*xi))
}
#plot z with current parameters
#curve(z, xlab = "Time", ylab = "z(t)",add=FALSE,xlim=c(0,100),type ="l")


##################### Contour Plotting

# step_size <- 15 #this whole method is wrong bc forgot abt Q_1 and Q_3
# Rs <- matrix(nrow=step_size,ncol=step_size) #matrix to contain R_0s for contour plotting
# F1_copy <- F1
# V1_copy <- V1
# for(i in 1:step_size){
#   V1_copy[4,1] <- -i/step_size #update -theta_1
#   for(j in 1:step_size){
#     V1_copy[4,3] <- -j/step_size #"" -theta_2
#     temp_eigs <- eigen(F1_copy %*% inv(V1_copy))
#     Rs[i,j] <- max(temp_eigs$values) #R_0 for that pair of thetas
#   }
# }

library(ContourFunctions)
F1_copy <- F1
V1_copy <- V1
# R_contour_func1 <- function(x){ #this one uses bad formula from paper
#   Q_1 = mu + k + x[1]
#   Q_3 = mu + r_1 + x[2]
#   PHI_1 = (sig*r_1)/(Q_2*Q_3)
#   PHI_2 = (k*rho)/(Q_1*Q_2)
#   PHI_3 = (sig*r_2*x[2])/(Q_2*Q_3*Q_4)
#   PHI_4 = (rho*r_2*x[1])/(Q_1*Q_2*Q_4)
#   PHI_5 = (r_3*x[1])/(Q_1*Q_4)
#   PHI_6 = (k*sig*r_3*x[2])/(Q_1*Q_2*Q_3*Q_4)
#   PI = 1/(1-(PHI_1+PHI_2+PHI_3+PHI_4+PHI_5*(1-PHI_1)+PHI_6))
#   return((beta_u/Q_1)*((1-PHI_1-PHI_3)/(1-(PHI_1+PHI_2+PHI_3+PHI_4))) +
#          (eta*beta_u/Q_2)*((k/Q_1)+(r_2*x[1])/(Q_1*Q_4))*PI +
#          ((phi*beta_hats/Q_2)*((k/Q_1) + (r_2*x[1]/(Q_1*Q_4))) +
#                   (beta_hats/Q_1)*(1-PHI_1-PHI_3))*PI)
# }

R_contour_func2 <- function(x){ #this one uses NGM method with good results
  V1_copy[4,1] <- -x[1] #update theta_1
  V1_copy[4,3] <- -x[2] #update theta_2
  V1_copy[1,1] <- mu + k + x[1] #update Q_1
  V1_copy[3,3] <- mu + r_1 + x[2] #update Q_3
  temp_ngm_eigs <- eigen(F1_copy %*% inv(V1_copy))
  return(temp_ngm_eigs$values[1])
}

library(ggplot2)
#these plot x[2] against x[1]

plot1 <- cf(R_contour_func2,
            main="(Proportional Policing, r_2 = 0.001, r_3 = 0.001)", 
            color.palette=heat.colors,bar=T, with_lines=T,
            cex.main=1.5)

# F1_copy <- F1
# V1_copy <- V1
# V1_copy[5,5] <- z_max/k_m
# 
# plot_sat <- cf(R_contour_func2, 
#                main="R_0 as function of Theta_2 vs. Theta_1 (Saturating Policing)",
#                color.palette=heat.colors, xlab="theta_1", bar=T, with_lines=T,
#                cex.main=1.5)

# F1_copy <- F1
# V1_copy <- V1
# V1_copy[2,4]<- V1_copy[1,4] <- 0
# V1_copy[4,4] <- mu
# plot_no_relapse_reg <- cf(R_contour_func2,
#                           main="(Proportional Policing, no relapse)", 
#                           color.palette=heat.colors,bar=T, with_lines=T,
#                           cex.main=1.2)
# 
# F1_copy <- F1
# V1_copy <- V1
# V1_copy[2,4]<- V1_copy[1,4] <- 0
# V1_copy[4,4] <- mu
# V1_copy[5,5] <- z_max/k_m
# plot_no_relapse_sat <- cf(R_contour_func2,
#                           main="(Saturating Policing, no relapse)", 
#                           color.palette=heat.colors,bar=T, with_lines=T,
#                           cex.main=1.2)