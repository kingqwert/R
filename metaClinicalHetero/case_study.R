library(meta)
library(scatterplot3d)

# Data in the application section
data = data.frame(
	author = seq(1,24),
	year = rep(2000,24),
	event.e =c(3,2,7,1,0,0,2,1,6,1,0,7,8,1,11,7,4,1,6,7,3,6,12,2),
	n.e =c(57,9,27,15,7,15,9,14,9,37,10,32,15,7,40,34,14,29,34,18,12,20,116,16),
	event.c = c(3,3,0,1,2,1,0,0,6,0,1,5,4,2,3,3,1,0,4,6,2,4,6,1),
	n.c = c(84,11,25,14,17,15,9,6,8,32,10,31,16,7,39,36,13,30,33,22,12,20,103,15) 
)

##################################################
# Main part of our method for the application section
##################################################

# Run the conventional meta-analysis with REML method
sim = metabin(event.e, n.e, event.c, n.c,
	data=data, 
	sm="RR", method="I",
	studlab=paste(author))

# Forest plot for Figure 1
forest(sim)

# Function for our propose method
re.reml = function(y,sigma.i,a,b,par){
	sigmaREREML =par[1]
	beta = par[2]

	0.5*(sum(log(sigmaREREML^2+sigma.i)+(y-beta)^2/(sigmaREREML^2+sigma.i))+log(sum((sigmaREREML^2+sigma.i)^(-1))))-a*exp(-length(y)^2/b)*log(sigmaREREML^2)
}

# Optimize the likelihood
# this$par[1]: estimated between-study variance
# this$par[2]: estimated overall effect
this = optim(par=c(0.1,1),fn=function(r){re.reml(y=sim$TE,sigma.i=sim$seTE^2,a=0.1,b=1000,par=r)},method="BFGS")

# Calculate Boundary criteria defined in the expression (8)
BC = sum((sim$TE-sim$TE.fixed)^2/sim$seTE^2)-sum(1/sim$seTE)+(sum(sim$seTE^(-2)))/(sum(sim$seTE^(-1)))

# 95% CI
upper.ci = exp(this$par[2]-qnorm(0.025)*sqrt(1/sum(1/(sim$seTE^2+this$par[1]^2))))
lower.ci = exp(this$par[2]+qnorm(0.025)*sqrt(1/sum(1/(sim$seTE^2+this$par[1]^2))))

# Calculate I^2 and H^2 indices for heterogeneity
vi.avg = (sim$k-1) / (sum(1/sim$seTE) - sum(1/sim$seTE^2)/sum(1/sim$seTE))
I2 = 100 * this$par[1]^2 / (vi.avg + this$par[1]^2)
H2 = this$par[1]^2 / vi.avg + 1


##################################################
# Sensitivity check by changing the parameters a,b
##################################################
result_beta = matrix(NA,nrow=100,ncol=100)
result_sigma = matrix(NA,nrow=100,ncol=100)

a = seq(0.01,1,by=0.01)
b = seq(100,10000,length=100)


for(i in 1:100){
	for(j in 1:100){
		ai = a[i]
		bi = b[j]
		
opt = optim(par=c(0.1,1),fn=function(r){re.reml(y=sim$TE,sigma.i=sim$seTE^2,a=ai,b=bi,par=r)},method="BFGS")$par

		result_beta[i,j]=exp(opt[2])
		result_sigma[i,j]=opt[1]
	}
}
