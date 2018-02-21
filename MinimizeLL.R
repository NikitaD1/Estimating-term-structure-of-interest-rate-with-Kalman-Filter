MinimizeLL = function()  {

clear all
for (i in 1:50) {
# Inputs for term structure generator, and minimizer
tau=c(3/12,6/12,2,5);
Ya=SimulateCIR(0.10,0.05,0.075,-0.4,1/12,0.06,120,tau);
Y=Ya;
nrow= dim(Y)[1];
ncol= dim(Y)[2];
#lb=[0.0001,0.0001,0.0001,-1, 0.00001*ones(1,ncol)];
#ub=[ones(1,8)];

para =c(0.10,0.05,0.075,-0.4,0.1*rnorm(ncol))
# Minimizer
CIR200=matrix(,200,8)
x=nlm(LogLikelihoodCIR,para) #,[],[],[],[],lb,ub,[],[],Y, tau, nrow, ncol)
CIR200[i,]=x$estimate

}
mean(CIR200)
sd(CIR200)

}
