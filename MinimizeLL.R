# Main function  
MinimizeLL = function()  {
CIR200=matrix(,200,8)
tau = c(3/12,6/12,2,5);
for (i in 1:200) {

# Simulating zero coupon rates Ya
Ya = SimulateCIR(0.10,0.05,0.075,-0.4,1/12,0.06,120,tau);
Y=Ya;
nrow = dim(Y)[1];
ncol= dim(Y)[2];

para=c(0.10,0.05,0.075,-0.4,0.1*rnorm(ncol))
# Minimizing negative loglikelihood of CIR parameters given Ya

x=nlm(LogLikelihoodCIR,para,Ya, tau,nrow,ncol) 
CIR200[i,]=x$estimate

}
mean(CIR200)
sd(CIR200)

}

