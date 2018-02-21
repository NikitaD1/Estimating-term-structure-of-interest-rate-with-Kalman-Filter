SimulateCIR = function(theta,kappa,sigma,lambda,dt,ratestart,months,tau) {

#CIR term structure simulation, inputs below


theta=0.10;
kappa=0.05;
sigma=0.075;
lambda=-0.4;
dt=1/12;
months=120;
ratestart=0.10;
tau= c(3/12,6/12,2,5);
srt <- vector(mode="numeric", length=480)
Rt = matrix(,120,length(tau))
i = 1 
 # Short Rate Dynamics
srt[1]=ratestart;
for (i in 1:480)

srt[i+1]=srt[i]+kappa*(theta-srt[i])*dt/4+sqrt(srt[i])*sigma*sqrt(dt/4)*rnorm(1);

# Term Structure Dynamics
for (i in 1:months) {

rttemp=srt[i*4-3];
#rttemp1[i]=rt[i*4-3];
for (j in 1:length(tau)) {

AffineG=sqrt((kappa+lambda)^2+2*sigma^2);                           
AffineB=2*(exp(AffineG*tau[j])-1)/((AffineG+kappa+lambda)
                                   *(exp(AffineG*tau[j])-1)+2*AffineG);                            
AffineA=2*kappa*theta/(sigma^2)*log(2*AffineG*
                                    exp((AffineG+kappa+lambda)*tau[j]/2)/((AffineG+kappa+lambda)*
                                                                          (exp(AffineG*tau[j])-1)+2*AffineG));                            
A=-AffineA/tau[j];       
B=AffineB/tau[j];        
Rt[i,j]=A+B*rttemp;
}
}

#figure(2)
#surface3d(tau,1:120,Rt)
#persp(x,y,z, border= "black")
#persp(x,y,z,theta=30, phi=30, col=rainbow(25, start=0.5, end=0.8))
#persp(1:120,tau,Rt,theta=30, phi=30, expand=0.6, col='lightblue', shade=0.75, ltheta=120, ticktype='detailed')
#persp(1:120,tau,Rt,theta=30, phi=30, col=rainbow(25, start=0.5, end=0.8) )
#xlabel('Time to Maturity (Fixed)')
#ylabel('Months Passed')
title('A Given Simulation of the Term Structure (CIR)')
return(Rt)
}
