SimulateCIR = function(theta,kappa,sigma,lambda,dt,ratestart,months,tau) {

#This function simulates one factor CIR short rate and returns the zero coupon rate
# CIR parameters taken are specified below
theta=0.10;
kappa=0.05;
sigma=0.075;
lambda=-0.4;
dt=1/12;
months=120; 
ratestart=0.10; #the first assumed short rate
tau= c(3/12,6/12,2,5); # tau = T- t , the time  until maturity
  
# Srt is the instantaneous short term interest rate(unobserved state variable) sampled every week for 10 years.
srt <- vector(mode="numeric", length=480)

# Rt is the monthly zero coupon rate (observed variable) observed every month for 10 years 
Rt = matrix(,120,length(tau))
i = 1 
 # Simulation of the shortrate following CIR model
srt[1]=ratestart;
for (i in 1:480)

srt[i+1]=srt[i]+kappa*(theta-srt[i])*dt/4+sqrt(srt[i])*sigma*sqrt(dt/4)*rnorm(1);

# Term Structure Dynamics
for (i in 1:months) {

rttemp=srt[i*4-3];
#rttemp1[i]=rt[i*4-3];
for (j in 1:length(tau)) {
 # P(t,T)= exp(A(t)-B(t)r(t)) 
Gamma=sqrt((kappa+lambda)^2+2*sigma^2);  #Eq 29                         
CIRB=2*(exp(Gamma*tau[j])-1)/((Gamma+kappa+lambda)*(exp(Gamma*tau[j])-1)+2*Gamma);                            
CIRA=2*kappa*theta/(sigma^2)*log(2*Gamma*exp((Gamma+kappa+lambda)*tau[j]/2)/((Gamma+kappa+lambda)*(exp(Gamma*tau[j])-1)+2*Gamma));
 # Measurement Equation for Kalman filter : Rt =  A + B*srt                                                                   
A=-CIRA/tau[j]       
B=CIRB/tau[j]        
Rt[i,j]=A+B*rttemp
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
#title('A Given Simulation of the Term Structure (CIR)')
return(Rt)
}
