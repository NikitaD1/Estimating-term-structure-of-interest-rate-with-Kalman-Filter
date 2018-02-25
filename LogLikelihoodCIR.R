LogLikelihoodCIR = function(para,Y, tau, nrow, ncol) {
# This function returns the log likelihood value of CIR parameters given the zero coupon rates Y

theta=para[1]; 
kappa=para[2];
sigma=para[3];
lambda=para[4];
sigmai=para[5:length(para)];

R=diag(ncol);
for (i in 1:ncol)              
R[i,i]=sigmai[i]^2; # R is measurement noise

dt=1/12;

# System Matrices Initialization
C=theta*(1-exp(-kappa*dt));   # Eq 63. C & F are co-efficent matrix of transition equation
F=exp(-kappa*dt);                   
A=as.matrix(t(rep(0,ncol)));  # A & H are co-efficient matrix of measurement equation                   
H=as.matrix(A);                                   

# Create A and B
for (j in 1:ncol) {                                     
  Gamma=sqrt((kappa+lambda)^2+2*sigma^2);       # Eq 29                    
  CIRB=2*(exp(Gamma*tau[j])-1)/((Gamma+kappa+lambda)*(exp(Gamma*tau[j])-1)+2*Gamma);                            
  CIRA=2*kappa*theta/(sigma^2)*log(2*Gamma*exp((Gamma+kappa+lambda)*tau[j]/2)/((Gamma+kappa+lambda)*(exp(Gamma*tau[j])-1)+2*Gamma));
  A=-CIRA/tau[j]       
  B=CIRB/tau[j]   # Eq. 65
}

# Kalman Filter 
# Step 1 Initializing the state vector

AdjS=as.matrix(theta);   # Unconditional mean of state variable 
VarS=as.matrix(sigma^2*theta/(2*kappa));  # Unconditional variance
LL=t(rep(0,nrow)); # Intializing likelihood vector
for (i in 1:nrow) {
  
# Updating state variable & its variance with transition eq.  
PredS=C+F*AdjS;   #Eq 75
Q=theta*sigma*sigma*(1-exp(-kappa*dt))^2/(2*kappa)+sigma*sigma/kappa*(exp(-kappa*dt)-exp(-2*kappa*dt))*AdjS;    #Eq 63
VarS= F*VarS*F+ Q;                             # eqn 76

# Step 2 Forecasting the Measurement equation
PredY= A+H*c(PredS);              # Eq 69. Conditional Forecast of the measurement equation 
VarY=(t(H)* c(VarS))%*%H+R;             # Eq 70. Associated Conditional Variance eqn 
# Step 3 Updating the inference about the state vector
PredError=as.matrix(Y[i,]-PredY);                       # Eq 71. Prediction error on z 
KalmanGain=(c(VarS)*H)%*% solve(VarY);                   # Eq 73. Kalman Gain Matrix 
AdjS=PredS+KalmanGain%*%t(PredError);             # Eq 72. Update inference about the unodserved transition system eqn c.7
VarS=c(VarS) %*% (c(1)-KalmanGain%*%c(t(H)))                  # Eq 74 Update conditional variance eqn c.9

# Step 5 Construct the likelihood function
DetY=det(VarY);
LL[i]=-(ncol/2)*log(2*pi)-0.5*log(DetY)-0.5*PredError%*%solve(VarY)%*%t(PredError); # eqn 77
}
sumll=-sum(LL);
return(sumll)
}
