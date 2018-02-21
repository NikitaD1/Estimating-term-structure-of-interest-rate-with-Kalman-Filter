LogLikelihoodCIR = function(para,Y, tau, nrow, ncol) {
# initialize the parameters for CIR model
tau=c(3/12,6/12,2,5);
ncol=4;
para= c(0.10,0.05,0.075,-0.4,0.1*rnorm(ncol));
Ya=SimulateCIR(0.10,0.05,0.075,-0.4,1/12,0.06,120,tau);
Y=Ya;
nrow= dim(Y)[1];


lb=c(0.0001,0.0001,0.0001,-1, 0.00001*rep(1,ncol));
ub=rep(1,8);
para =c(0.10,0.05,0.075,-0.4,0.1*rnorm(ncol))

theta=para[1]; 
kappa=para[2];
sigma=para[3];
lambda=para[4];
sigmai=para[5:length(para)];

R=diag(ncol);
for (i in 1:ncol)              
R[i,i]=sigmai[i]^2;

dt=1/12;

# System Matrices Initialization
C=theta*(1-exp(-kappa*dt));           # eqn b.3
F=exp(-kappa*dt);                     # eqn b.3
A=as.matrix(t(rep(0,ncol)));                      
H=as.matrix(A);                                   

# Create A and B
for (i in 1:ncol) {                  # System matrices are made for each tau                    
AffineG=sqrt((kappa+lambda)^2+2*sigma^2);                           # eqn a.10
AffineB=2*(exp(AffineG*tau[i])-1)/((AffineG+kappa+lambda)
                                   *(exp(AffineG*tau[i])-1)+2*AffineG);                           # eqn a.9
AffineA=2*kappa*theta/(sigma^2)*log(2*AffineG*
                                    exp((AffineG+kappa+lambda)*tau[i]/2)/((AffineG+kappa+lambda)*
                                                                          (exp(AffineG*tau[i])-1)+2*AffineG));                            # eqn a.8
A[i]=-AffineA/tau[i];       # eqn b.1
H[i]=AffineB/tau[i];       # eqn b.1
}

# Kalman Filter 
# Step 1 Initializing the state vector
initx=theta;                       # Unconditional mean eqn c.1  
initV=sigma^2*theta/(2*kappa);    # Unconditional variance eqn c.3
# Starting values
AdjS=as.matrix(initx);
VarS=as.matrix(initV);
LL=t(rep(0,nrow));
for (i in 1:nrow) {
PredS=C+F*AdjS;                              # eqn c.10
Q=theta*sigma*sigma*(1-exp(-kappa*dt))^2/(2*kappa)+sigma*sigma/kappa*(exp(-kappa*dt)-exp(-2*kappa*dt))*AdjS;    # eqn b.4
VarS= F*VarS*F+ Q;                             # eqn c.11
# Step 2 Forecasting the Measurement equation
PredY= A+H*c(PredS);              # Conditional Forecast of the measurement equation eqn c.4
VarY=(t(H)* c(VarS))%*%H+R;             # Associated Conditional Variance eqn c.5
# Step 3 Updating the inference about the state vector
PredError=as.matrix(Y[i,]-PredY);                       # Prediction error on z eqn c.6
KalmanGain=(c(VarS)*H)%*% solve(VarY);                   # Kalman Gain Matrix eqn c.8
AdjS=PredS+KalmanGain%*%t(PredError);             # Update inference about the unodserved transition system eqn c.7
VarS=c(VarS) %*% (c(1)-KalmanGain%*%c(t(H)))                  # Update conditional variance eqn c.9

# Step 5 Construct the likelihood function
DetY=det(VarY);
LL[i]=-(ncol/2)*log(2*pi)-0.5*log(DetY)-0.5*PredError%*%solve(VarY)%*%t(PredError); # eqn c.12
}
sumll=-sum(LL);
return(sumll)
}
