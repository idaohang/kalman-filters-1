function [x_post,P_post]=UKF(fstate,hmeas,x,P,z,Q,R)
% UKF   Unscented Kalman Filter 
%         
% Inputs:   f: function handle for f(x) -- process  update
%           h: function handle for h(x) -- measurements update
%           x: "a priori" state estimate
%           P: "a priori" estimated state covariance
%           z: current state measurement
%           Q: process noise covariance 
%           R: measurement noise covariance

% Output:   x: "a posteriori" state estimate
%           P: "a posteriori" state covariance


%Latent Variable and Observed Variables
d=numel(x);                                 %dimension(number) of latent state space
m=numel(z);                                 %dimension(number) of observed state space

%Tunable Scaling Parameters for Sigma Sampling Points
alpha=1e-3;                                 %tunable
kappa=0;                                    %tunable
beta=2;                                     %tunable

[X, wM, wC] = SigmaPoints(x, P, alpha, beta, kappa); %Generating Sigma Points from Prior Gaussian
                           
[x1,X1,P1,X_c]=UT(fstate,X,wM,wC,Q);       %Unscented Transformation of Process

[z1,~,P2,Z_c]=UT(hmeas,X1,wM,wC,R);      %Unscented Transformation of Measurements

P12=X_c*diag(Wc)*Z_c';                    %Transformed Cross-Covariance

%K = P12*inv(P2);   %Kalman gain
K = P12/P2; %Efficient Computation of Kalman Gain
x_post=x1+K*(z-z1);                              %Latent State Update
P_post=P1-K*P12*K';                                %Latent State Covariance update

function [mu, xPts, P, xPts_devs]=UT(f, X, wM, wC, Q)
%Unscented Transformation
%Input:
%        f: nonlinear map
%        X: sigma points
%       wM: weights for mean
%       wC: weights for covariance
%        Q: additive covariance
%Output:
%        mu: transformed mean
%        trans_sig: transformed sampling points
%        Y_c: transformed deviations
%        P: transformed covariance

num_points = size(X,2);
d = size(X, 1);
mu=zeros(d, 1);
xPts=zeros(d, num_points);

for k=1:num_points                   
    xPts(:,k)=f(X(:,k));       
    mu=mu+wM(k)*xPts(:,k);       
end
xPts_devs = xPts-repmat(mu, 1, d);
P=xPts_devs*diag(wC)*xPts_devs'+Q;          
