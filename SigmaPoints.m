function [xPts, wM, wC] = SigmaPoints(mu,P,alpha,beta,kappa)

% Returns scaled symmetric sigma point distribution (deterministic sample
% from MVN)
%	 mu:      MVN mean
%	 P:	     MVN covariance
%    alpha:     tunable scaling parameter 
%    beta:        tunable parameter to upweight on zeroth sample point
%	 kappa:	      tunable scaling parameter (usually default 0)
%
%    sig_pts sigma points
%    wM mean weights
%    wC covariance weights

% Number of sigma points and scaling terms
d = size(mu(:),1); %Dimension of the Space
N = 2*d+1;     %Number of Sigma Points

% Recalculate kappa according to scaling parameters
kappa_rs = alpha^2*(d+kappa)-d;

% Calculate matrix square root of weighted covariance matrix
Psqrtm=(chol((d+kappa_rs)*P))';  

% Array of the sigma points
xPts=[zeros(size(mu,1),1) -Psqrtm Psqrtm]+repmat(mu,1,N);

% Array of the mean/covariance weights for each sigma point 
wM=[kappa_rs 0.5*ones(1,N-1)]/(d+kappa_rs);
wC=[0 0.5*ones(1,N-1)]/(d+kappa_rs);
wC(1) = wM(1) + (1-alpha^2) + beta;



