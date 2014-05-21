
%First Attempts at a Run File

q=0.1;    %std of process 
r=0.1;    %std of measurement
Q=q^2*eye(n); % covariance of process
R=r^2;        % covariance of measurement  
f=@(x)[x(2);x(3);0.05*x(1)*(x(2)+x(3))];  % nonlinear state equations
h=@(x)x(1);                               % measurement equation
s=[0;0;1];                                % initial state
x=s+q*randn(3,1); %initial state          % initial state with noise
P = eye(n);                               % initial state covariance
N=20;                                     % total dynamic steps
xV = zeros(len(x),N);          %estmate        % allocate memory
sV = zeros(len(x),N);          %actual
zV = zeros(1,N);
for k=1:N
  z = h(s) + r*randn;                     % measurements
  sV(:,k)= s;                             % save actual state
  zV(k)  = z;                             % save measurment
  [x, P] = ukf(f,x,P,h,z,Q,R);            % ekf 
  xV(:,k) = x;                            % save estimate
  s = f(s) + q*randn(3,1);                % update process 
end



%{


%x = [X, Y, alpha, phi, v]

x = [.1; .2; .3; .4; 1]; 
P = eye(5);
z = [.45; 1.1];

alpha = 1;
beta = 0;
kappa = 2;
