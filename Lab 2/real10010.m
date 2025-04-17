function [Z,X,PHI,Q] = generate_process(n,P0,m0,F,G,W,dt,H,R)
% generates n+1 data points from a stationary state space model
% x(k) = PHI*x(k-1)+w(k-1), z(k)=H*x(k)+v(k)
% where
% x(p,1), x(0)~N(m0,P0), w(q,1)~N(0,Q), v(r,1)~N(0,R)
% output
% X = [x(0), x(1), ..., x(n)] Z = [z(0), z(1), ..., z(n)]
% PHI and Q are computed from F, W, G in dx/dt=F*x+G*f
[p,q] = size(G); [r,~] = size(R);
X = zeros(p,n+1); Z = zeros(r,n+1);
% calculate PHI and Q
A = [-F G*W*G'; zeros(size(F)) F'];
B = expm(A.*dt);
B12 = B(1:p,p+1:2*p);
B22 = B(p+1:2*p,p+1:2*p);
PHI = B22';
Q = PHI*B12;
% compute the square roots of P0, Q, R
sqrtP0 = chol(P0)'; sqrtR = chol(R)'; sqrtQ = chol(Q)';
% generate random process
x = m0+sqrtP0*randn(p,1); z = H*x+sqrtR*randn(r,1);
X(:,1) = x; Z(:,1) = z;
for i=1:1:n
x = PHI*x+sqrtQ*randn(p,1);
z = H*x+sqrtR*randn(r,1);
X(:,i) = x; Z(:,i) = z;
end
end

% set seed value for random number generator
rng(0)
% number of samples
n = 100;
% define parameters for a Markov Process
% where the continuous form is:
% dx/dt = -beta*x+sqrt(2*sigs*beta)*f
sigs = 100;
beta = 1;
P0 = sigs; % initial covariance
m0 = 0; % initial mean
F = -beta; % state transition matrix
G = sqrt(2*sigs*beta); % distribution matrix
W = 1; % psd matrix
dt = .2; % time step
R = .0001; % measurement covariance
H = 1; % measurement model
% generate a random process (n+1 points)
figure(1)
hold on
for n = 1:100
    [Z,X,PHI,Q] = generate_process(n,P0,m0,F,G,W,dt,H,R);
    % plot the random process
    
    plot(0:1:n,Z(:),'*',0:1:n,X(1,:),'-')
    xlabel('samples'); ylabel('X,Z')
end