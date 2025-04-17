clear all
close all 
clc

% problem parameters
nSamples = 100; 
m0 = 0; 
P0 = 4;
PHI = 1; 
Q = 0; 
H = 1; 
R = 0.25;
[p,q] = size(PHI);   
[r,~] = size(R);
 
% simulate random process and generate sensor measurements
%rng('default')
randn('seed',0)
sqrtP0 = sqrt(P0); 
sqrtR = sqrt(R); 
sqrtQ = 0;
x = m0+sqrtP0*randn(p,1);
z = H*x+sqrtR*randn(r,1);
X = zeros(p,nSamples);  
Z = zeros(r,nSamples);
X(:,1) = x;  Z(:,1) = z;
for i=2:1:nSamples
    x = PHI*x+sqrtQ*randn(p,1);
    z = H*x+sqrtR*randn(r,1);
    X(:,i) = x;  Z(:,i) = z;
end

% plot simulated state and measurements
figure(1)
hold on
plot(1:nSamples, X, 'b','linewidth',2)
plot(1:nSamples, Z, 'r*')
xlabel('Samples (k)')
ylabel('State (x_1)')
legend('state (x_1)','measurement (z)')
set(gca,'fontsize',14)




