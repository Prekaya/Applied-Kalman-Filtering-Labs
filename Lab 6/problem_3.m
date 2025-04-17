%% Problem 3
clc; close all;

% rng(0,"v4");
set(groot,'defaultTextInterpreter','latex') 

% SIM variables
tMax = 100;  % Max time, s
dt=1;
dtMeas=1;
tNextMeas = dtMeas;
meas_available = 0;

% EKF Parameters
var_init_uncertainty = [1; 1];
var_meas_noise = 0;
var_process_noise = [1; 0];
 
P = diag(var_init_uncertainty);% nxn Covariance matrix
R = diag(var_meas_noise); % mxm Measurement noise matrix
Q = diag(var_process_noise);

m = length(R);
n = length(P);
 
x0_true = [1; 0.9];
x0_est  = [0; 0];

H = [1, 0];

% For storing process and t=0 values
saveVars = {"T", "X_true", "X_est", "Z_true", "Z_est", "P_est", "P_plot", "step","P_lim","K_lim", "L_lim","info", "makeplot"};
xNames={'$\hat{X}$', '$\hat{\phi}$'};

T = 0:dt:tMax;
t_length = length(T);
X_true = nan(n,t_length); % True state vectors (n x steps)
X_est = nan(n,t_length); % Estimate state vectors (n x steps)
Z_true = nan(m,t_length); % True measurement vectors (m x steps)
Z_est = nan(m,t_length); % Estimate measurement vectors (m x steps)
P_est = nan(n,n,t_length); % Estimate variance vectors (n x n steps)
P_plot= nan(n,t_length);


xtrue = x0_true;
xest = x0_est;
ztrue = H*xtrue + sqrt(R)*randn(m,1);
zest = H*xest;

% A = eye(n) - F;
% [P_lim,K_lim,L_lim,info] = dare(A,H',Q,R,[],[]);

for i = 1:t_length
    % Store current time info
    t = T(i);
    X_true(:,i) = xtrue;
    X_est(:,i)  = xest;
    Z_true(:,i) = ztrue;
    Z_est(:,i) =  zest;
    P_est(:,:,i) = P;
    P_plot(:,i) = diag(P);
    
    % Propagate foward
    t = t+dt;

    % Sim update
    xtrue = [xtrue(2);1].*xtrue +  sqrt(Q)*randn(n,1);

    % Propagation of previous state estimate to current time    
    L = eye(n);
    Fest = jacobian_F(xest);
    P = Fest*P*Fest' + L*Q*L';
    xest = [xest(2);1].*xest;

    ztrue = H*xtrue + sqrt(R)*randn(m,1);
    zest = H*xest;

    if t>=tNextMeas
        tNextMeas = t+dtMeas;
        meas_available=1;
    else
        meas_available=0;
     end

    if meas_available == 1

        % Gain Matrix
        M = eye(m);
        K = (P*H')/(H*P*H' + M*R*M');

        % States updated with measurement information
        xest = xest + K*(ztrue - zest);

        % Covariance matrix updated with measurement information
        P = (eye(n) - K*H)*P;
    end
end
figure()
t = tiledlayout();
xlabel(t,'Time (s)')
title(t,'System Estimation')
nexttile
hold on
plot(T,X_true(2,:), 'k-');
plot(T,X_est(2,:), 'r-');
ylabel(xNames(2));

figure()
t = tiledlayout();
xlabel(t,'Time (s)')
title(t,'System Errors')
nexttile
hold on
plot(T,X_est(2,:)-X_true(2,:), 'r-');
plot(T, sqrt(P_plot(2,:)),'b-',T,- sqrt(P_plot(2,:)),'b-')
ylabel(xNames(2));

clearvars('-except',saveVars{:})

%-------------------------------------------------------------------------
function F = jacobian_F(x)
F = [x(2), x(1); 0, 1];
end
%-------------------------------------------------------------------------
function H = jacobian_h(x,N,E)
  
    H=[(x(1)-N(1))/sqrt((x(1)-N(1))^2 + (x(2)-E(1))^2), (x(2)-E(1))/sqrt((x(1)-N(1))^2 + (x(2)-E(1))^2), 0, 0;
       (x(1)-N(2))/sqrt((x(1)-N(2))^2 + (x(2)-E(2))^2), (x(2)-E(2))/sqrt((x(1)-N(2))^2 + (x(2)-E(2))^2), 0, 0];
end

%-------------------------------------------------------------------------
function z=h(x,N,E)
    z = [sqrt((x(1)-N(1))^2 + (x(2)-E(1))^2);
         sqrt((x(1)-N(2))^2 + (x(2)-E(2))^2)];
end