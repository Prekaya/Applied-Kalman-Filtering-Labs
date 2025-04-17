clc;clearvars;
dt = 1;   % Sim time step, s
tMax = 5;  % Max time, s

sd_init_uncertainty = [30; 30]; % nx1 (units of states)
sd_meas_noise = 30;

A = [0 1 0; 0 0 1; 0 0 0]; % nxn (function of xhat, not x)

H = [1 0; 0 1];
R = diag(sd_meas_noise.^2); % mxm

P_a = diag(sd_init_uncertainty.^2); % nxn (Initial P)

tHistory = []; % Time vector
pHistory = []; % P diagonal elements (n columns)
xHistory = []; % State vector

xhat = [0;0;0];

% Time Loop
t=0;
while t<tMax
    % Calculate posteriori
    K = P_a*H'/(H*P_a*H' + R);
    P_b = P_a - K*H*P_a;

    % record posteriori
    tHistory(end+1,:) = t';
    pHistory(end+1,:) = diag(P_b)';
    
    % Calculate priori
    t = t + dt;
    F = expm(A*t);
    xhat = F*xhat;

    P_a = F*P_b*F';
   
    % record priori
    tHistory(end+1,:) = t';
    pHistory(end+1,:) = diag(P_a)';

    
end

plot(tHistory,pHistory(:,1))