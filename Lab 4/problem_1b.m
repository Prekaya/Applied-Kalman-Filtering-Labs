clc;
dt = 0.1;   % Sim time step, s
tMax = 10;  % Max time, s
q = 2;

var_init_uncertainty = [10; 1]; % nx1 (units of states)
var_meas_noise = 10;

A = [1 dt; 0 1]; % nxn (function of xhat, not x)
Q = [dt^3/3, dt^2/2; dt^2/2, dt].*q;  % nxn

H = [0 1];
R = diag(var_meas_noise); % mxm

P = diag(var_init_uncertainty); % nxn (Initial P)

tHistory = []; % Time vector
pHistory = []; % P diagonal elements (n columns)

% Time Loop
t=0;
while t<tMax
    K = P*H'/(H*P*H' + R);
    I = eye(size(P));
    P = (I - K*H)*P;
    P = real(.5*P + .5*P'); % Make sure P stays real and symmetric
    tHistory(end+1,:) = t';
    pHistory(end+1,:) = diag(P)';

    t = t + dt;
    P = A*P*A' + Q;
    
end

tcl = tiledlayout(2,1);

nexttile;
plot(tHistory,sqrt(pHistory(:,1)),'.-');
ylim([0.5 2])
ylabel('\surdP(1,1)')

nexttile;
plot(tHistory,sqrt(pHistory(:,2)),'.-');
% ylim([0.95 1.15])
ylabel('\surdP(2,2)')
xlabel(tcl,'Times (s)')
title(tcl, 'Covariance Plots Position')