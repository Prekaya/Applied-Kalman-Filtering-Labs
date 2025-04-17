clc;clearvars

%true position and time bias being tested against
x_true = [1132049; -4903445; 3905453; 85000];

% measurements from satellites
measurements = [15764733, -1592675, 21244655;
                6057534, -17186958, 19396689;
                4436748, -25771174, 1546041;
                -9701586,-19687467, 15359118;
                23617496, -11899369, 1492340;
                14540070, -12201965, 18352632]; 

%generating pseudo distances using true measurement
sim_pseudo = vecnorm( measurements - x_true(1:3)',2,2) + x_true(4);
measurements(:,4) = sim_pseudo;

[estimates, vars] = find_estimates(measurements);

estimates

function [estimates, vars] = find_estimates(measurements)

    num_meas = size(measurements,1);

    sigma_MeasNoise_est = 25; % noise of satellite distances
    sigma_InitUncertainty_est = 10^18 * [1, 1, 1, 1]; % large uncertainty to simulate infinity

    R = sigma_MeasNoise_est*eye(num_meas); %mxm
    P = diag(sigma_InitUncertainty_est.^2); %nxn

    x_est=[0;0;0;0]; % States: [init_pos; init_vel; accel]
    z_meas = measurements(:,end); % simulated pseudo_range measurements

    tolerance = 10; % trying to get withing a 100m of true position
    distance = inf;
    i = 1;
    while distance > tolerance % looping till tolerance is reached

        % pseudo_range estimate based on current measurement estimates
        z_est = vecnorm(measurements(:,1:3) - x_est(1:3)',2,2) + x_est(4);     

        mags = vecnorm(x_est(1:3) - measurements(:,1:3)',2);
        
        % derivative of pseuudo_range formula to get measurement matrix
        H = [(x_est(1:3)' - measurements(:,1:3))./mags', ones(length(mags),1)];
        
       
        L = (P*H')/(H*P*H' + R); % gain calculation
        x_est = x_est + L*(z_meas - z_est); % updating new estimate
        P = (eye(size(P)) - L*H)*P;
        P = real(.5*P + .5*P'); % keeping it real and symmetric

       estimates(i,:) = x_est';
       vars(i,:) = diag(P);
       distance = norm(z_meas - z_est);
       i = i+1;
    end
end
