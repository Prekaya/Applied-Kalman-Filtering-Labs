clc;clearvars
load("sensor_data_w_noise_19.mat")
format long

times = T;
z = z_n;

l = length(times);

H = [ones(l,1), times,times.^2];% Measurement matrix

sigma_MeasNoise = 16; % Measurement noise of radar
R = sigma_MeasNoise*eye(l); 

vars = eye(size(H,2))/(H'/R*H);
estimates = vars*H'/R*z;