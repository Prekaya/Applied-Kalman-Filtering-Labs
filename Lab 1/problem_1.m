% R = RANDN(N) returns an N-by-N matrix containing pseudo-random
% values drawn from a normal distribution with mean zero and standard
% deviation of 1
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);
R_0_1 = Normal(0,1,100000);
disp([mean(R_0_1), cov(R_0_1)]);

R_0_16 = Normal(0,16,100000);
disp([mean(R_0_16), sqrt(var(R_0_16))]);

R_5_25 = Normal(5,25,100000);
disp([mean(R_5_25), sqrt(var(R_5_25))]);

function R = Normal(mean,covariance,n)
    R = mean + covariance*randn(n,1);
end