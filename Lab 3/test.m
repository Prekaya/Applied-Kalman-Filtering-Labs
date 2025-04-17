x = [1;2;3];

y = [1,2,3;
    3,2,3;
    1,5,3];
r = vecnorm(x-y',2)
l = [4,2,1];
b = [(x-y')./l', ones(length(r),1)]

estimate(1:4,1) = 3;
