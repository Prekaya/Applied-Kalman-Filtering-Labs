s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

P = [6 2;
     2 2];

[D, Q, V] = svd(P);

w = D*sqrt(Q);
xx = w*randn(2,5000);
r = cov(xx');
scatter(xx(1,:),xx(2,:),2,'filled');
title(mat2str(r,4))