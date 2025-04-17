s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);
P = [6 2;
          2 2];
L = chol(P);
xx = L'*randn(2,5000);
r = cov(xx');

figure()
scatter(xx(1,:),xx(2,:),2,'filled');
title(mat2str(r,4))
grid on

white_xx = L\xx;
white_r = cov(white_xx');

figure()
scatter(white_xx(1,:),white_xx(2,:),2,'filled')
title(mat2str(white_r,4))
grid on