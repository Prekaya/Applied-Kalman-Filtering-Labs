s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);
set1_P = [6 2;
          2 2];

set2_P = [6 0;
          0 2];

figure()
scatter_sigmas(set1_P);
grid on

figure()
scatter_sigmas(set2_P);
grid on

function scatter_sigmas(P)
    L = chol(P);
    xx = L'*randn(2,5000);
    r = cov(xx');
    hold on
    scatter(xx(1,:),xx(2,:),2,'filled');
    sigma_1 = plot_n_sigma_ellipse(r,1);
    sigma_1.Color = "red";
    
    sigma_2 = plot_n_sigma_ellipse(r,2);
    sigma_2.Color = "black";
    
    sigma_3 = plot_n_sigma_ellipse(r,3);
    sigma_3.Color = "green";
    title(mat2str(r,4))
    hold off
end