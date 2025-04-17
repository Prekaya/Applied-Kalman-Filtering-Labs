P = [4, 1;
     1, 1];

hold on
sigma_1 = plot_n_sigma_ellipse(P,1);
sigma_1.Color = "blue";

sigma_2 = plot_n_sigma_ellipse(P,2);
sigma_2.Color = "red";

sigma_3 = plot_n_sigma_ellipse(P,3);
sigma_3.Color = "black";
hold off
legend('1\sigma', '2\sigma', '3\sigma')
grid on