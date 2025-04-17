function ellipse_plot = plot_n_sigma_ellipse(cov_mat, n, center)
    arguments
        cov_mat;
        n;
        center = [0; 0];
    end
    
    phi = linspace(0,2*pi,100);
    [eig_vec, eig_value] = eig(cov_mat, 'vector');

    lambda = sort(n*sqrt(eig_value));
    xy = lambda.*[cos(phi);sin(phi)];

    theta = -atan2(eig_vec(1,2), eig_vec(2,2));
    Sxy = center + [cos(theta) -sin(theta);sin(theta) cos(theta)]*xy;
    ellipse_plot = plot(Sxy(1,:),Sxy(2,:));
end