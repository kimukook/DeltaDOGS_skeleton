n = 2;
xE = [0.2500    ,0.7500    ,0.5000         ,0    ,1.0000    ,1.0000         ,0;
    0.2500   , 0.5000  ,  0.8500    ,1.0000       ,  0    ,1.0000     ,    0];
bounds = [zeros(n, 1), ones(n, 1)];
surrogate_ = 'constant';
func_eval = @(x) -sum((2*x) .* sin(sqrt(500*abs(x))));
xmin = 0.8419 * ones(n, 1);
y0 = -1.6759 * n;
K = 3;
obj_lim = [-50, 50];

ddogs = DelaunayTriangulationSearch;
ddogs.initial(n, bounds, surrogate_, func_eval, xE, xmin, y0, K, obj_lim);
ddogs.DelaunaySearch;


% check the circumcenter project for this simplex
n = 2;
simplex =[0.7500,    0.7500 ,   1.0000;
    0.9250 ,   0.5000   , 1.0000];
center =  [   0.9500;
    0.7125];
simplex_center =  simplex * ones(n + 1, 1) / (n + 1);

figure;
plot_simplex = [simplex, simplex(:, 1)];
plot(plot_simplex(1,:), plot_simplex(2,:));
hold on ; grid on
scatter(center(1), center(2), 'r', 'filled')
scatter(simplex_center(1), simplex_center(2), 'b', 'filled')
line([center(1), simplex_center(1)],[center(2), simplex_center(2)])
scatter(points(1, 1), points(2, 1), 'g', 'filled')