clear, clc, close all
n = 2;
xE = [0.2500    ,0.7500    ,0.5000         ,0    ,1.0000    ,1.0000         ,0;
      0.2500   , 0.5000  ,  0.87500    ,1.0000       ,  0    ,1.0000     ,    0];
bounds = [zeros(n, 1), ones(n, 1)];
surrogate_ = 'adaptive'; % or use 'constant' to switch to constant surrogate model
func_eval = @(x) -sum((2*x) .* sin(sqrt(500*abs(x))));
xmin = 0.8419 * ones(n, 1);
y0 = -1.6759 * n;
K = 3;
num_mesh_size = 5;

ddogs = DelaunayTriangulationSearch;
ddogs.initial(n, bounds, surrogate_, func_eval, xE, xmin, y0, K, num_mesh_size);
ddogs.DeltaDogsOptimize;