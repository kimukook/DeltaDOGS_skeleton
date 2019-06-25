function func_values_calculator(self)
% Calculate the search function values, uncertainty function values and the
% interpolation function values at the centers (circumcenter or its projection, 
% or inscribe center) for each Delaunay simplices.
% 
% Author:   Muhan Zhao
% Date  :   May. 21, 2019
N = size(self.DT_centers, 1);
self.interpolation_values = zeros(N, 1);
self.uncertainty_values   = zeros(N, 1);
self.search_values        = zeros(N, 1);
for i = 1 : N
    simplex  = self.xE(:, self.tri(i, :));
    center   = self.DT_centers(i, :)';
    [xc, R2] = circhyp(simplex, self.n);
    % calculate the function values
    self.interpolation_values(i) = self.inter_par.interpolate_eval(center);
    self.uncertainty_values(i)   = R2 - norm(center - xc)^2;
    if self.surrogate_type == 'constant'
        self.search_values(i)    = self.interpolation_values(i) - ...
                                       self.K * self.uncertainty_values(i);
    else
        self.search_values(i)    = (self.interpolation_values(i) - self.y0)/...
                                        self.uncertainty_values(i);
    end
end
end