function [search_values, uncertainty_values, interpolation_values] = search_uncertainty_func(self)
% Calculate the search function values and uncertainty function values at
% the circumcenters (or, its projection) for each Delaunay simplices.
% 
% Author:   Muhan Zhao
% Date  :   May. 21, 2019
N = size(self.DT_circumcenters, 1);
search_values = zeros(N, 1);
uncertainty_values = zeros(N, 1);
interpolation_values = zeros(N, 1);
self.incenter = self.DT.incenter;
for i = 1 : N
    simplex = self.xE(:, self.tri(i, :));
    % Use circumcenter of DT
    center = self.DT_circumcenters(i, :)';
    % % Use incenter of DT
%     center = self.incenter(i, :)';
    [xc, R2] = circhyp(simplex, self.n);
    interpolation_values(i) = self.inter_par.interpolate_eval(center);
    uncertainty_values(i) = R2 - norm(center - xc)^2;
    search_values(i) = interpolation_values(i) - self.K * uncertainty_values(i);
end
end