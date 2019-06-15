function [A, b] = search_simplex_bounds(xi)
% For n dimensional simplex xi, return the n+1 simplex bounds A and b
% such that Ax >= b.
% Input :        xi -> The dimension should be n by (n+1), dim-by-num-point.
% Output:        A and b such that Ax >= b is within the simplex.
%
% Author:        Muhan Zhao
% Institute    :  Mechanical and Aerospace Engineering, UC San Diego
% Data  :        May. 17, 2019
n = size(xi, 1); N = size(xi, 2);
A = zeros(N, n); b = zeros(N, 1);
for i = 1 : N
    direction_pointer = xi(:, i);
    extract_points = 1 : N;
    extract_points(i) = [];
    plane_points = xi(:, extract_points);
    [A(i, :), b(i)] = hyperplane_equations(plane_points, n);
    if A(i, :) * direction_pointer < b(i)
        A(i, :) = -A(i, :);
        b(i) = -b(i);
    end
end
end