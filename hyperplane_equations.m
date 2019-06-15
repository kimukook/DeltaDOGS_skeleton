function [a, b] = hyperplane_equations(points, n)
% points should be the n-by-n matrix, num_dim-by-num_points.
% n is the dimension
basepoint = points(:, end);
matrix = (points - repmat(basepoint, [1, n]));
matrix = matrix(:, 1:end-1)';
a = zeros(1, n);
for j = 1 : n
    extract_points = 1 : n;
    extract_points(j) = [];
    block = matrix(:, extract_points);
    a(1, j) = (-1)^(j) * det(block);
end
b = a * basepoint;
end