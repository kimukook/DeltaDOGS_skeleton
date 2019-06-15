function [center_points] = circumcenter_projection(self)
center_points = self.DT_circumcenters;
N = size(center_points, 1);
for i = 1 : N
    center = center_points(i, :)';
    % First check if the circumcenter is within the DT simplex or not.
    simplex = self.xE(:, self.tri(i, :));
    [A, b] = search_simplex_bounds(simplex);
    if min((A * center - b) >= 0) == 0
        % all simplex constraints satisfied
        center_points(i, :) = outside_circumcenter_projections(center, simplex, A, b);
    else
        continue
    end
end
end

function [point] = outside_circumcenter_projections(center, simplex, A, b)
n = size(simplex, 1);
N = size(A, 1);
simplex_center = simplex * ones(n + 1, 1) / (n + 1);
coef = simplex_center - center;
syms t;
line_vector = coef * t + center;
% dis: distance from circumcenter to the projected points
dis = zeros(N, 1);
% feasb: feasibility of the projected points, if they are outside of the
% boundary then omit.
feasb = zeros(N, 1);
% feasible checker: Ain and bin: Ain * x >= bin
Ain = [-eye(n); eye(n)]; bin = [-ones(n, 1); zeros(n, 1)];
% points: projected points
points = zeros(n, N);
for i = 1 : N
    % First project each point onto the hyperplane, 
    sol = solve(A(i, :) * line_vector - b(i), t);
    points(:, i) = (coef * sol + center)';
    % check the feasibility
    if min(Ain * points(:, i) - bin) >= 0
        feasb(i) = 1;
    else
        feasb(i) = 0;
    end
    dis(i) = norm(center - points(:, i));
end
[~, ind] = min(dis(feasb==1));
point = points(:, ind);
end