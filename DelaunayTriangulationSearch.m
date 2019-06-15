% Name         :  Delaunay triangulation skeleton code
% Functionality:  Create the Delaunay triangulation class;
%              :  Compare the search function values at the center of 
%                 each Delaunay simplex;
%              :  Evaluate the circumcenter with the lowest search values
%                 for next function sampling.
%              :  Plot the result for each iteration.
%
% Author       :  Muhan Zhao
% Institute    :  Mechanical and Aerospace Engineering, UC San Diego
% Date         :  May. 16, 2019

classdef DelaunayTriangulationSearch < handle
    % This is skeleton code for Delaunay Triangulation search
    % Notice that we have '< handle' to change private properties of DTS at
    % each sub-function call. Otherwise this is a value class and do not
    % change the inherit properties.
    properties
        % sites for constructing Delaunay triangulation
        xE
        yE
        % There is no support points because all corners of boundary are
        % evaluated.
        % xU
        
        % objective function
        func_eval
        % physical bounds
        physical_ub
        physical_lb
        % Mesh size info
        mesh_size
        num_mesh_refinement
        % Iteration information
        iter_type
        iter
        % Upper and lower bound
        ub
        lb
        % Delaunay triangulation info
        DT
        tri
        DT_circumcenters
        incenter
        % General information
        n
        surrogate_type
        inter_par
        xmin
        y0
        K
        % Delaunay Triangulation functions' info
        search_func_values
        uncertainty_func_values
        interpolation_values
        % Iterative sampling at each iteration
        xc
        yc
        xc_eval
        yc_eval
        refine_trigger
        % plot parameters: 
        % obj_lim: is the objective function values range
        obj_lim
    end
    methods
        function initial(self, n, bounds, surrogate_type, func_eval, xE, xmin, y0, K, obj_lim, num_mesh_size, ms)
            if nargin < 12
                ms = 8;
            end
            if nargin < 11
                num_mesh_size = 8;
            end
            self.iter = 0;
            self.n = n;
            self.func_eval = func_eval;
            self.xE = xE;
            self.physical_lb = bounds(:, 1);
            self.physical_ub = bounds(:, 2);
            self.lb = zeros(n, 1);
            self.ub = ones(n, 1);
            self.yE = zeros(1, size(self.xE, 2));
            self.obj_lim = obj_lim;
            self.surrogate_type = surrogate_type;
            self.num_mesh_refinement = num_mesh_size;
            self.mesh_size = ms;
            if self.surrogate_type == 'constant'
                self.K = K;
            else
                self.y0 = y0;
            end
            self.xmin = xmin;
            self.yE = zeros(1, size(self.xE, 2));
            for i = 1 : size(self.xE, 2)
                self.yE(i) = func_eval(self.xE(:, i));
            end
        end
        function DeltaDogsOptimize(self)
            for k = 1 : 30
                self.DelaunaySearch
            end
        end
        function DelaunaySearch(self)
            % This is the Delaunay simple search at each iteration:
            % Evaluate the p(x), e(x), s(x) at each circumcenter (or, its
            % projection), of each Delaunay simplices.
            self.iter = self.iter + 1;
            % Create the interpolation
            inter = NPSInterpolation;
            self.inter_par = inter.interpolateparametarization(self.xE, self.yE);
            % Construct the Delaunay triangulation and all functions'
            % values.
            self = DelaunayConstruction(self);
            % Generate the plot
            plot_maker2D(self);
            % Update the evaluated data and proceeding
            DelaunayUpdate(self);
        end
        function self = DelaunayConstruction(self)
            % Construct the Delaunay triangulation
            self.DT = delaunayTriangulation(self.xE.');
            self.tri = self.DT.ConnectivityList;
            self.DT_circumcenters = self.DT.circumcenter;
            % Project the circumcenters if necesssary.
            self.DT_circumcenters = circumcenter_projection(self);
            % Build up the function values at each circumcenters.
            [self.search_func_values, self.uncertainty_func_values, ...
                self.interpolation_values] = search_uncertainty_func(self);
            % Find the next data to evaluate
            [self.yc, func_ind_min] = min(self.search_func_values);
            % Use circumcenter of DT
            self.xc = self.DT_circumcenters(func_ind_min, :)';
            % Use incenter of DT
%             self.xc = self.incenter(func_ind_min, :)';
            % TODO:
            % quantize the minimizer of search function values inside 
            % the incenter list.
            self.xc_eval = round(self.xc * self.mesh_size) / self.mesh_size;
            self.yc_eval = self.func_eval(self.xc_eval);
            dis = norm2_dis(self.xc_eval, self.xE);
            if min(dis) < 1e-6
                self.refine_trigger = 1;
            else
                self.refine_trigger = 0;
            end 
        end
        function plot_maker2D(self)
            N = size(self.DT_circumcenters, 1);
            X = zeros(N, 1);
            Y = zeros(N, 1);
            Z = zeros(N, 1);
            for i = 1 : N
                % Use circumcenter of DT
                point = self.DT_circumcenters(i, :)';
                % Use incenter of DT
%                 point = self.incenter(i, :)';
                X(i) = point(1);
                Y(i) = point(2);
                Z(i) = self.search_func_values(i);
            end
            figure; 
            scatter3(X, Y, Z, 30, 'filled')
            
            zlim(self.obj_lim);
            hold on; grid on;
            [~, ind_min] = min(Z);
            scatter3(X(ind_min), Y(ind_min), Z(ind_min), 'r', 'filled')
            scatter3(X, Y, self.obj_lim(1) * ones(N, 1), 15, 'filled')
            for i = 1 : N
                line([X(i); X(i)], [Y(i); Y(i)], [Z(i), self.obj_lim(1)], 'Color', 'black', 'Linestyle', '--');
            end
            scatter3(self.xmin(1), self.xmin(2), self.obj_lim(1), 'r*')
            scatter3(self.xc_eval(1), self.xc_eval(2), self.obj_lim(1), 'g', 'filled')
            view(45, 45)
            DelaunayTriangulationPlot2D(self);
            if self.K > 20
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
            end
            set(gca,'ZTickLabel',[])
            set(gca, 'xtick', [0:(1/self.mesh_size):1])
            set(gca, 'ytick', [0:(1/self.mesh_size):1])
            saveas(gca, ['figures/2D_DT_', num2str(self.iter), '.png'])
        end
        function DelaunayTriangulationPlot2D(self)
            N = size(self.DT_circumcenters, 1);
            comb = nchoosek(1:3, 2);
            for i = 1 : N
                simplex = self.xE(:, self.tri(i, :));
                for j = 1 : size(comb, 1)
                    line(simplex(1, comb(j, :)), simplex(2, comb(j, :)), [self.obj_lim(1), self.obj_lim(1)], 'Color', 'green');
                end
            end
        end
        function DelaunayUpdate(self)
            if self.refine_trigger ~= 1
                self.xE = [self.xE, self.xc_eval];
                self.yE = [self.yE, self.yc_eval];
            else
                self.mesh_size = self.mesh_size * 2;
                self.K = self.K * 2;
            end
        end
    end
end




