% Name         :  Delaunay triangulation skeleton code
% Functionality:  - Create the data structure DelaunayTriangulationSearch 
%                 class to store information gathered during Delaunay 
%                 triangulation optimization;
% 
%                 - Compare the search function values at the inscribe 
%                 center of each Delaunay simplex;
%
%                 - Evaluate the function values at the inscribe center 
%                 which has the smallest values of the search model
%                 (cosntant K version for simple illustration).
%
%                 - Plot the result for each iteration. All figures will be
%                 stored in the 'figures' folder under the current
%                 directory.
%
% Example      :  To use the following data structure, try the following 
%                 1D schwefel example:
%                 ========================================================
%                 n = 2;
%                 xE = [0.2500    ,0.7500    ,0.5000         ,0    ,1.0000    ,1.0000         ,0;
%                       0.2500   , 0.5000  ,  0.87500    ,1.0000       ,  0    ,1.0000     ,    0];
%                 bounds = [zeros(n, 1), ones(n, 1)];
%                 surrogate_ = 'constant';
%                 func_eval = @(x) -sum((2*x) .* sin(sqrt(500*abs(x))));
%                 xmin = 0.8419 * ones(n, 1);
%                 y0 = -1.6759 * n;
%                 K = 3;
%                 obj_lim = [-50, 50];
% 
%                 ddogs = DelaunayTriangulationSearch;
%                 ddogs.initial(n, bounds, surrogate_, func_eval, xE, xmin, y0, K, obj_lim);
%                 ddogs.DeltaDogsOptimize;
%                 ========================================================
% Author       :  Thomas R. Bewley & Muhan Zhao
% Institute    :  Mechanical and Aerospace Engineering, UC San Diego
% Date         :  May. 16, 2019

classdef DelaunayTriangulationSearch < handle
    % This is skeleton code for Delaunay Triangulation search
    % Notice that we have '< handle' to change private properties of DTS at
    % each sub-function call. Otherwise this is a value class and do not
    % change the inherit properties.
    properties
        % sites for constructing Delaunay triangulation
        xE  % Evaluated points.
        yE  % Function values associated with 'xE'.
        % There is no support points because we assume that all corners of 
        % the boundary domain are evaluated.
        % xU
        
        % Objective function
        func_eval  % The function evaluator
        
        % Physical upper and lower bounds
        physical_ub
        physical_lb
        
        % Mesh size info
        mesh_size            % The current mesh size
        num_mesh_refinement  % The maximum times of mesh refinement
        
        % Iteration information
        iter_type  % The type of iteration during the optimization process
        iter       % The number of iteration
        
        % Upper and lower bound of the normalized box domain
        ub  % all to be 1
        lb  % all to be 0
        
        % Delaunay triangulation info
        DT   % The data structure that stores the Delaunay triangulation
        tri  % The index of Delauany triangulation
        DT_center_type  % Declare the type of DT center, 'inscribecenter' or 'circumcenter'
        DT_centers      % The desired centers for each Delaunay triangular
        
        % General information
        n  % The size of the input parameter space
        % The type of the surrogate model, 'constant' or 'adaptive'
        surrogate_type
        
        inter_par  % stores the parameters for interpolation
        xmin       % The global minimizer of the synthetic test problem
        y0         % The target value
        K          % The tradeoff parameter for constant K search model
        
        % Stores the associated function values at the center of Delaunay 
        % Triangulation
        search_values       % values of the search model
        uncertainty_values  % values of the uncertainty function
        interpolation_values     % values of the interpolation function
        
        % Iterative sampling at each iteration
        xc  % minimizer of search model
        yc  % search model value of xc
        xc_eval  % point to evaluate (quantized xc)
        yc_eval  % search model value of xc_eval
        refine_trigger  % trigger of mesh refinement
        
        % plot parameters: 
        % obj_lim: is the objective function values range
        obj_lim
    end
    
    methods
        function initial(self, n, bounds, surrogate_type, func_eval, xE, xmin, y0, K, num_mesh_size, ms)
            % This is the initialization of Delaunay triangulation
            % optimization, the input information is declared by users 
            % given as follows:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % n                      :  The size of input parameter space 
            % bounds                 :  The physical upper and lower bounds
            % surrogate_type         :  The surrogate model type, default to be ''
            % func_eval              :  The function evaluator
            % xE                     :  The initial points to start,
            %                           including corners of the box domain
            % xmin                   :  The global minimizer
            % y0                     :  The target value
            % K                      :  The tradeoff parameter of constantK
            % obj_lim(optional)      :  The objective function limit for
            %                           illustration
            % num_mesh_size(optional):  The maximum time of mesh refinement
            % ms(optional)           :  The initial mesh size
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin < 11
                ms = 8;  % If no mesh size inputed, initialize it with 8
            end
            if nargin < 10
                % If no number of mesh refinement inputed, initialize it with 8
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
            
            % Evaluate the initial points
            for i = 1 : size(self.xE, 2)
                self.yE(i) = func_eval(self.xE(:, i));
            end
            % The following line controls which center to be evaluated,
            % inscribe center = 'inscribecenter', or
            % circumcenter = 'circumcenter'.
            self.DT_center_type = 'inscribecenter';
            
            % create the folder that stores all the figures
            if ~exist('figures', 'dir')
                mkdir('figures')
            end
        end
        
        function DeltaDogsOptimize(self)
            % This is the main function called to optimize 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for kk = 1 : self.num_mesh_refinement
                for k = 1 : 100
                    self.DelaunaySearch
                    if self.refine_trigger == 1
                        break
                    end 
                end
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
            % Each row denotes a center.
            if self.DT_center_type == 'inscribecenter'
                % =====  Using inscribe centers  =====
                self.DT_centers = self.DT.incenter;
            else
                % =====  Using circumcenters  =====
                % Project the circumcenters back into the Delaunay triangulation 
                % if necesssary.
                self.DT_centers = self.DT.circumcenter_projection(self);
            end
            % Build up the function values at each circumcenters.
            func_values_calculator(self);
            
            % Find the next data to evaluate
            [self.yc, func_ind_min] = min(self.search_values);
            self.xc = self.DT_centers(func_ind_min, :)';
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
            self.obj_lim = [min(self.search_values), max(self.search_values)];
            N = size(self.DT_centers, 1);
            X = zeros(N, 1);
            Y = zeros(N, 1);
            Z = zeros(N, 1);
            for i = 1 : N
                point = self.DT_centers(i, :)';
                X(i) = point(1);
                Y(i) = point(2);
                Z(i) = self.search_values(i);
            end
            figure; 
            
            % scatter plot the center of each Delaunay triangulation
            scatter3(X, Y, Z, 30, 'filled')
            
            zlim(self.obj_lim);
            hold on; grid on;
            [~, ind_min] = min(Z);
            
            % scatter plot the minimizer of all centers
            scatter3(X(ind_min), Y(ind_min), Z(ind_min), 'r', 'filled')
            scatter3(X, Y, self.obj_lim(1) * ones(N, 1), 15, 'filled')
            
            % plot the vertical line of each center, connecting from the
            % search model function values to the obj_lim lower bound
            for i = 1 : N
                line([X(i); X(i)], [Y(i); Y(i)], [Z(i), self.obj_lim(1)], 'Color', 'black', 'Linestyle', '--');
            end
            
            % plot the global minimum
            scatter3(self.xmin(1), self.xmin(2), self.obj_lim(1), 'r*')
            
            % plot the quantized point of minimizer of all centers
            scatter3(self.xc_eval(1), self.xc_eval(2), self.obj_lim(1), 'g', 'filled')
            view(45, 45)
            DelaunayTriangulationPlot2D(self);
            if self.mesh_size > 16
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
            end
            set(gca,'ZTickLabel',[])
            set(gca, 'xtick', [0:(1/self.mesh_size):1])
            set(gca, 'ytick', [0:(1/self.mesh_size):1])
            saveas(gca, ['figures/2D_DT_', num2str(self.iter), '.png'])
        end
        function DelaunayTriangulationPlot2D(self)
            % plot the boundary lines of each Delaunay triangulation
            N = size(self.DT_centers, 1);
            comb = nchoosek(1:3, 2);
            for i = 1 : N
                simplex = self.xE(:, self.tri(i, :));
                for j = 1 : size(comb, 1)
                    line(simplex(1, comb(j, :)), simplex(2, comb(j, :)), [self.obj_lim(1), self.obj_lim(1)], 'Color', 'green');
                end
            end
        end
        function DelaunayUpdate(self)
            % update the evaluated point at the end of each iteration
            if self.refine_trigger ~= 1
                self.xE = [self.xE, self.xc_eval];
                self.yE = [self.yE, self.yc_eval];
            else
                self.mesh_size = self.mesh_size * 2;
                % if using constant K surrogate
                self.K = self.K * 2;
            end
        end
    end
end




