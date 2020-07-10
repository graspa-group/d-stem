%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% D-STEM - Distributed Space Time Expecation Maximization              %
%%%                                                                      %
%%% Author: Francesco Finazzi                                            %
%%% E-mail: francesco.finazzi@unibg.it                                   %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Information and                    %
%%%              Production Engineering                                  %
%%% Author website: http://www.unibg.it/pers/?francesco.finazzi          %
%%%                                                                      %
%%% Author: Yaqiong Wang                                                 %
%%% E-mail: yaqiongwang@pku.edu.cn                                       %
%%% Affiliation: Peking University,                                      %
%%%              Guanghua school of management,                          %
%%%              Business Statistics and Econometrics                    %
%%%                                                                      %
%%% Author: Alessandro Fassò                                             %
%%% E-mail: alessandro.fasso@unibg.it                                    %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Information and                    %
%%%              Production Engineering                                  %
%%% Author website: http://www.unibg.it/pers/?alessandro.fasso           %
%%%                                                                      %
%%% Code website: https://github.com/graspa-group/d-stem                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% This file is part of D-STEM.
% 
% D-STEM is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 2 of the License, or
% (at your option) any later version.
% 
% D-STEM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with D-STEM. If not, see <http://www.gnu.org/licenses/>.

classdef stem_EM_options
    
    %PROPERTIES
    %Each class property or method property is defined as follows
    %
    %"Name"="Default value";    %["type"]    "dimension"     "description" 
    %
    %DIMENSION NOTATION
    %(1 x 1) is a scalar
    %(N x 1) is a Nx1 vector
    %(N x T) is a NxT matrix
    %(N x B x T) is a NxBxT array
    %{q} is a cell array of length q
    %{q}{p} is a cell array of length q, each cell is a cell array of length p
    %{q}(NxT) is a cell array of length q, each cell is a NxT matrix
    %
    %CONSTANTS
    %B - dimension of vector of block dimensions
    
    properties
        exit_tol_par=0.0001;                        %[double >0]  (1x1) the EM algorithm stops if the maximum relative norm of the model parameters between two consecutive iterations is below exit_tol_par
        exit_tol_loglike=0.0001;                    %[double >0]  (1x1) the EM algorithm stops if the relative norm of the log-likelihood between two consecutive iterations is below exit_tol_loglike

        max_iterations=100;                         %[integer >0] (1x1) the EM algorithm stops if the number of iterations exceed max_iterations
        fminsearch_max_iterations=100;              %[integer >0] (1x1) the fminsearch method used for updating model parameters stops if the number of iterations exceed fminsearch_max_iterations
       
        workers=1;                                  %[integer>0]  (1x1) the number of matlab workers used for EM algorithm

        partitions=0;                               %[integer >=0](1x1)|(Bx1) the dimension each partition or the Bx1 vector partition sizes. If equal to 0, partitioning is NOT enabled
        compute_logL_at_all_steps=1;                %[boolean]    (1x1) 1: the observed data log-likelihood is evaluated at each iteration of the EM algorithm
        
        verbose=1;                                  %[boolean]    (1x1) 1: all the intermediate operations of the EM algorithm are displayed

        wpar_estimation=1;                          %[boolean]    (1x1) 0: parameters \theta and V are not estimated; 1: parameters are estimated
        mstep_system_size=1500;                     %[integer >0] (1x1) if N_r(N_g)>mstep_system_size then theta_r and v_r (theta_g and v_g) are optimized by considering diagonal blocks of maximum dimension mstep_system_size
    end
    
    methods
        function obj = stem_EM_options()
            %DESCRIPTION: object constructor
            %
            %OUTPUT
            %obj - [stem_EM_options object] (1x1)
        end
        
        %Class set methods
        function obj = set.exit_tol_par(obj,exit_tol_par)
            if not(isempty(exit_tol_par))
                if exit_tol_par<=0
                    error('exit_tol_par must be >0');
                end
                obj.exit_tol_par=exit_tol_par;
            end
        end
        
        function obj = set.exit_tol_loglike(obj,exit_tol_loglike)
            if not(isempty(exit_tol_loglike))
                if exit_tol_loglike<=0
                    error('exit_tol_loglike must be >0');
                end
                obj.exit_tol_loglike=exit_tol_loglike;
            end
        end
        
        function obj = set.max_iterations(obj,max_iterations)
            if not(isempty(max_iterations))
                if max_iterations<=0
                    error('max_iterations must be >0');
                end
                obj.max_iterations=max_iterations;
            end
        end
        
        function obj = set.fminsearch_max_iterations(obj,fminsearch_max_iterations)
            if not(isempty(fminsearch_max_iterations))
                if fminsearch_max_iterations<=0
                    error('max_iterations for estimating parameter v and theta must be >0');
                end
                obj.fminsearch_max_iterations=fminsearch_max_iterations;
            end
        end
        
        function obj = set.mstep_system_size(obj,mstep_system_size)
            if not(isempty(mstep_system_size))
                if mstep_system_size<=0
                    error('mstep_system_size must be >0');
                end
                obj.mstep_system_size=mstep_system_size;
            end
        end
        
        function obj = set.compute_logL_at_all_steps(obj,compute_logL_at_all_steps)
            if not(isempty(compute_logL_at_all_steps))
                if not(compute_logL_at_all_steps==0)&&not(compute_logL_at_all_steps==1)
                    error('compute_logL_at_all_steps must be either 0 or 1');
                end
                obj.compute_logL_at_all_steps=compute_logL_at_all_steps;
            end
        end
        
        function obj = set.verbose(obj,verbose)
            if not(isempty(verbose))
                if not(verbose==0)||not(verbose==1)
                    error('compute_logL_at_all_steps must be either 0 or 1');
                end
                obj.verbose=verbose;
            end
        end
       
        function obj = set.partitions(obj,partitions)
            if length(partitions)==1
                if partitions<0
                    error('partitions must be >=0');
                end
            else
                if sum(partitions<=0)>1
                    error('The elements of partitions must be > 0');
                end
            end
            obj.partitions=partitions;
        end
        
        function obj = set.workers(obj,workers)
            if not(isempty(workers))
                if workers<=0
                    error('workers must be >0');
                elseif rem(workers,1)~=0
                    error('workers must be integer');
                end
                obj.workers=workers;
            end
        end
  
    end
end