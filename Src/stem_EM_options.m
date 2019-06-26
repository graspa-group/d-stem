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
        exit_toll=0.0001;                           %[double >0]  (1x1) the EM algorithm stops if the relative norm between two consecutive iterations is below exit_toll
        max_iterations=100;                         %[integer >0] (1x1) the EM algorithm stops if the number of iterations exceed max_iterations
        fminsearch_max_iterations=100;              %[integer >0] (1x1) the fminsearch method used for updating model parameters stops if the number of iterations exceed fminsearch_max_iterations
       
        workers=1;                                  %[integer>0]  (1x1) the number of matlab workers used for EM algorithm

        block_tapering_block_size=0;                %[integer >=0](1x1)|(Bx1) the dimension of blocks in block tapering or the Bx1 vector of block dimensions. If equal to 0, block tapering is NOT enabled
        compute_logL_at_all_steps=1;                %[boolean]    (1x1) 1: the observed data log-likelihood is evaluated at each iteration of the EM algorithm
        
        verbose=1;                                  %[boolean]    (1x1) 1: all the intermediate operations of the EM algorithm are displayed

        wpar_estimation=1;                          %[boolean]    (1x1) 0: parameters \theta and V are not estimated; 1: parameters are estimated
        mstep_system_size=1500;                     %[integer >0] (1x1) if N_r(N_g)>mstep_system_size then theta_r and v_r (theta_g and v_g) are optimized by considering diagonal blocks of maximum dimension mstep_system_size
    end
    
    methods
        function obj = stem_EM_options(exit_toll,max_iterations)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %exit_toll          - [double >0]              (1x1) the EM algorithm stops if: 1) the relative norm of the model parameter vector between two iterations is below exit toll; 2) the relative norm of the observed data log-likelihood between two iterations is below exit toll (if computed)
            %max_iterations     - [integer >0]             (1x1) the EM algorithm stops if the number of iterations exceed max_iterations
            %
            %OUTPUT
            %obj                - [stem_EM_options object] (1x1)
            
            if nargin<2
                error('exit_toll and max_iterations must be provided when creating objects of class stem_EM_options');
            end
            
            obj.exit_toll=exit_toll;
            obj.max_iterations=max_iterations;
        end
        
        %Class set methods
        function obj = set.exit_toll(obj,exit_toll)
            if not(isempty(exit_toll))
                if exit_toll<=0
                    error('The exit_toll must be >0');
                end
                obj.exit_toll=exit_toll;
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
                if fmins  earch_max_iterations<=0
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
       
        function obj = set.block_tapering_block_size(obj,block_tapering_block_size)
            if length(block_tapering_block_size)==1
                if block_tapering_block_size<0
                    error('block_tapering_block_size must be >=0');
                end
            else
                if sum(block_tapering_block_size<=0)>1
                    error('The elements of block_tapering_block_size must be > 0');
                end
            end
            obj.block_tapering_block_size=block_tapering_block_size;
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