%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef stem_EM_options
    properties
        exit_toll=0.0001;            %[double >0] (1x1) the EM algorithm stops if the relative norm between two consecutive iterations is below exit_toll
        max_iterations=100;          %[integer >0](1x1) the EM algorithm stops if the number of iterations exceed max_iterations
        numeric_opt_type='single';   %[string]    (1x1) 'single': then elements of the V_i matrices are numerically estimated one-by-one; 'full': the elements are jointly estimated.
        mstep_system_size=3500;      %[integer >0](1x1) if N_r(N_g)>mstep_system_size then theta_r and v_r (theta_g and v_g) are optimized by considering diagonal blocks of maximum dimension mstep_system_size
        compute_logL_at_all_steps=0; %[boolean]   (1x1) 1: the observed data log-likelihood is evaluated at each iteration of the EM algorithm
        verbose=1;                   %[boolean]   (1x1) 1: all the intermediate operations of the EM algorithm are displayed
        pathparallel=[];             %[string]    (1x1) full or relative path of the folder to use for parallel computation
    end
    
    methods
        function obj = stem_EM_options(exit_toll,max_iterations,numeric_opt_type,mstep_system_size,compute_logL_at_all_steps,verbose,pathparallel)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %<exit_toll>                 - [double >0]         (1x1) (default: 0.0001) the EM algorithm stops if the relative norm between two consecutive iterations is below exit_toll
            %<max_iterations>            - [integer >0]        (1x1) (default: 1000)  the EM algorithm stops if the number of iterations exceed max_iterations
            %<numeric_opt_type>          - [string]            (1x1) (default: 'single') 'single': then elements of the V_i matrices are numerically estimated one-by-one; 'full': the elements are jointly estimated.
            %<mstep_system_size>         - [integer >0]        (1x1) (default: 3500) if N_r(N_g)>mstep_system_size then theta_r and v_r (theta_g and v_g) are optimized by considering diagonal blocks of maximum dimension mstep_system_size
            %<compute_logL_at_all_steps> - [boolean]           (1x1) (dafault: 0) 1: the observed data log-likelihood is evaluated at each iteration of the EM algorithm
            %<verbose>                   - [boolean]           (1x1) (default: 0) 1: all the intermediate operations of the EM algorithm are displayed
            %<pathparallel>              - [string]            (1x1) (default: []) full or relative path of the folder to use for parallel computation
            %
            %
            %OUTPUT
            %obj       - [stem_EM_options object] (1x1)
            if nargin>0
                obj.exit_toll=exit_toll;
            end
            if nargin>1
                obj.max_iterations=max_iterations;
            end
            if nargin>2
                obj.numeric_opt_type=numeric_opt_type;
            end
            if nargin>3
                obj.mstep_system_size=mstep_system_size;
            end
            if nargin>4
                obj.compute_logL_at_all_steps=compute_logL_at_all_steps;
            end
            if nargin>5
                obj.verbose=verbose;
            end
            if nargin>6
                obj.pathparallel=pathparallel;
            end
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
        
        function obj = set.numeric_opt_type(obj,numeric_opt_type)
            if not(isempty(numeric_opt_type))
                if not(strcmp(numeric_opt_type,'single'))&&not(strcmp(numeric_opt_type,'full'))
                    error('numeric_opt_type must be either ''single'' or ''full''');
                end
                obj.numeric_opt_type=numeric_opt_type;
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
        
        function obj = set.pathparallel(obj,pathparallel)
            if not(isempty(pathparallel))
                if not(ischar(pathparallel))
                    error('pathparallel must be a string');
                else
                    obj.pathparallel=pathparallel;
                end
            end
        end
    end
end