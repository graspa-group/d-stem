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
        exit_toll=0.0001;            %relative exit tolerance on both parameter norm and log-likelihood
        max_iterations=100;          %maximum number of EM iterations
        numeric_opt_type='single';
        mstep_system_size=2500;
        compute_logL_at_all_steps=0;
        verbose=1;
        pathparallel=[];
    end
    
    methods
        function obj = stem_EM_options(exit_toll,max_iterations,numeric_opt_type,mstep_system_size,compute_logL_at_all_steps,verbose,pathparallel)
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