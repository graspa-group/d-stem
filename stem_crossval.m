%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef stem_crossval < handle
    
    properties
        variable_name=[];
        idx_step=[];
        type=[];
        
        stem_varset=[];         %stem_varset object for crossvalidation
        stem_gridlist=[];       %gridlist object for crossvalidation
        stem_krig_result=[];
        min_distance=[];         %for each cross-validation site contains the distance to the nearest non-cross-validation site (of the same variable)
        mse=[];
        relative_mse=[];
        avg_relative_mse=[];
        avg_relative_mse_higher50km=[];
        res=[];
        res_backtransformed=[];
        relative_res=[];
        y_hat_back=[];
        y_back=[];
    end

    
    methods
        function obj = stem_crossval(type,variable_name,idx_step)
            if nargin<3
                error('Not enough input arguments');
            end
            obj.type=type;
            obj.variable_name=variable_name;
            obj.idx_step=idx_step;
        end
       
        function set.type(obj,type)
            if strcmp(type,'ground')==0
                error('The cross-validation is available only for ground variables. The type must be equal to ground.');
            end
            obj.type=type;
        end
        
        function set.idx_step(obj,idx_step)
            if idx_step<2
                error('idx_step must be >=2');
            end
            obj.idx_step=idx_step;
        end
        
        function set.stem_varset(obj,stem_varset)
           if not(isa(stem_varset,'stem_varset'))
               error('stem_varset must be of class stem_varset');
           end
           obj.stem_varset=stem_varset;
        end
        
        function set.stem_gridlist(obj,stem_gridlist)
            if not(isa(stem_gridlist,'stem_gridlist'))
                error('stem_gridlist must be of class stem_gridlist');
            end
            
            for i=1:length(stem_gridlist.grid)
                if not(size(stem_gridlist.grid{i}.coordinate,1)==size(obj.stem_varset.Y{i},1))
                    error('The number of coordinates in the grid{i} must be equal to the number of rows of Y{i}');
                end
                if not(strcmp(stem_gridlist.grid{i}.site_type,'point'))
                    error('Only point data are supported in stem_gridlist');
                end
            end            
            
            obj.stem_gridlist=stem_gridlist;
        end            
    end
    
end