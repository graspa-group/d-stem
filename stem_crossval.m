%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% D-STEM - Distributed Space Time Expecation Maximization              %
%%%                                                                      %
%%% Author: Francesco Finazzi                                            %
%%% E-mail: francesco.finazzi@unibg.it                                   %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Economics and Quantitative Methods %
%%% Author website: http://www.unibg.it/pers/?francesco.finazzi          %
%%% Code website: https://code.google.com/p/d-stem/                      %
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



classdef stem_crossval < handle
    
    %CONSTANTS
    %
    %dN - the number of cross-validation sites
    %T  - the number of time steps
    
    properties
        variable_name=[];           %[string]                  (1x1) the name of the cross-validation variable
        type=[];                    %[string]                  (1x1) 'point': cross-validation is on point-data
        indices=[];                 %[integer >0]              (dNx1) the indices of the cross-validation sites 
        
        stem_varset=[];             %[stem_varset object]      (1x1)  the subset of data used for cross-validation
        stem_gridlist=[];           %[strm_gridlist object]    (1x1)  this object is needed for storing the coordinates of the cross-validation sites
        stem_krig_result=[];        %[strm_krig_result object] (1x1)  the prediction over the cross-validation sites is obtained through kriging
        
        y_back=[];                  %[double]                  (dNxT) the original data (back-transformed) of the cross-validation sites
        y_hat_back=[];              %[double]                  (dNxT) the estimated data (back-transformed) for the cross-validation sites
        res=[];                     %[double]                  (dNxT) cross-validation residuals        
        res_back=[];                %[double]                  (dNxT) back-transformed cross-validation residuals (if the original data are log-transformed and/or standardized)
     
        mse=[];                     %[double >=0]              (dNx1) mean squared error for each cross-validation site
        mse_time=[];                %[double >=0]              (dNx1) mean squared error for each time step
        relative_mse=[];            %[double >=0]              (dNx1) relative mean squared error for each cross-validation site
        relative_mse_time=[];       %[double >=0]              (dNx1) relative mean squared error for each time step
        min_distance=[];            %[double >=0]              (dNx1) the distance from each cross-validation site to the nearest non cross-validation site
    end

    
    methods
        function obj = stem_crossval(variable_name,type,indices)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %variable_name - [string]               (1x1) the name of the cross-validation variable
            %type          - [string]               (1x1) 'point': cross-validation is on point-data
            %indices       - [integer >0]           (dNx1) the indices of the cross-validation sites 
            %
            %OUTPUT
            %obj           - [stem_crossval object] (1x1)    
            
            if nargin<3
                error('Not enough input arguments');
            end
            obj.variable_name=variable_name;
            obj.type=type;            
            obj.indices=indices;
        end
       
        %Class set methods
        function set.type(obj,type)
            if not(strcmp(type,'point'))
                error('The cross-validation is available only for point variables. The type must be equal to point.');
            end
            obj.type=type;
        end
        
        function set.indices(obj,indices)
            if min(indices)<1
                error('The indices vector cannot contain negative values');
            end
            obj.indices=indices;
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