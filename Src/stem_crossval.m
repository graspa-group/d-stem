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

classdef stem_crossval < handle
    
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
    %dq - number of cross-validation variables
    %dN - number of cross-validation sites
    
    properties
        variable_name={};           %[string]                      {dqx1} the names of the cross-validation variables
        type={};                    %[string]                      {dqx1} 'point': cross-validation is on point-data
        indices={};                 %[integer >0]                  {dq}x(dNx1) the indices of the cross-validation sites 
        
        stem_varset={};             %[stem_varset objects]         {dqx1} the subset of data used for cross-validation
        stem_gridlist={};           %[stem_gridlist objects]       {dqx1} this object is needed for storing the coordinates of the cross-validation sites

        nn_size=0;                  %[integer>=0]                  (1x1) the size of the nearest neighbor set for each kriging site in cross-validation
    end

    methods
        function obj = stem_crossval(variable_name,indices,nn_size,type)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %variable_name - [string]               {1x1} the names of the cross-validation variables
            %indices       - [integer >0]           {dq}x(dNx1) the indices of the cross-validation sites for each variable
            %nn_size=0;    - [integer>=0]           (1x1) the size of the nearest neighbor set for each kriging site in cross-validation
            %type          - [string]               {dqx1} 'point': cross-validation is on point-data
            %
            %OUTPUT
            %obj           - [stem_crossval object] (1x1)    
            
            if nargin<2
                error('Not enough input arguments');
            end
            
            obj.variable_name=variable_name;
            obj.indices=indices;
            
            if nargin>2
                obj.nn_size=nn_size;
            end
            if nargin>3
                obj.type=type;            
            else
                obj.type='point';
            end  
        end
       
        %Class set methods
        function set.type(obj,type)
            if iscell(type)
                for i=1:length(type)
                    if not(strcmp(type{i},'point'))
                        error('The cross-validation is available only for point variables. The type must be equal to ''point.''');
                    end
                end
            else
                if not(strcmp(type,'point'))
                    error('The cross-validation is available only for point variables. The type must be equal to ''point.''');
                end
            end
            obj.type=type;
        end
        
        function set.indices(obj,indices)
            if iscell(indices)
                for i=1:length(indices)
                    if min(indices{i})<1
                        error('The indices vector cannot contain negative values');
                    end
                end
            else   
                if min(indices)<1
                    error('The indices vector cannot contain negative values');
                end
            end 
            obj.indices=indices;
        end         
        
        function set.nn_size(obj,nn_size)
            if nn_size>=0
                obj.nn_size=nn_size;
            else
                error('nn_size must be >=0');
            end
        end
    end
    
end