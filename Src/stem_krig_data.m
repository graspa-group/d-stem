%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% D-STEM - Distributed Space Time Expecation Maximization              %
%%%                                                                      %
%%% Author: Francesco Finazzi                                            %
%%% E-mail: francesco.finazzi@unibg.it                                   %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Economics and Quantitative Methods %
%%% Author website: http://www.unibg.it/pers/?francesco.finazzi          %
%%% Author: Yaqiong Wang                                                 %
%%% E-mail: yaqiongwang@pku.edu.cn                                       %
%%% Affiliation: Peking University,                                      %
%%%              Guanghua school of management,                          %
%%%              Business Statistics and Econometrics                    %
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

classdef stem_krig_data < handle
    properties
        grid=[];                    %[stem_grid object] (1x1)   a stem_grid object for the kriging coordinates
        X=[];                       %[double]           (NxBxT) a matrix with all the loading coefficients used in model estimation. The dimension is NxBxT where N is the number of kriging locations and B the number of loading coefficients
        X_names={};                 %[string]           {Bx1}   a cell-array of covariate names. Must be the same used in model estimation
        mask=[];                    %[integer]          (Nx1)   a mask for the subset of spatial locations where the kriging is needed. Must be a vector of NaNs and 1s of the same lenght of the vector of coordinates
    end
    
    methods
        function obj = stem_krig_data(grid,X,X_names,mask)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %grid              %[stem_grid object] (1x1)   a stem_grid object for the kriging coordinates
            %X                 %[double]           (NxBxT) a matrix with all the loading coefficients used in model estimation. The dimension is NxBxT where N is the number of kriging locations and B the number of loading coefficients
            %X_names           %[string]           {Bx1}   a cell-array of covariate names. Must be the same used in model estimation
            %mask              %[integer]          (Nx1)   a mask for the subset of spatial locations where the kriging is needed. Must be a vector of NaNs and 1s of the same lenght of the vector of coordinates
            %
            %OUTPUT
            %obj       - [stem_krig_options] (1x1)
            if nargin<1
                error('The input arguments grid must be provide');
            end
            obj.grid=grid;
            
            if nargin>1
                if nargin<3
                    error('X_names must be provided');
                end
                if not(isempty(X))
                    if not(size(X,1)==size(grid.coordinate,1))
                        error('The number of rows of X must be equal to the number of coordinates in the stem_grid object');
                    end
                    obj.X=X;
                end
                if not(isempty(X_names))
                    if isempty(X)
                        error('X_names does not have to be provided when X is empty');
                    end
                    if not(size(X,2)==length(X_names))
                        error('The length of X_names must be equal to the second dimension of X');
                    end
                    obj.X_names=X_names;
                end
            end
            
            if nargin>3
                if not(isempty(mask))
                    if not(length(mask)==size(grid.coordinate,1))
                        error('mask must be a column vector with length equal to kriging locations');
                    end
                    if not(strcmp(grid.grid_type,'regular'))
                        error('mask can be provided only in the case of a regular kriging grid');
                    end
                    obj.mask=mask;
                end
            end
        end
        
        function time_average(obj,n_steps)
            %DESCRIPTION: computes time averages of n_steps for the matrix X
            %
            %INPUT
            %obj        - [stem_krig_data object]   (1x1) the stem_krig_data object
            %n_steps    - [integer >0]              (1x1) the number of temporal steps to average
            %
            %OUTPUT
            %
            %none: matrix X is updated
            
            if nargin<2
                error('n_steps must be provided');
            end
            if n_steps<=1
                error('n_steps must be greater than 1');
            end
            if round(n_steps)~=n_steps
                error('n_steps must be an integer value');
            end
            if isempty(obj.X)
                error('X is empty. Computing time averages is not possible');
            end
            if n_steps>size(obj.X,3)
                error('n_steps cannot be higher than the number of time steps in X');
            end
            
            disp('Time averaging of kriging matrix X started...');
            indices=0:n_steps:size(obj.X,3);
            
            if not(indices(end)==size(obj.X,3))
                disp(['The last ',num2str(size(obj.X,3)-indices(end)),' time steps will be discarded when taking time averages']);
            end
            
            X_temp=zeros(size(obj.X,1),size(obj.X,2),length(indices)-1);
            for j=1:length(indices)-1
                X_temp(:,:,j)=nanmean(obj.X(:,:,indices(j)+1:indices(j+1)),3);
            end
            obj.X=X_temp;
            
            disp('Time averaging of kriging matrix X ended.');
        end
        
        %Class set methods
        function obj = set.grid(obj,grid)
            if not(isa(grid,'stem_grid'))
                error('grid must be of class stem_grid');
            end
            obj.grid=grid;
        end
        
        function obj = set.mask(obj,mask)
            if size(mask,2)>1
                error('mask must be a column vector');
            end
            n_ones=sum(mask==1);
            n_nans=sum(isnan(mask));
            if (n_ones+n_nans<length(mask))
                error('All the elements of mask must be either 1 or NaN');
            end
            obj.mask=mask;
        end
        
        function obj = set.X(obj,X)
            if sum(isnan(X(:)))>1
                error('X cannot include NaNs');
            end 
            obj.X=X;
        end
        
        function obj = set.X_names(obj,X_names)
            if not(iscell(X_names))
                error('X_names must be a cell-array');
            end
            obj.X_names=X_names;
        end

    end
end