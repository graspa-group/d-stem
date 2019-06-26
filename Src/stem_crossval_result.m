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



classdef stem_crossval_result < handle
    
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
    %dN - the number of cross-validation sites
    %T  - the number of time steps
    %H  - the dimension of vector h
    
    properties
        
        y_back=[];                  %[double]                  (dNxT) the original data (back-transformed) of the cross-validation sites
        y_hat_back=[];              %[double]                  (dNxT) the estimated data (back-transformed) for the cross-validation sites
        res=[];                     %[double]                  (dNxT) cross-validation residuals
        res_back=[];                %[double]                  (dNxT) back-transformed cross-validation residuals (if the original data are log-transformed and/or standardized)
        
        cv_mse_t=[];                %[double >=0]              (Tx1) mean squared error with respect to time t 
        cv_R2_t=[];                 %[double >=0]              (Tx1) R-square with respect to time t 
        cv_mse_s=[];                %[double >=0]              (dNx1) mean squared error with respect to site s 
        cv_R2_s=[];                 %[double >=0]              (dNx1) R-square with respect to time t 
        
        cv_h=[];                    %[double >=0]              (Hx1) the h values if f-HDGM
        cv_mse_h=[];                %[double >=0]              (Hx1) mean squared error with respect to h if f-HDGM
        cv_R2_h=[];                 %[double >=0]              (Hx1) R-square with respect to h if f-HDGM
        
        min_distance=[];            %[double >=0]              (dNx1) the distance from each cross-validation site to the nearest non cross-validation site
    end
    
    methods
        function obj = stem_crossval_result()
            %DESCRIPTION: object constructor
            %
            %INPUT
            %no inputs required
            %
            %
            %OUTPUT
            %obj           - [stem_crossval_result object] (1x1)
        end
        
    end
end