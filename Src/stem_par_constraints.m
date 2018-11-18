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

classdef stem_par_constraints
   
    properties
        time_diagonal=1;        %[boolean]	(1x1) 1: matrix G and sigma_eta are diagonal; 0: matrix G and sigma_eta are full
        pixel_correlated=0;       %[boolean]  (1x1) 1: when model_name is 'DCM', pixel variables are cross-correlated; 0: pixel variables are NOT cross-correlated
    end
    
    methods
        function obj = stem_par_constraints()
            %DESCRIPTION: object constructor
            %
            %INPUT
            %
            %OUTPUT
            %obj - [stem_par_constraints object] (1x1)
        end
        
        %Class set methods
        function obj = set.time_diagonal(obj,time_diagonal)
            if not(time_diagonal==0||time_diagonal==1)
                error('time_diagonal must be either 0 or 1');
            end        
            obj.time_diagonal=time_diagonal;
        end
        
        function obj = set.pixel_correlated(obj,pixel_correlated)
            if (pixel_correlated<0)||(pixel_correlated>1)
                error('The pixel_correlated input argument must be either 0 or 1');
            end
            obj.pixel_correlated=pixel_correlated;
        end
        
    end
end