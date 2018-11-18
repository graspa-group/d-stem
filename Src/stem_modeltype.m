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



classdef stem_modeltype < handle
    
    properties
        model_name;             %[string]       (1x1) the name of the model either 'DCM' (Dynamic Coregionalization Model), 'HDGM' (Hidden Dynamic Geostatistical Model), 'f-HDGM' (functional Hidden Dynamic Geostatistical Model), 'MBC' (Model-based clustering) or 'Emulator'
        clustering_type;        %[string]       (1x1) the type of clusterig if model_name is 'MBC'. Can be either 'Monopoles' or 'Dipoles'.
        clustering_error_type;  %[string]       (1x1) the way the model parameter sigma_eps is computer. Can be either 'Shared', 'Proportional' or 'Dynamic'.
    end

    methods
        function obj = stem_modeltype(model_name)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %model_name             -   [string]    (1x1) the name of the model
            %
            %OUTPUT
            %obj                    -   [stem_modeltype object] (1x1)    
            
            if nargin<1
                error('Not enough input arguments');
            end
            
            obj.model_name=model_name;
            
            if strcmpi(model_name,'MBC')
                obj.clustering_type='Monopoles';
                obj.clustering_error_type='Shared';
            end
        end
        
        function res = is(obj,model_name)
            if iscell(model_name)
                res=false;
                for i=1:length(model_name)
                    if strcmpi(obj.model_name,model_name{i})
                        res=true;
                    end
                end
            else
                res = strcmpi(obj.model_name,model_name);
            end
        end
        
        function res = clustering_type_is(obj,clustering_type)
             res = strcmpi(obj.clustering_type,clustering_type);
        end
        
        function res = clustering_error_type_is(obj,clustering_error_type)
             res = strcmpi(obj.clustering_error_type,clustering_error_type);
        end
       
        %Class set methods
        function set.model_name(obj,model_name)
            if sum(strcmpi(model_name,{'DCM','HDGM','f-HDGM','MBC','Emulator'}))==0
                error('The model_name input argument must be either ''DCM'', ''HDGM'', ''f-HDGM'', ''MBC'' or ''Emulator''');
            end
            obj.model_name=model_name;
        end
        
        function set.clustering_type(obj,clustering_type)
            if sum(strcmpi(clustering_type,{'Monopoles','Dipoles'}))==0
                error('The clustering_type input argument must be either ''Monopoles'' or ''Dipoles''');
            end
            obj.clustering_type=clustering_type;
        end 
        
        function set.clustering_error_type(obj,clustering_error_type)
            if sum(strcmpi(clustering_error_type,{'Shared','Proportional','Dynamic'}))==0
                error('The clustering_error_type input argument must be either ''Shared'', ''Proportional'' or ''Dynamic''');
            end
            obj.clustering_error_type=clustering_error_type;
        end
    end
    
end