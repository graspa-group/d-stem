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

classdef stem_krig_result < handle
    
    %CONSTANTS
    %NN - number of kriging sites
    %T - number of temporal steps    
    
    properties
        variable_name=[];       %[string]                (1x1) the name of the kriged variable
        stem_grid=[];           %[stem_grid object]      (1x1) a stem_grid object with the information of the kriging grid
        stem_datestamp=[];      %[stem_datestamp object] (1x1) a stem_datestamp object with te information on the date stamps
        shape=[];               %[struct]                (1x1) boundary of the geographic region loaded from a shape file
        
        y_hat=[];               %[double]                (NNxT)   the kriging result  
        diag_Var_y_hat=[];      %[double]                (NNxT)   the variance of the kriging result (only the variance and no covariance)
        E_wp_y1=[];             %[double]                (NNxTxK) the estimated latent variable w_g 
        diag_Var_wp_y1=[];      %[double]                (NNxTxK) the variance of the estimated latent variable w_g 
    end
    
    methods
        function obj = stem_krig_result(variable_name,stem_grid,shape)
            %DESCRIPTION: the constructor of the class stem_krig_result
            %
            %INPUT
            %
            %variable_name      - [string]                (1x1) the name of the kriged variable
            %stem_grid          - [stem_grid object]      (1x1) a stem_grid object with the information of the kriging grid
            %<shape>            - [stem_varset object]    (1x1) (default: []) boundary of the geographic region loaded from a shape file
            %
            %OUTPUT
            %obj                - [stem_krig_result object]  (1x1) the stem_krig_result object
                        
            if nargin<2
                error('Not enough input arguments');
            end
            obj.variable_name=variable_name;
            obj.stem_grid=stem_grid;
            if nargin>2
                obj.shape=shape;
            end
        end
        
        function fig_h = plot(obj,time_step,type)
            %DESCRIPTION: plot the kriging result
            %
            %INPUT
            %obj                - [stem_krig_result object]  (1x1) the stem_krig_result object 
            %time_step          - [integer >=0]              (1x1) the time step to plot. If time_step=0 the average with respect to time is plotted
            %<type>             - [string]                   (1x1)(default: 'both') 'variable': only the kriged variable is plotted; 'std': only the standard deviation of the kriged variable is plotted; 'both': both of the previous are plotted
            %
            %OUTPUT
            %fig_h              - [integer]                  (1x1) the handle of the figure
            if nargin<2
                error('Not enough input arguments');
            end
            if nargin<3
                type='both';
            end
            if not(strcmp(type,'variable')||strcmp(type,'std')||strcmp(type,'both'))
                error('type can be either ''variable'', ''std'' or ''both''');
            end
            if strcmp(type,'std')&&isempty(obj.diag_Var_y_hat)
                error('The variance of the kriged variabled has not been evaluated and cannot be plotted');
            end
            if strcmp(type,'both')&&isempty(obj.diag_Var_y_hat)
                type='variable';
                disp('WARNING: the variance of the kriged variabled has not been evaluated and cannot be plotted. Only the variable is plotted');
            end

            if (time_step<0)||(time_step>size(obj.y_hat,3))
                error('time_step out of bounds');
            end

            lat=reshape(obj.stem_grid.coordinate(:,1),obj.stem_grid.grid_size(1),obj.stem_grid.grid_size(2));
            lon=reshape(obj.stem_grid.coordinate(:,2),obj.stem_grid.grid_size(1),obj.stem_grid.grid_size(2));
            
            if strcmp(type,'both')
                subplot(2,1,1);
            end
            if strcmp(type,'variable')||strcmp(type,'both')
                hold on
                if not(isempty(obj.shape))
                    mapshow(obj.shape);
                end
                if time_step>0
                    temp=obj.y_hat(:,:,time_step);
                    title([obj.variable_name,' - ',datestr(obj.stem_datestamp.stamp(time_step))]);
                else
                    temp=mean(obj.y_hat,3);
                    title(['Average ',obj.variable_name,' from ',datestr(obj.stem_datestamp.stamp(1)),' to ',datestr(obj.stem_datestamp.stamp(end))]);
                end
                h = mapshow(lon,lat,temp,'DisplayType','texture');
                set(h,'FaceColor','flat');
                axis equal
                xlim([min(lon(:)),max(lon(:))]);
                ylim([min(lat(:)),max(lat(:))]);
                if strcmp(obj.stem_grid.unit,'deg')
                    xlabel('Longitude');
                    ylabel('Latitude');
                else
                    xlabel(obj.stem_grid.unit);
                    ylabel(obj.stem_grid.unit);
                end
                colorbar;
                grid on;
            end

            if strcmp(type,'both')
                subplot(2,1,2);
            end
            if strcmp(type,'std')||strcmp(type,'both')
                hold on
                if not(isempty(obj.shape))
                    mapshow(obj.shape);
                end
                if time_step>0
                    temp=sqrt(obj.diag_Var_y_hat(:,:,time_step));
                    title(['Std of ',obj.variable_name,' - ',datestr(obj.stem_datestamp.stamp(time_step))]);
                else
                    temp=mean(sqrt(obj.diag_Var_y_hat),3);
                    title(['Average std of',obj.variable_name,' from ',datestr(obj.stem_datestamp.stamp(1)),' to ',datestr(obj.stem_datestamp.stamp(end))]);
                end
                h = mapshow(lon,lat,temp,'DisplayType','texture');
                set(h,'FaceColor','flat');
                axis equal
                xlim([min(lon(:)),max(lon(:))]);
                ylim([min(lat(:)),max(lat(:))]);
                if strcmp(obj.stem_grid.unit,'deg')
                    xlabel('Longitude');
                    ylabel('Latitude');
                else
                    xlabel(obj.stem_grid.unit);
                    ylabel(obj.stem_grid.unit);
                end                
                colorbar;
                grid on;
            end
            if nargout>0
                fig_h=h;
            end
        end
    end
end