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

classdef stem_grid
    
    %CONSTANTS
    %
    %d    - the dimension of the space
    %ni_g - the number of point sites for the i-th variable
    %ni_r - the number of pixel sites for the i-th variable
    
    
    properties
        coordinate=[];              %[double]    (ni_g|ni_rxd)  matrix of spatial coordinates composed of either x_1,...,x_p coordinates or latitude and longitude
        unit='';                    %[string]    (1x1)          unit of measure of the coordinates. 'deg': degree; 'km': kilometers; 'm': meters; 'none': none (to use when model_name is 'emulator');
        grid_type=[];               %[string]    (1x1)          'sparse': sparse spatial locations; 'regular': spatial locations at fixed intervals in space
        grid_size=[];               %[integer >0](1xd)          if grid_type='regular', then grid_size is the number of elements in each of the d dimensions of the grid
        site_type=[];               %[string]    (1x1)          'point': the spatial coordinates relate to point data; 'pixel': the spatial coordinates related to pixel data
        pixel_shape='';             %[string]    (1x1)          'square': square pixels; 'rectangular': rectangular pixels
        pixel_side_w=[];            %[double]    (1x1)          the width of the pixel (in the same unit of measure of the unit property)
        pixel_side_h=[];            %[double]    (1x1)          the height of the pixel (in the same unit of measure of the unit property)
        duplicated_sites=[];        %[integer]   (hx1)          indices of the duplicated spatial locations in the coordinate property 
    end
    
    properties (Dependent, SetAccess = private)
       box=[];                      %[double]    (4x1) the bounding box of the geographic area covered by the grid [lat_min,lat_max,lon_min,lon_max] 
    end
    
    methods
        function obj = stem_grid(coordinate,unit,grid_type,site_type,grid_size,pixel_shape,pixel_side_w,pixel_side_h)
            %DESCRIPTION: is the constructor of the class stem_data
            %
            %INPUT
            %
            %coordinate         - [double]              (ni_g|ni_rx2)  matrix of spatial coordinates composed of either x and y vectors or latitude and longitude vectors
            %unit='';           - [string]              (1x1)          unit of measure of the coordinates. 'deg': degree; 'km': kilometers; 'm': meters; 'none': none;
            %grid_type=[];      - [string]              (1x1)          'sparse': sparse spatial locations; 'regular': spatial locations at fixed intervals in space
            %grid_size=[];      - [integer >0]          (1xd)          if grid_type='regular', then grid_size is the number of elements in each dimension of the grid
            %site_type=[];      - [string]              (1x1)          'point': the spatial coordinates relate to point data; 'pixel': the spatial coordinates related to pixel data
            %pixel_shape='';    - [string]              (1x1)          'square': square pixels; 'rectangular': rectangular pixels
            %pixel_side_w=[];   - [double]              (1x1)          the width of the pixel (in the same unit of measure of the unit property)
            %pixel_side_h=[];   - [double]              (1x1)          the height of the pixel (in the same unit of measure of the unit property)
            %
            %OUTPUT
            %obj                - [stem_grid object]    (1x1)
            if nargin<4
                error('Not enough input arguments');
            end
            if strcmp(site_type,'pixel')&&(nargin<8)
                error('pixel_shape, pixel_side_w and pixel_side_h must be provided');
            end
            if nargin==4
                grid_size=[];
            end
            
            %the other of the following 3 lines is important and cannot change
            obj.unit=unit;   
            obj.grid_type=grid_type; 
            obj.coordinate=coordinate;
            
            obj.site_type=site_type;
            obj.grid_size=grid_size;
            if nargin>5
                if not(strcmp(site_type,'pixel'))
                    warning('pixel_shape in ignored');
                else
                    obj.pixel_shape=pixel_shape;
                end
            end
            if nargin==7
                error('Also pixel_side_h must be provided');
            end
            if nargin>=8
                if not(strcmp(site_type,'pixel'))
                    warning('pixel_side is ignored');
                else
                    obj.pixel_side_w=pixel_side_w;
                    obj.pixel_side_h=pixel_side_h;
                end                
            end
        end
        
        function obj = permute(obj,indices)
            obj.coordinate=obj.coordinate(indices,:);
        end
        
        function obj = sorted_by_lat(obj)
            %DESCRIPTION: reorder coordinates with respect to the first coordinate
            %
            %INPUT
            %none
            %
            %OUTPUT
            %
            %none
            [~,idx]=sort(obj.coordinate(:,1));
            obj.coordinate=obj.coordinate(idx,:);
        end
        
        % Class set methods
        function box = get.box(obj)
            if not(isempty(obj.coordinate))
                box=zeros(size(obj.coordinate,2)*2,1);
                counter=1;
                for i=1:size(obj.coordinate,2)
                    box(counter)=min(obj.coordinate(:,i));
                    counter=counter+1;
                    box(counter)=max(obj.coordinate(:,i));
                    counter=counter+1;
                end
            else
                box=0;
            end
        end
        
        function obj = set.coordinate(obj,coordinate)
            if strcmp(obj.unit,'deg')
                if not(size(coordinate,2)==2)
                    error('coordinate must be a Nx2 matrix since unit is ''deg''');
                end
            end
            
            if not(strcmp(obj.grid_type,'regular'))
                if length(coordinate)<100000
                    obj.duplicated_sites=[];
                    for i=1:size(coordinate,1)-1
                        temp=coordinate(i,:);
                        temp2=coordinate((i+1):end,:);
                        temp=repmat(temp,[size(temp2,1),1]);
                        c=all(temp==temp2,2);
                        if sum(c)>0
                            obj.duplicated_sites=[obj.duplicated_sites;find(c,1)+i];
                            disp(['WARNING: coordinate ',num2str(i),' equal to coordinate ',num2str(find(c,1)+i)]);
                        end
                    end
                else
                    disp(['WARNING: too many spatial locations. The test for duplicate spatial locations is skipped']);
                end
            end

            obj.coordinate=coordinate;
        end
        
        function obj = set.unit(obj,unit)
            if not(strcmp(unit,'deg') || strcmp(unit,'m') || strcmp(unit,'km') || strcmp(unit,'none'))
                error('unit must be ''deg'' or ''m'' or ''km'' or ''none''');
            end
            
            obj.unit=unit;
        end
        
        function obj = set.grid_type(obj,grid_type)
            if not(strcmp(grid_type,'regular')||strcmp(grid_type,'sparse'))
                error('The grid type must be either regular or sparse');
            end
            obj.grid_type=grid_type;
        end
        
        function obj = set.site_type(obj,site_type)
            if not(strcmp(site_type,'point')||strcmp(site_type,'pixel'))
                error('The grid type must be either point or pixel');
            end
            obj.site_type=site_type;            
        end
        
        function obj = set.grid_size(obj,grid_size)
            if (not(isempty(grid_size)))&&(strcmp(obj.grid_type,'sparse'))
                error('The grid size must be provided only for regular grid');
            end
            if not(isempty(grid_size))
                if not(isvector(grid_size))
                    error('The grid_size must be a 1xd vector with d the number of grid dimensions');
                end
                if not(size(grid_size,1)==1)
                    error('grid_size must be a 1xd vector');
                end
                if not(prod(grid_size)==size(obj.coordinate,1))
                    error('The grid size is not compatible with the grid');
                end
            end
            obj.grid_size=grid_size;
        end
        
        function obj = set.pixel_shape(obj,pixel_shape)
            if not(strcmp(pixel_shape,'square') || strcmp(pixel_shape,'rectangular'))
                error('pixel_shape can be only ''square'' or ''rectangular''');
            end
            obj.pixel_shape=pixel_shape;
        end
        
        function obj = set.pixel_side_w(obj,pixel_side_w)
            if pixel_side_w<0
                error('pixel_side_w must be >0');
            end
            obj.pixel_side_w=pixel_side_w;
        end
        
        function obj = set.pixel_side_h(obj,pixel_side_h)
            if pixel_side_h<0
                error('pixel_side_h must be >0');
            end
            obj.pixel_side_h=pixel_side_h;
        end 
    end
    
    methods (Static)
        
    end
        
end