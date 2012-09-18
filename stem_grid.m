%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef stem_grid
    
    properties
        coordinate=[];
        unit='';
        grid_size=[];
        grid_type=[];
        site_type=[];
        pixel_shape='';
        pixel_side_w=[];
        pixel_side_h=[];
        duplicated_sites=[];
    end
    
    properties (SetAccess = private)
        box=[];
    end
    
    methods
        function obj = stem_grid(coordinate,unit,grid_type,site_type,grid_size,pixel_shape,pixel_side_w,pixel_side_h)
            if nargin<4
                error('Not enough input arguments');
            end
            if strcmp(site_type,'pixel')&&(nargin<8)
                error('pixel_shape, pixel_side_w and pixel_side_h must be provided');
            end
            if nargin==4
                grid_size=[];
            end
            obj.grid_type=grid_type; %the other of this two lines is important and cannot change
            obj.coordinate=coordinate;
            obj.unit=unit;
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
       
        function google_map(obj)
            datafile=[obj.coordinate(:,1),obj.coordinate(:,2)];
            csvwrite('..\Data\google_bridge\gridplot.csv',datafile);
            winopen('..\Data\google_bridge\open_gridplot.bat');
        end
        
        %set functions
        function obj = set.coordinate(obj,coordinate)
            if not(size(coordinate,2)==2)
                error('coordinate must be a Nx2 matrix');
            end
            
%             if not(strcmp(obj.grid_type,'regular'))
%                 obj.duplicated_sites=[];
%                 for i=1:length(coordinate)-1
%                     temp=coordinate(i,:);
%                     temp2=coordinate((i+1):end,:);
%                     temp_lat=temp2(:,1);
%                     temp_lon=temp2(:,2);
%                     a=temp_lat==temp(1);
%                     b=temp_lon==temp(2);
%                     c=a&b;
%                     if sum(c)>0
%                         obj.duplicated_sites=[obj.duplicated_sites;find(c,1)+i];
%                         disp(['WARNING: coordinate ',num2str(i),' equal to coordinate ',num2str(find(c,1)+i)]);
%                     end
%                 end
%             end

            obj.coordinate=coordinate;
            obj.box(1)=min(coordinate(:,1));
            obj.box(2)=max(coordinate(:,1));
            obj.box(3)=min(coordinate(:,2));
            obj.box(4)=max(coordinate(:,2));
        end
        
        function obj = set.unit(obj,unit)
            if not(strcmp(unit,'deg') || strcmp(unit,'m') || strcmp(unit,'km'))
                error('unit must be ''deg'' or ''m'' or ''km''');
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
                    error('The grid_size must be a 2x1 vector');
                end
                if not(grid_size(1)*grid_size(2)==size(obj.coordinate,1))
                    error('The grid size is not compatible with the grid');
                end
            end
            obj.grid_size=grid_size;
        end
        
        function obj = set.pixel_shape(obj,pixel_shape)
            if not(strcmp(pixel_shape,'square') || strcmp(pixel_shape,'rectangular') || strcmp(pixel_shape,'irregular'))
                error('pixel_shape can be only ''square'', ''rectangular'' or ''irregular''');
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
end