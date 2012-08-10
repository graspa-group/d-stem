%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef stem_krig_result < handle
    
    properties
        stem_grid=[];
        variable_name=[];
        y_hat=[];
        var_y_hat=[];
        E_wg_y1=[];
        shape=[];
    end
    
    methods
        function obj = stem_krig_result(variable_name,stem_grid,shape)
            if nargin<3
                error('Not enough input arguments');
            end
            obj.variable_name=variable_name;
            obj.stem_grid=stem_grid;
            obj.shape=shape;
        end
        
        function google_map(obj,type)
            if sum(strcmp(type,{'measure','std'}))==0
                error('type must be either measure or std');
            end
            if strcmp(type,'measure')
                datafile=[1,size(obj.y_hat,3),obj.stem_grid.pixel_side_w];
                csvwrite('..\Data\google_bridge\parameters_kriging.csv',datafile);
                for t=1:size(obj.y_hat,3)
                    t
                    min_value=nanmin(nanmin(obj.y_hat(:,:,t)));
                    max_value=nanmax(nanmax(obj.y_hat(:,:,t)));                    
                    data=obj.y_hat(:,:,t);
                    data=data(:);
                    value=(data-min_value)/(max_value-min_value);
                    %color=stem_krig_result.toColor(value);
                    datafile=[obj.stem_grid.coordinate(:,1),obj.stem_grid.coordinate(:,2),value,data];
                    csvwrite(['..\Data\google_bridge\kriging',num2str(t),'.csv'],datafile);
                end
            else
                datafile=[1,size(obj.var_y_hat,3),obj.stem_grid.pixel_side_w];
                csvwrite('..\Data\google_bridge\parameters_kriging.csv',datafile);
                for t=1:size(obj.var_y_hat,3)
                    t
                    min_value=nanmin(nanmin(sqrt(obj.var_y_hat(:,:,t))));
                    max_value=nanmax(nanmax(sqrt(obj.var_y_hat(:,:,t))));
                    data=sqrt(obj.var_y_hat(:,:,t));
                    data=data(:);
                    value=(data-min_value)/(max_value-min_value);
                    %color=stem_krig_result.toColor(value);
                    datafile=[obj.stem_grid.coordinate(:,1),obj.stem_grid.coordinate(:,2),value,data];
                    csvwrite(['..\Data\google_bridge\kriging',num2str(t),'.csv'],datafile);
                end
                winopen('..\Data\google_bridge\open_kriging.bat');
            end
        end
        
        function plot(obj,time_step,type)
            if nargin<1
                error('Not enough input arguments');
            end
            if (time_step<0)||(time_step>size(obj.y_hat,3))
                error('time_step out of bounds');
            end

            lat=reshape(obj.stem_grid.coordinate(:,1),obj.stem_grid.grid_size(1),obj.stem_grid.grid_size(2));
            lon=reshape(obj.stem_grid.coordinate(:,2),obj.stem_grid.grid_size(1),obj.stem_grid.grid_size(2));
            
            if strcmp(type,'both')
                subplot(1,2,1);
            end
            if strcmp(type,'measure')||strcmp(type,'both')
                axesm('MapProjection','sinusoid','MapLatLimit',[44 54],'MapLonLimit',[7 14],...
                    'MLineLocation',1,'PLineLocation',1,'GColor',[0 0 0],'ParallelLabel','on','MeridianLabel','on','Grid','on')
                hold on
                %obj.shape.X(isnan(obj.shape.X))=0;
                %obj.shape.Y(isnan(obj.shape.Y))=0;
                %indices=find(obj.shape.X==0);
                
                %for j=1:length(indices)-1
                %    geoshow(obj.shape.Y(indices(j)+1:indices(j+1)-1),obj.shape.X(indices(j)+1:indices(j+1)-1),'DisplayType','polygon'); %,'FaceColor','none'
                %end
                %geoshow(obj.shape);
                if time_step>0
                    temp=obj.y_hat(:,:,time_step);
                else
                    temp=mean(obj.y_hat,3);
                end
                h = geoshow(lat,lon,temp,'DisplayType','texture');
                %set(h,'FaceColor','flat');
                colorbar
            end

            if strcmp(type,'both')
                subplot(1,2,2);
            end
            if strcmp(type,'variance')||strcmp(type,'both')
                axesm('MapProjection','sinusoid','MapLatLimit',[54.5 61],'MapLonLimit',[-8 0],...
                    'MLineLocation',1,'PLineLocation',1,'GColor',[0 0 0],'ParallelLabel','on','MeridianLabel','on','Grid','on')
                obj.shape.X(isnan(obj.shape.X))=0;
                obj.shape.Y(isnan(obj.shape.Y))=0;
                indices=find(obj.shape.X==0);
                
                for j=1:length(indices)-1
                    geoshow(obj.shape.Y(indices(j)+1:indices(j+1)-1),obj.shape.X(indices(j)+1:indices(j+1)-1),'DisplayType','polygon');
                end
                if time_step>0
                    temp=obj.var_y_hat(:,:,time_step);
                else
                    temp=sqrt(mean(obj.var_y_hat,3));
                end
                h = geoshow(lat,lon,temp,'DisplayType','texture');
                set(h,'FaceColor','flat');
                colorbar;
                colormap('gray');
            end
        end
    end
    
    methods (Static)
        function color = toColor(vector)
            if (min(vector)<0)||(max(vector)>1)
                error('The elements of vector must be between 0 and 1');
            end
            for i=1:length(vector)
                if not(isnan(vector(i)))
                    if vector(i)<0.25
                        val=round(vector(i)*4*255);
                        s=dec2hex(val);
                        if length(s)==1
                            s=['0',s];
                        end
                        string=['#00',s,'FF'];
                        color(i,:)=['#00',dec2hex(val),'FF'];
                    end
                    if (vector(i)>=0.25)&&(vector(i)<0.50)
                        val=round((1-(vector(i)-0.25)*4)*255);
                        s=dec2hex(val);
                        if length(s)==1
                            s=['0',s];
                        end                        
                        string=['#00FF',s];
                        color(i,:)=string;
                    end
                    if (vector(i)>=0.50)&&(vector(i)<0.75)
                        val=round((vector(i)-0.5)*4*255);
                        s=dec2hex(val);
                        if length(s)==1
                            s=['0',s];
                        end                        
                        string=['#',s,'FF00'];
                        color(i,:)=string;
                    end
                    if (vector(i)>=0.75)
                        val=round((1-(vector(i)-0.75)*4)*255);
                        s=dec2hex(val);
                        if length(s)==1
                            s=['0',s];
                        end                        
                        string=['#FF',s,'00'];
                        color(i,:)=string;
                    end
                else
                    color(i,:)='#000000';
                end
            end
        end
    end
    
end