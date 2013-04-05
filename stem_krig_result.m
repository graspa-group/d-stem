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
    
    %CONSTANTS
    %NN - number of kriging sites
    %T - number of temporal steps    
    
    properties
        variable_name=[];       %[string]           (1x1)    the name of the kriged variable
        stem_grid=[];           %[stem_grid object] (1x1)    a stem_grid object with the information of the kriging grid
        shape=[];               %[struct]           (1x1)    boundary of the geographic region loaded from a shape file
        
        y_hat=[];               %[double]           (NNxT)   the kriging result  
        diag_Var_y_hat=[];      %[double]           (NNxT)   the variance of the kriging result (only the variance and no covariance)
        E_wg_y1=[];             %[double]           (NNxTxK) the estimated latent variable w_g 
        diag_Var_wg_y1=[];      %[double]           (NNxTxK) the variance of the estimated latent variable w_g 
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
        
        %         function google_map(obj,type)
        %             if sum(strcmp(type,{'measure','std'}))==0
        %                 error('type must be either measure or std');
        %             end
        %             if strcmp(type,'measure')
        %                 datafile=[1,size(obj.y_hat,3),obj.stem_grid.pixel_side_w];
        %                 csvwrite('..\Data\google_bridge\parameters_kriging.csv',datafile);
        %                 for t=1:size(obj.y_hat,3)
        %                     t
        %                     min_value=nanmin(nanmin(obj.y_hat(:,:,t)));
        %                     max_value=nanmax(nanmax(obj.y_hat(:,:,t)));
        %                     data=obj.y_hat(:,:,t);
        %                     data=data(:);
        %                     value=(data-min_value)/(max_value-min_value);
        %                     %color=stem_krig_result.toColor(value);
        %                     datafile=[obj.stem_grid.coordinate(:,1),obj.stem_grid.coordinate(:,2),value,data];
        %                     csvwrite(['..\Data\google_bridge\kriging',num2str(t),'.csv'],datafile);
        %                 end
        %             else
        %                 datafile=[1,size(obj.var_y_hat,3),obj.stem_grid.pixel_side_w];
        %                 csvwrite('..\Data\google_bridge\parameters_kriging.csv',datafile);
        %                 for t=1:size(obj.var_y_hat,3)
        %                     t
        %                     min_value=nanmin(nanmin(sqrt(obj.var_y_hat(:,:,t))));
        %                     max_value=nanmax(nanmax(sqrt(obj.var_y_hat(:,:,t))));
        %                     data=sqrt(obj.var_y_hat(:,:,t));
        %                     data=data(:);
        %                     value=(data-min_value)/(max_value-min_value);
        %                     datafile=[obj.stem_grid.coordinate(:,1),obj.stem_grid.coordinate(:,2),value,data];
        %                     csvwrite(['..\Data\google_bridge\kriging',num2str(t),'.csv'],datafile);
        %                 end
        %                 winopen('..\Data\google_bridge\open_kriging.bat');
        %             end
        %         end
 
    end
end