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

classdef stem_krig_result < handle
    
    %CONSTANTS
    %NN - number of kriging sites
    %T - number of temporal steps    
    %C - number of basis functions (when the model is f-HDGM)
    
    properties
        variable_name=[];       %[string]                (1x1) the name of the kriged variable
        stem_grid=[];           %[stem_grid object]      (1x1) a stem_grid object with the information of the kriging grid
        stem_datestamp=[];      %[stem_datestamp object] (1x1) a stem_datestamp object with te information on the date stamps
        shape=[];               %[struct]                (1x1) boundary of the geographic region loaded from a shape file
        stem_grid_sites=[];     %[stem_grid object]      (1x1) a stem_grid object with the sampling sites of the kriged variable
        
        y_hat=[];               %[double]                (NNxT)   the kriging result  
        diag_Var_y_hat=[];      %[double]                (NNxT)   the variance of the kriging result (only the variance and no covariance)
        E_wp_y1=[];             %[double]                (NNxTxK) the estimated latent variable w_g 
        diag_Var_wp_y1=[];      %[double]                (NNxTxK) the variance of the estimated latent variable w_g 
        zk_s=[];                %[double]                (NNxTxC) the kalman smoother output (only when model is f-HDGM)
        diag_Pk_s=[];           %[double]                (NNxTxC) the variance of the kalman smoother output (only when model is f-HDGM)
        
        coord_output_block={};  %[integer>0]             {B}(NBx2) a cell array with the coordinates of the output blocks
        coord_cond_block={}     %[integer>0]             {B}(NBx2) a cell array with the coordinates of the conditioning blocks
    
        %isSplineCoeff=0;        %[boolean]               (1x1) 0: the stem_krig_result object includes classic kriging results; 1: the stem_krig_result object includes the kriging of the spline coefficients when model_type is f-HDGM
       
    end
    %Yaqiong
    properties (Hidden = true)
        %alpha_z=[];             %[double]                (px1) alpha parameters related to the z latent variable when model_name is 'HDGM' or 'f-HDGM'
        %beta=[];                %[double]                (N_bx1) beta parameters
        stem_fda=[];                 %[stem_fda object]       (1x1) the stem_fda object of an estimated model
        stem_par=[];                 %[stem_par object]        (1x1) the stem_par object of an estimated model
        %stem_data=[];                %[stem_data object]      (1x1) the stem_data object of an estimated model
    end
    
    methods
        function obj = stem_krig_result(variable_name,stem_grid,stem_grid_sites,shape)
            %DESCRIPTION: the constructor of the class stem_krig_result
            %
            %INPUT
            %
            %variable_name      - [string]                (1x1) the name of the kriged variable
            %stem_grid          - [stem_grid object]      (1x1) a stem_grid object with the information of the kriging grid
            %stem_grid_sites    - [stem_grid object]      (1x1) a stem_grid object with the sampling sites of the kriged variable
            %<shape>            - [stem_varset object]    (1x1) (default: []) boundary of the geographic region loaded from a shape file
            %
            %OUTPUT
            %obj                - [stem_krig_result object]  (1x1) the stem_krig_result object
                        
            if nargin<3
                error('Not enough input arguments');
            end
            obj.variable_name=variable_name;
            obj.stem_grid=stem_grid;
            obj.stem_grid_sites=stem_grid_sites;
            if nargin>3
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
                if time_step>0
                    temp=obj.y_hat(:,:,time_step);
                    title([obj.variable_name,' - ',datestr(obj.stem_datestamp.stamp(time_step))],'FontSize',16);
                else
                    temp=mean(obj.y_hat,3);
                    title(['Average ',obj.variable_name,' from ',datestr(obj.stem_datestamp.stamp(1)),' to ',datestr(obj.stem_datestamp.stamp(end))],'FontSize',16);
                end
                h = mapshow(lon,lat,temp,'DisplayType','texturemap');
                set(h,'FaceColor','flat');
                if not(isempty(obj.shape))
                    mapshow(obj.shape,'FaceColor','none');
                end
                axis equal
                xlim([min(lon(:)),max(lon(:))]);
                ylim([min(lat(:)),max(lat(:))]);
                if strcmp(obj.stem_grid.unit,'deg')
                    xlabel('Longitude','FontSize',16);
                    ylabel('Latitude','FontSize',16);
                else
                    xlabel(obj.stem_grid.unit);
                    ylabel(obj.stem_grid.unit);
                end
                mapshow(obj.stem_grid_sites.coordinate(:,2),obj.stem_grid_sites.coordinate(:,1),'DisplayType','multipoint','Marker','*','MarkerEdgeColor','w');
                colormap jet;
                colorbar;
                grid on;
                box on;
                set(gca,'FontSize',16);
                set(gcf, 'renderer', 'zbuffer');
            end

            if strcmp(type,'both')
                subplot(2,1,2);
            end
            if strcmp(type,'std')||strcmp(type,'both')
                hold on
                if time_step>0
                    temp=sqrt(obj.diag_Var_y_hat(:,:,time_step));
                    title(['Std of ',obj.variable_name,' - ',datestr(obj.stem_datestamp.stamp(time_step))],'FontSize',16);
                else
                    temp=mean(sqrt(obj.diag_Var_y_hat),3);
                    title(['Average std of ',obj.variable_name,' from ',datestr(obj.stem_datestamp.stamp(1)),' to ',datestr(obj.stem_datestamp.stamp(end))],'FontSize',16);
                end
                h = mapshow(lon,lat,temp,'DisplayType','texturemap');
                set(h,'FaceColor','flat');
                if not(isempty(obj.shape))
                    mapshow(obj.shape,'FaceColor','none');
                end
                axis equal
                xlim([min(lon(:)),max(lon(:))]);
                ylim([min(lat(:)),max(lat(:))]);
                if strcmp(obj.stem_grid.unit,'deg')
                    xlabel('Longitude','FontSize',16);
                    ylabel('Latitude','FontSize',16);
                else
                    xlabel(obj.stem_grid.unit,'FontSize',16);
                    ylabel(obj.stem_grid.unit,'FontSize',16);
                end    
                mapshow(obj.stem_grid_sites.coordinate(:,2),obj.stem_grid_sites.coordinate(:,1),'DisplayType','multipoint','Marker','*','MarkerEdgeColor','w');
                colormap jet;
                colorbar;
                grid on;
                box on;
                set(gca,'FontSize',16);
                set(gcf, 'renderer', 'zbuffer');
            end
            if nargout>0
                fig_h=h;
            end
        end
        
        function [y_hat,diag_Var_y_hat]=surface_plot(obj,f,t,X_beta)
            %DESCRIPTION: plot the kriging result in the surface
            %
            %INPUT
            %obj                - [stem_krig_result object]  (1x1) the stem_krig_result object 
            %f                  - [integer >=0]              (1x1) the domain of functional data to predict the the kriged variable, with 0 indicating we want the average over the f at specific date t
            %t                  - [integer >=0]              (1x1) the date when we get the kriged variable, with 0 indicating we want average over the period at specific f
            %X_beta             - [double]                   (NxMxTxp) the value of X_beta at f
            %
            %OUTPUT
            %y_hat              - [integer]                  (1x1) the kriged variable
            %diag_Var_y_hat     - [integer]                  (1x1) the standard deviation of the kriged variable
            % the kriged variable is plotted; and the standard deviation of the kriged variable is plotted
            
            if isempty(obj.stem_fda)
                error('This function is only available for the model type f-HDGM ')
            end
            %{
            if not(isempty(f_range))&&f~=0
                warning('The f_range is not used when f is not 0')
            end
            %}
            if f==0 && t==0
                error('f and t cannot be = 0 at the same time (out of memory)');
            end
            if obj.stem_fda.flag_beta_spline==1
                k=obj.stem_fda.spline_nbasis_beta;
                %k=getnbasis(obj.stem_fda.spline_basis_beta);
                if size(X_beta,4)~=(length(obj.stem_par.beta)/k)
                    warning('You must provide the p dimension covariates X_beta ')
                end
            else
                if size(X_beta,4)~=length(obj.stem_par.beta)
                    warning('You must provide the p dimension covariates X_beta ')
                end
            end
            
            %X_beta=nan(size(z,1),size(z,2),t,p);
            %X_beta=ones(size(z,1),size(z,2),size(z,3),p);
            %alpha_z=obj.stem_par.alpha_z;
            %alpha_z_2=alpha_z.^2;
          
            z=obj.zk_s;
            var_z=obj.diag_Pk_s; 
            if length(size(X_beta))~=length(size(z))
                error('You must provide the X_beta(h) with 4D(N*M*T*p)')
            end
            %Yaqiong!  
            if obj.stem_fda.flag_beta_spline==1
                q = size(X_beta,4);
                k=getnbasis(obj.stem_fda.spline_basis_beta);
                X_beta_tmp = nan(size(z,1),size(z,2),size(z,3),q*k);
                basis_beta = obj.stem_fda.spline_basis_beta;
                b_beta=full(eval_basis(f,basis_beta));
                for lat_i=1:obj.stem_grid.grid_size(1)
                    for lon_i = 1:obj.stem_grid.grid_size(2)
                        for i = 1:q
                            X_beta_tmp(lat_i,lon_i,:,(i-1)*k+(1:k)) = squeeze(X_beta(lat_i,lon_i,:,i)).*b_beta;
                        end
                    end
                end 
                X_beta = X_beta_tmp;
            end
            
     
            basis = obj.stem_fda.spline_basis_z;
            b=full(eval_basis(f,basis));
            s=nan(size(z,1),size(z,2),size(z,3));
            v=nan(size(z,1),size(z,2),size(z,3));
            Xbeta = nan(size(z,1),size(z,2),size(z,3));
            y_hat = nan(size(z,1),size(z,2),size(z,3));
            for tt=1:size(z,3)
                for i=1:size(z,1)
                    %s(i,:,tt)=(squeeze(z(i,:,tt,:)).*repmat(alpha_z',[size(z,2),1]))*b';
                    s(i,:,tt)=squeeze(z(i,:,tt,:))*b';
                    %v(i,:,tt)=(squeeze(var_z(i,:,tt,:)).*repmat(alpha_z_2',[size(z,2),1]))*(b.^2)';
                    v(i,:,tt)=squeeze(var_z(i,:,tt,:))*(b.^2)';
                    Xbeta(i,:,tt)=squeeze(X_beta(i,:,tt,:))*obj.stem_par.beta;  
                    y_hat(i,:,tt)=Xbeta(i,:,tt)+s(i,:,tt);
                end
            end
            diag_Var_y_hat = v;
           
            %figure
            lat=reshape(obj.stem_grid.coordinate(:,1),obj.stem_grid.grid_size(1),obj.stem_grid.grid_size(2));
            lon=reshape(obj.stem_grid.coordinate(:,2),obj.stem_grid.grid_size(1),obj.stem_grid.grid_size(2));
            subplot(1,2,1);   
            hold on
            if t>0
                temp=y_hat(:,:,t);
                title([obj.variable_name,' - ',datestr(obj.stem_datestamp.stamp(t))],'FontSize',16);
            else
                temp=mean(y_hat,3);
                title({['Average ',obj.variable_name{:}, '@domain ', num2str(f)] ;[' from ',datestr(obj.stem_datestamp.stamp(1)),' to ',datestr(obj.stem_datestamp.stamp(end))]},'FontSize',16);
                 %temp1=strsplit(obj.variable_name,'_');
                %title(['Average ',temp1{1},'_{',temp1{2},'}'],'FontSize',16)
                %' from ',datestr(obj.stem_datestamp.stamp(1)),' to ',datestr(obj.stem_datestamp.stamp(end))],'FontSize',16);
            end
            h = mapshow(lon,lat,temp,'DisplayType','texturemap');
            set(h,'FaceColor','flat');
            if not(isempty(obj.shape))
                mapshow(obj.shape);
            end
            axis equal
            xlim([min(lon(:)),max(lon(:))]);
            ylim([min(lat(:)),max(lat(:))]);
            if strcmp(obj.stem_grid.unit,'deg')
                xlabel('Longitude','FontSize',16);
                ylabel('Latitude','FontSize',16);
            else
                xlabel(obj.stem_grid.unit);
                ylabel(obj.stem_grid.unit);
            end
            %Yaqiong
            %mapshow(c(:,2),c(:,1),'DisplayType','point');
            mapshow(obj.stem_grid_sites.coordinate(:,2),obj.stem_grid_sites.coordinate(:,1),'DisplayType','multipoint','Marker','*','MarkerEdgeColor','w');
            colormap jet;
            colorbar;
            grid on;
            box on;
            set(gca,'FontSize',16);
            set(gcf, 'renderer', 'zbuffer');
            
            subplot(1,2,2);
            hold on
            if t>0
                temp=sqrt(diag_Var_y_hat(:,:,t));
                title(['Std of ',obj.variable_name,' - ',datestr(obj.stem_datestamp.stamp(t))],'FontSize',16);
            else
                temp=mean(sqrt(diag_Var_y_hat),3);
                title({['Average std of ',obj.variable_name{:}, '@domain ', num2str(f)] ;[' from ',datestr(obj.stem_datestamp.stamp(1)),' to ',datestr(obj.stem_datestamp.stamp(end))]},'FontSize',16);
                %title({['Average std of ',obj.variable_name{:}, '@', num2str(f)] ;[' from ',datestr(obj.stem_datestamp.stamp(1)),' to ',datestr(obj.stem_datestamp.stamp(end)+1/24)]},'FontSize',16);
                %temp1=strsplit(obj.variable_name,'_');
                %title(['Average ',temp1{1},'_{',temp1{2},'}'],'FontSize',16)
                %title(['Average std of ',obj.variable_name,' from ',datestr(obj.stem_datestamp.stamp(1)),' to ',datestr(obj.stem_datestamp.stamp(end))],'FontSize',16);
            end
            h = mapshow(lon,lat,temp,'DisplayType','texturemap');
            set(h,'FaceColor','flat');
            if not(isempty(obj.shape))
                mapshow(obj.shape);
            end
            axis equal
            xlim([min(lon(:)),max(lon(:))]);
            ylim([min(lat(:)),max(lat(:))]);
            if strcmp(obj.stem_grid.unit,'deg')
                xlabel('Longitude','FontSize',16);
                ylabel('Latitude','FontSize',16);
            else
                xlabel(obj.stem_grid.unit,'FontSize',16);
                ylabel(obj.stem_grid.unit,'FontSize',16);
            end    
            mapshow(obj.stem_grid_sites.coordinate(:,2),obj.stem_grid_sites.coordinate(:,1),'DisplayType','multipoint','Marker','*','MarkerEdgeColor','w');
            colormap jet;
            colorbar;
            grid on;
            box on;
            set(gca,'FontSize',16);
            set(gcf, 'renderer', 'zbuffer');
        end
        
        function [y_hat,diag_Var_y_hat] = profile_plot(obj,f_range,lon0,lat0,t,X_beta)
           
            alpha_z=obj.stem_par.alpha_z;
            alpha_z_2=alpha_z.^2;
            z=obj.zk_s;
            var_z=obj.diag_Pk_s; 
            basis = obj.stem_fda.spline_basis_z;
            %p=50:0.25:1000;
            b=eval_basis(f_range,basis);
            %mapshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
            mapshow(obj.shape);
            hold on
            
            plot(obj.stem_grid_sites.coordinate(:,2),obj.stem_grid_sites.coordinate(:,1),'.r');
            lat=reshape(obj.stem_grid.coordinate(:,1),obj.stem_grid.grid_size(1),obj.stem_grid.grid_size(2));
            lon=reshape(obj.stem_grid.coordinate(:,2),obj.stem_grid.grid_size(1),obj.stem_grid.grid_size(2));
            xlim([min(lon(:)),max(lon(:))]);
            ylim([min(lat(:)),max(lat(:))]);
            
            grid_size=obj.stem_grid.grid_size;
            [min_value,idx]=min(abs(obj.stem_grid.coordinate(:,1)-lat0)+abs(obj.stem_grid.coordinate(:,2)-lon0));
            [lon_idx,lat_idx]=ind2sub(grid_size,idx);
            
            %t=1;
            profile=b*(squeeze(z(lon_idx,lat_idx,t,:)).*alpha_z);
            v_coeff=diag(alpha_z_2.*squeeze(var_z(lon_idx,lat_idx,t,:)));
            
            if obj.stem_fda.flag_beta_spline==1
            end
            profile1 = profile+X_beta*obj.stem_par.beta;
            y_hat = profile1;
            
            v_spline=b*v_coeff*b';
            std_spline=sqrt(diag(v_spline));
            diag_Var_y_hat = std_spline;
            
            figure;
            subplot(1,2,1);
            plot(profile1,f_range);
            hold on
            plot(profile1-2*std_spline,f_range,'--r');
            plot(profile1+2*std_spline,f_range,'--r');
            xlabel('K');
            ylabel('hPa');
            set(gca,'Ydir','reverse')
            %ylim([50,1000]);
            title(['lat=',num2str(lat0),'? lon=',num2str(lon0),'? t=',num2str(t)]);

            subplot(1,2,2);
            imagesc(f_range,f_range,v_spline);
            title(['f-varcov lat=',num2str(lat0),'? lon=',num2str(lon0),'? t=',num2str(t)]);
            xlabel('hPa');
            ylabel('hPa');
            axis equal
            axis tight
            colorbar
            
        end
            
    end
end