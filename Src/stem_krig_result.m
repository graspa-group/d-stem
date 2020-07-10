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
%%% Author: Alessandro Fassò                                             %
%%% E-mail: alessandro.fasso@unibg.it                                    %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Information and                    %
%%%              Production Engineering                                  %
%%% Author website: http://www.unibg.it/pers/?alessandro.fasso           %
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

classdef stem_krig_result < handle
    
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
    %NN - number of kriging sites
    %T - number of temporal steps    
    %C - number of basis functions (when the model is f-HDGM)
    %K - the number of loading vectors related to the latent variable w_p
    
    properties
        variable_name=[];       %[string]                (1x1) the name of the kriged variable
        Y_unit={};              %[string]                (1x1) unit of variable y

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
    end

    properties (Hidden = true)
        stem_fda=[];                 %[stem_fda object]       (1x1) the stem_fda object of an estimated model
        stem_par=[];                 %[stem_par object]       (1x1) the stem_par object of an estimated model
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
            %shape              - [struct]                (1x1) boundary of the geographic region loaded from a shape file
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
                    title([obj.variable_name,' - ',datestr(obj.stem_datestamp.stamp(time_step))],'FontSize',14);
                else
                    temp=mean(obj.y_hat,3);
                    title(['Average ',obj.variable_name,' from ',datestr(obj.stem_datestamp.stamp(1)),' to ',datestr(obj.stem_datestamp.stamp(end))],'FontSize',14);
                end
                h = mapshow(lon,lat,temp,'DisplayType','texturemap');
                set(h,'FaceColor','flat');
                if not(isempty(obj.shape))
                    geoshow(obj.shape);
                end
                axis equal
                xlim([min(lon(:)),max(lon(:))]);
                ylim([min(lat(:)),max(lat(:))]);
                if strcmp(obj.stem_grid.unit,'deg')
                    xlabel('Longitude','FontSize',14);
                    ylabel('Latitude','FontSize',14);
                else
                    xlabel(obj.stem_grid.unit);
                    ylabel(obj.stem_grid.unit);
                end
                geoshow(obj.stem_grid_sites.coordinate(:,1),obj.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','w');
                
                colorbar;
                grid on;
                box on;
                set(gca,'FontSize',14);
                set(gcf, 'renderer', 'zbuffer');
            end

            if strcmp(type,'both')
                subplot(2,1,2);
            end
            if strcmp(type,'std')||strcmp(type,'both')
                hold on
                if time_step>0
                    temp=sqrt(obj.diag_Var_y_hat(:,:,time_step));
                    title(['Std of ',obj.variable_name,' - ',datestr(obj.stem_datestamp.stamp(time_step))],'FontSize',14);
                else
                    temp=mean(sqrt(obj.diag_Var_y_hat),3);
                    title(['Average std of ',obj.variable_name,' from ',datestr(obj.stem_datestamp.stamp(1)),' to ',datestr(obj.stem_datestamp.stamp(end))],'FontSize',14);
                end
                h = mapshow(lon,lat,temp,'DisplayType','texturemap');
                set(h,'FaceColor','flat');
                if not(isempty(obj.shape))
                    geoshow(obj.shape);
                end
                axis equal
                xlim([min(lon(:)),max(lon(:))]);
                ylim([min(lat(:)),max(lat(:))]);
                if strcmp(obj.stem_grid.unit,'deg')
                    xlabel('Longitude','FontSize',14);
                    ylabel('Latitude','FontSize',14);
                else
                    xlabel(obj.stem_grid.unit,'FontSize',14);
                    ylabel(obj.stem_grid.unit,'FontSize',14);
                end    
                geoshow(obj.stem_grid_sites.coordinate(:,1),obj.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','w');
                
                colorbar;
                grid on;
                box on;
                set(gca,'FontSize',14);
                set(gcf, 'renderer', 'zbuffer');
            end
            if nargout>0
                fig_h=h;
            end
        end
        
        function [y_hat,diag_Var_y_hat] = surface_plot(obj,h,t,X_beta)
            %DESCRIPTION: surface plot of the kriged variable and its standard deviation 
            %
            %INPUT
            %obj                - [stem_krig_result object]  (1x1) the stem_krig_result object 
            %h                  - [double]                   (1x1|1x2) the value of h or the domain interval
            %t                  - [integer >=0]              (1x1|1x2) the date or the date interval
            %X_beta             - [double]                   (NxMxb)   the value of X_beta at h and time t
            %
            %OUTPUT
            %y_hat              - [integer]                  (1x1) the kriged variable
            %diag_Var_y_hat     - [integer]                  (1x1) the standard deviation of the kriged variable

            if nargin<3
                error('Not enough input arguments');
            end
            
            if isempty(obj.stem_fda)
                error('This function is only available for the model type f-HDGM ')
            end
            
            if not(length(h)==1)
                error('h must be a scalar');
            end
            
            if not(length(t)==1)
                error('t must be a scalar');
            end
            
            if abs(obj.stem_grid.box(4)-obj.stem_grid.box(3))>180
                global_idx=1;
            else
                global_idx=0;
            end

            
            if nargin<4
                
                warning('The surface plot is about the \phi(h)z(s,t).')

                z=obj.zk_s;
                var_z=obj.diag_Pk_s; 
                
                basis = obj.stem_fda.spline_basis_z;
                range_z = getbasisrange(basis);
                if h<range_z(1)||h>range_z(end)
                    error('The argument h is out of the basis range.');
                end
                
                b=full(eval_basis(h,basis));
                s=nan(size(z,1),size(z,2));
                v=nan(size(z,1),size(z,2));
               
                for i=1:size(z,1)
                    s(i,:)=squeeze(z(i,:,t,:))*b';
                    v(i,:)=squeeze(var_z(i,:,t,:))*(b.^2)';
                end
                surface1=s;

                if nargout==1
                    y_hat = surface1;
                end

                if nargout==2
                    y_hat = surface1;
                    diag_Var_y_hat = v;
                end
                  
                lat=reshape(obj.stem_grid.coordinate(:,1),obj.stem_grid.grid_size(1),obj.stem_grid.grid_size(2));
                lon=reshape(obj.stem_grid.coordinate(:,2),obj.stem_grid.grid_size(1),obj.stem_grid.grid_size(2));

                if global_idx

                    f1=figure;
                    axesm eckert4; 
                    framem; gridm;
                    axis off

                    title(['\phi(h)''z(s,t) on ',datestr(obj.stem_datestamp.stamp(t)),' @ h=',num2str(h)],'FontSize',14);

                    hh = geoshow(lat,lon,surface1,'DisplayType','texturemap');
                    set(hh,'FaceColor','flat');
                    if not(isempty(obj.shape))
                        try
                            geoshow(obj.shape,'FaceColor','none');
                        catch
                            geoshow(obj.shape);
                        end
                    end

                    if strcmp(obj.stem_grid.unit,'deg')
                        xlabel('Longitude [deg]','FontSize',14);
                        ylabel('Latitude [deg]','FontSize',14);
                    else
                        xlabel(obj.stem_grid.unit);
                        ylabel(obj.stem_grid.unit);
                    end

                    geoshow(obj.stem_grid_sites.coordinate(:,1),obj.stem_grid_sites.coordinate(:,2),...,
                        'DisplayType','multipoint','Marker','*','MarkerEdgeColor','k');

                    
                    if min(surface1(:))*max(surface1(:))>0
                        colormap(f1,flipud(stem_misc.get_s_colormap()))
                        cl = colorbar;
                        cl.Limits=[min(surface1(:)) max(surface1(:))]; 
                    else
                        colormap(f1,stem_misc.get_d_colormap())
                        cl = colorbar;
                        caxis( [-max(abs(surface1(:))) max(abs(surface1(:)))] ); 
                    end
                    ylabel(cl,['[',obj.Y_unit,']']);
                    
                    grid on;
                    box on;
                    set(gca,'FontSize',14);
                    set(gcf, 'renderer', 'zbuffer');

                    f2=figure;
                    axesm eckert4; 
                    framem; gridm;
                    axis off

                    temp=sqrt(v);
                    title(['Std of \phi(h)''z(s,t) on ',datestr(obj.stem_datestamp.stamp(t)),' @ h=',num2str(h)],'FontSize',14);

                    hh = geoshow(lat,lon,temp,'DisplayType','texturemap');
                    set(hh,'FaceColor','flat');
                    if not(isempty(obj.shape))
                        try
                            geoshow(obj.shape,'FaceColor','none');
                        catch
                            geoshow(obj.shape);
                        end
                    end

                    if strcmp(obj.stem_grid.unit,'deg')
                        xlabel('Longitude [deg]','FontSize',14);
                        ylabel('Latitude [deg]','FontSize',14);
                    else
                        xlabel(obj.stem_grid.unit,'FontSize',14);
                        ylabel(obj.stem_grid.unit,'FontSize',14);
                    end    
                    geoshow(obj.stem_grid_sites.coordinate(:,1),obj.stem_grid_sites.coordinate(:,2),...
                        'DisplayType','multipoint','Marker','*','MarkerEdgeColor','k');

                    if min(temp(:))*max(temp(:))>0
                        colormap(f2,flipud(stem_misc.get_s_colormap()))
                        cl = colorbar;
                        cl.Limits=[min(temp(:)) max(temp(:))]; 
                    else
                        colormap(f2,stem_misc.get_d_colormap())
                        cl = colorbar;
                        caxis( [-max(abs(temp(:))) max(abs(temp(:)))] ); 
                    end
                    ylabel(cl,['[',obj.Y_unit,']']);
                    
                    grid on;
                    box on;
                    set(gca,'FontSize',14);
                    set(gcf, 'renderer', 'zbuffer');
                else

                    figure;
                    ax1=subplot(1,2,1); 

                    hold on
                    title(['\phi(h)''z(s,t) on ',datestr(obj.stem_datestamp.stamp(t)),' @ h=',num2str(h)],'FontSize',14);

                    hh = geoshow(lat,lon,surface1,'DisplayType','texturemap');
                    set(hh,'FaceColor','flat');
                    if not(isempty(obj.shape))
                        try
                            geoshow(obj.shape,'FaceColor','none');
                        catch
                            geoshow(obj.shape);
                        end
                    end

                    axis equal
                    xlim([min(lon(:)),max(lon(:))]);
                    ylim([min(lat(:)),max(lat(:))]);

                    if strcmp(obj.stem_grid.unit,'deg')
                        xlabel('Longitude [deg]','FontSize',14);
                        ylabel('Latitude [deg]','FontSize',14);
                    else
                        xlabel(obj.stem_grid.unit);
                        ylabel(obj.stem_grid.unit);
                    end

                    geoshow(obj.stem_grid_sites.coordinate(:,1),obj.stem_grid_sites.coordinate(:,2),...
                        'DisplayType','multipoint','Marker','*','MarkerEdgeColor','k');

                    if min(surface1(:))*max(surface1(:))>0
                        colormap(ax1,flipud(stem_misc.get_s_colormap()))
                        cl = colorbar;
                        cl.Limits=[min(surface1(:)) max(surface1(:))]; 
                    else
                        colormap(ax1,stem_misc.get_d_colormap())
                        cl = colorbar;
                        caxis( [-max(abs(surface1(:))) max(abs(surface1(:)))] );
                    end
                    ylabel(cl,['[',obj.Y_unit,']']);
                    
                    grid on;
                    box on;
                    set(gca,'FontSize',14);
                    set(gcf, 'renderer', 'zbuffer');

                    ax2=subplot(1,2,2); 
                    hold on
                    temp=sqrt(v);
                    title(['Std of \phi(h)''z(s,t) on ',datestr(obj.stem_datestamp.stamp(t)),' @ h=',num2str(h)],'FontSize',14);

                    hh = geoshow(lat,lon,temp,'DisplayType','texturemap');
                    set(hh,'FaceColor','flat');
                    if not(isempty(obj.shape))
                        try
                            geoshow(obj.shape,'FaceColor','none');
                        catch
                            geoshow(obj.shape);
                        end
                    end

                    axis equal
                    xlim([min(lon(:)),max(lon(:))]);
                    ylim([min(lat(:)),max(lat(:))]);

                    if strcmp(obj.stem_grid.unit,'deg')
                        xlabel('Longitude [deg]','FontSize',14);
                        ylabel('Latitude [deg]','FontSize',14);
                    else
                        xlabel(obj.stem_grid.unit,'FontSize',14);
                        ylabel(obj.stem_grid.unit,'FontSize',14);
                    end    
                    geoshow(obj.stem_grid_sites.coordinate(:,1),obj.stem_grid_sites.coordinate(:,2),...
                        'DisplayType','multipoint','Marker','*','MarkerEdgeColor','k');

                    if min(temp(:))*max(temp(:))>0
                        colormap(ax2,flipud(stem_misc.get_s_colormap()))
                        cl = colorbar;
                        cl.Limits=[min(temp(:)) max(temp(:))]; 
                    else
                        colormap(ax2,stem_misc.get_d_colormap())
                        cl = colorbar;
                        caxis( [-max(abs(temp(:))) max(abs(temp(:)))] ); 
                    end
                    ylabel(cl,['[',obj.Y_unit,']']);
                    
                    grid on;
                    box on;
                    set(gca,'FontSize',14);
                    set(gcf, 'renderer', 'zbuffer');
                end
            else
                if obj.stem_fda.flag_beta_spline==1
                    k=obj.stem_fda.spline_nbasis_beta;
                    if size(X_beta,3)~=(length(obj.stem_par.beta)/k)
                        warning('You must provide the p dimension covariates X_beta ')
                    end
                else
                    if size(X_beta,3)~=length(obj.stem_par.beta)
                        warning('You must provide the p dimension covariates X_beta ')
                    end
                end

                z=obj.zk_s;
                var_z=obj.diag_Pk_s; 

                if obj.stem_fda.flag_beta_spline==1
                    b = size(X_beta,3);
                    k=getnbasis(obj.stem_fda.spline_basis_beta);
                    X_beta_tmp = nan(size(z,1),size(z,2),b*k);
                    basis_beta = obj.stem_fda.spline_basis_beta;
                    b_beta=full(eval_basis(h,basis_beta));
                    for lat_i=1:obj.stem_grid.grid_size(1)
                        for lon_i = 1:obj.stem_grid.grid_size(2)
                            for i = 1:b
                                X_beta_tmp(lat_i,lon_i,(i-1)*k+(1:k)) = squeeze(X_beta(lat_i,lon_i,i)).*b_beta;
                            end
                        end
                    end 
                    X_beta = X_beta_tmp;
                end

                basis = obj.stem_fda.spline_basis_z;
                b=full(eval_basis(h,basis));
                s=nan(size(z,1),size(z,2));
                v=nan(size(z,1),size(z,2));
                Xbeta = nan(size(z,1),size(z,2));
                surface1 = nan(size(z,1),size(z,2));
                for i=1:size(z,1)
                    s(i,:)=squeeze(z(i,:,t,:))*b';
                    v(i,:)=squeeze(var_z(i,:,t,:))*(b.^2)';
                    Xbeta(i,:)=squeeze(X_beta(i,:,:))*obj.stem_par.beta;  
                    surface1(i,:)=Xbeta(i,:)+s(i,:);
                end

                if nargout==1
                    y_hat = surface1;
                end

                if nargout==2
                    y_hat = surface1;
                    diag_Var_y_hat = v;
                end


                lat=reshape(obj.stem_grid.coordinate(:,1),obj.stem_grid.grid_size(1),obj.stem_grid.grid_size(2));
                lon=reshape(obj.stem_grid.coordinate(:,2),obj.stem_grid.grid_size(1),obj.stem_grid.grid_size(2));

                if global_idx

                    f1=figure;
                    axesm eckert4; 
                    framem; gridm;
                    axis off

                    temp1=strsplit(obj.variable_name,'_');
                    title([temp1{1},' on ',datestr(obj.stem_datestamp.stamp(t)),' @ h=',num2str(h)],'FontSize',14);

                    hh = geoshow(lat,lon,surface1,'DisplayType','texturemap');
                    set(hh,'FaceColor','flat');
                    if not(isempty(obj.shape))
                        try
                            geoshow(obj.shape,'FaceColor','none');
                        catch
                            geoshow(obj.shape);
                        end
                    end

                    if strcmp(obj.stem_grid.unit,'deg')
                        xlabel('Longitude [deg]','FontSize',14);
                        ylabel('Latitude [deg]','FontSize',14);
                    else
                        xlabel(obj.stem_grid.unit);
                        ylabel(obj.stem_grid.unit);
                    end

                    geoshow(obj.stem_grid_sites.coordinate(:,1),obj.stem_grid_sites.coordinate(:,2),'DisplayType','multipoint','Marker','*','MarkerEdgeColor','k');

                    if min(surface1(:))*max(surface1(:))>0
                        colormap(f1,flipud(stem_misc.get_s_colormap()))
                        cl = colorbar;
                        cl.Limits=[min(surface1(:)) max(surface1(:))]; 
                    else
                        colormap(f1,stem_misc.get_d_colormap())
                        colorbar;
                        caxis( [-max(abs(surface1(:))) max(abs(surface1(:)))] ); 
                    end
                    
                    grid on;
                    box on;
                    set(gca,'FontSize',14);
                    set(gcf, 'renderer', 'zbuffer');


                    f2=figure;
                    axesm eckert4; 
                    framem; gridm;
                    axis off

                    temp=sqrt(v);
                    temp1=strsplit(obj.variable_name,'_');
                    title(['Std of ',temp1{1},' on ',datestr(obj.stem_datestamp.stamp(t)),' @ h=',num2str(h)],'FontSize',14);

                    hh = geoshow(lat,lon,temp,'DisplayType','texturemap');
                    set(hh,'FaceColor','flat');
                    if not(isempty(obj.shape))
                        try
                            geoshow(obj.shape,'FaceColor','none');
                        catch
                            geoshow(obj.shape);
                        end
                    end

                    if strcmp(obj.stem_grid.unit,'deg')
                        xlabel('Longitude [deg]','FontSize',14);
                        ylabel('Latitude [deg]','FontSize',14);
                    else
                        xlabel(obj.stem_grid.unit,'FontSize',14);
                        ylabel(obj.stem_grid.unit,'FontSize',14);
                    end    
                    geoshow(obj.stem_grid_sites.coordinate(:,1),obj.stem_grid_sites.coordinate(:,2),'DisplayType','multipoint','Marker','*','MarkerEdgeColor','k');

                    if min(temp(:))*max(temp(:))>0
                        colormap(f2,flipud(stem_misc.get_s_colormap()))
                        cl = colorbar;
                        cl.Limits=[min(temp(:)) max(temp(:))]; 
                    else
                        colormap(f2,stem_misc.get_d_colormap())
                        colorbar;
                        caxis( [-max(abs(temp(:))) max(abs(temp(:)))] ); 
                    end
                    
                    grid on;
                    box on;
                    set(gca,'FontSize',14);
                    set(gcf, 'renderer', 'zbuffer');
                else
                    
                    figure;
                    ax1=subplot(1,2,1); 
                    hold on

                    temp1=strsplit(obj.variable_name,'_');
                    title([temp1{1},' on ',datestr(obj.stem_datestamp.stamp(t)),...
                        ' @ h=',num2str(h)],'FontSize',14);

                    hh = geoshow(lat,lon,surface1,'DisplayType','texturemap');
                    set(hh,'FaceColor','flat');
                    if not(isempty(obj.shape))
                        try
                            geoshow(obj.shape,'FaceColor','none');
                        catch
                            geoshow(obj.shape);
                        end
                    end

                    axis equal
                    xlim([min(lon(:)),max(lon(:))]);
                    ylim([min(lat(:)),max(lat(:))]);

                    if strcmp(obj.stem_grid.unit,'deg')
                        xlabel('Longitude [deg]','FontSize',14);
                        ylabel('Latitude [deg]','FontSize',14);
                    else
                        xlabel(obj.stem_grid.unit);
                        ylabel(obj.stem_grid.unit);
                    end

                    geoshow(obj.stem_grid_sites.coordinate(:,1),obj.stem_grid_sites.coordinate(:,2),...
                        'DisplayType','multipoint','Marker','*','MarkerEdgeColor','k');

                    if min(surface1(:))*max(surface1(:))>0
                        colormap(ax1,flipud(stem_misc.get_s_colormap()))
                        cl = colorbar;
                        cl.Limits=[min(surface1(:)) max(surface1(:))]; 
                    else
                        colormap(ax1,stem_misc.get_d_colormap())
                        colorbar;
                        caxis( [-max(abs(surface1(:))) max(abs(surface1(:)))] ); 
                    end
                    
                    grid on;
                    box on;
                    set(gca,'FontSize',14);
                    set(gcf, 'renderer', 'zbuffer');

                    ax2=subplot(1,2,2);
                    hold on

                    temp=sqrt(v);
                    temp1=strsplit(obj.variable_name,'_');
                    title(['Std of ',temp1{1},' on ',datestr(obj.stem_datestamp.stamp(t)),...
                        ' @ h=',num2str(h)],'FontSize',14);

                    hh = geoshow(lat,lon,temp,'DisplayType','texturemap');
                    set(hh,'FaceColor','flat');
                    if not(isempty(obj.shape))
                        try
                            geoshow(obj.shape,'FaceColor','none');
                        catch
                            geoshow(obj.shape);
                        end
                    end

                    axis equal
                    xlim([min(lon(:)),max(lon(:))]);
                    ylim([min(lat(:)),max(lat(:))]);

                    if strcmp(obj.stem_grid.unit,'deg')
                        xlabel('Longitude [deg]','FontSize',14);
                        ylabel('Latitude [deg]','FontSize',14);
                    else
                        xlabel(obj.stem_grid.unit,'FontSize',14);
                        ylabel(obj.stem_grid.unit,'FontSize',14);
                    end    
                    geoshow(obj.stem_grid_sites.coordinate(:,1),obj.stem_grid_sites.coordinate(:,2),...
                        'DisplayType','multipoint','Marker','*','MarkerEdgeColor','k');

                    if min(temp(:))*max(temp(:))>0
                        colormap(ax2,flipud(stem_misc.get_s_colormap()))
                        cl = colorbar;
                        cl.Limits=[min(temp(:)) max(temp(:))]; 
                    else
                        colormap(ax2,stem_misc.get_d_colormap())
                        colorbar;
                        caxis( [-max(abs(temp(:))) max(abs(temp(:)))] ); 
                    end
                    
                    grid on;
                    box on;
                    set(gca,'FontSize',14);
                    set(gcf, 'renderer', 'zbuffer');
                end
         
            end
        
        end
        
        function [y_hat,diag_Var_y_hat] = profile_plot(obj,h,lon,lat,t,X_beta,vertical)
            %DESCRIPTION: plot the profile of kriging result at target site 
            %
            %INPUT
            %obj                 - [stem_krig_result object]  (1x1) the stem_krig_result object 
            %h_range             - [double]                   (lx1) the range of domain of functional data to predict the the kriged variable
            %lon                 - [double]                   (1x1) the lontitude of the spatial site
            %lat                 - [double]                   (1x1) the latitude of the spatial site
            %X_beta              - [double]                   (lxb) the value of X_beta at day t, where l is the length of h vector, and b is the number of covariates
            %vertical            - [boolean]                  (1x1) indicating if the plot is vertical or not.
            %
            %OUTPUT
            %y_hat               - [integer]                  (1x1) the kriged variable
            %diag_Var_y_hat      - [integer]                  (1x1) the standard deviation of the kriged variable
            % 
            %Profile plot of the kriged variable and its standard deviation
            
            if nargin<6
                X_beta=[];
                vertical=0;
            end
            
            if nargin==6
                vertical=0;
                %error('Not enough input arguments');
            end
            
            if isempty(X_beta)
                
                warning('The surface plot is about the latent variable Z.')

                grid_size=obj.stem_grid.grid_size;
                [~,idx]=min(abs(obj.stem_grid.coordinate(:,1)-lat)+abs(obj.stem_grid.coordinate(:,2)-lon));
                [lon_idx,lat_idx]=ind2sub(grid_size,idx);

                z=obj.zk_s;
                var_z=obj.diag_Pk_s; 
                basis_z = obj.stem_fda.spline_basis_z;
                b_range =getbasisrange(basis_z);
                if b_range(1)<h(1)||b_range(2)>h(2)
                    warning('The begining/end of h is not covered')
                end
                basis_range = b_range(1):0.1:b_range(2);
                b_z=eval_basis(basis_range,basis_z);

                profile=b_z*squeeze(z(lon_idx,lat_idx,t,:));
                v_coeff=diag(squeeze(var_z(lon_idx,lat_idx,t,:)));

               
                v_spline=b_z*v_coeff*b_z';
                std_spline=sqrt(diag(v_spline));

                if nargout==1
                    y_hat = profile;
                end

                if nargout==2
                    y_hat = profile;
                    diag_Var_y_hat = std_spline;
                end
                quant = [0.005, 0.025, 0.05];
                
                
                figure;
                ax1=subplot(1,2,1);
                if vertical
                    plot(profile,basis_range);
                    for j=1:3
                        coef = -norminv(quant(j),0,1);
                        profile1_up = profile+coef*std_spline;
                        profile1_low = profile-coef*std_spline;  
                        patch([profile1_up' fliplr(profile1_low')],[basis_range fliplr(basis_range)],'r','FaceAlpha',0.15,'EdgeAlpha',0)
                    end
                    ylim([b_range(1), b_range(2)])

                    xlabel('\phi(h)''z(s,t)','FontSize',14);
                    ylabel('h','FontSize',14);
                    set(gca,'Ydir','reverse','FontSize',14);
                    ylim([b_range(1), b_range(2)]);
                    axis square

                    title(['Lat ',num2str(lat),' Lon ',num2str(lon),' on ',datestr(obj.stem_datestamp.stamp(t))],'FontSize',14);

                    ax2=subplot(1,2,2);
                    imagesc(basis_range,basis_range,v_spline);
                    title(['h-varcov(\phi(h)''z(s,t)) Lat ',num2str(lat),' Lon ',num2str(lon),' on ',datestr(obj.stem_datestamp.stamp(t))],'FontSize',14);
                    xlabel('h','FontSize',14);
                    ylabel('h','FontSize',14);
                    set(gca,'Xdir','reverse','FontSize',14);
                    colormap(ax2,flipud(stem_misc.get_s_colormap()))
                    colorbar
                    axis square
                    ax2.Position(2:4)=ax1.Position(2:4);
                    ax2.PlotBoxAspectRatio=ax1.PlotBoxAspectRatio;
                else
                    plot(basis_range,profile,'k');
                    for j=1:3
                        coef = -norminv(quant(j),0,1);
                        profile1_up = profile+coef*std_spline;
                        profile1_low = profile-coef*std_spline;  
                        patch([basis_range fliplr(basis_range)],[profile1_up' fliplr(profile1_low')],'r','FaceAlpha',0.15,'EdgeAlpha',0)
                    end
                    xlim([b_range(1), b_range(2)])
                    ylabel('\phi(h)''z(s,t)','FontSize',14);
                    xlabel('h','FontSize',14);
                    axis square

                    title(['Lat ',num2str(lat),' Lon ',num2str(lon),' on ',datestr(obj.stem_datestamp.stamp(t))],'FontSize',14);

                    ax2=subplot(1,2,2);
                    imagesc(basis_range,basis_range,v_spline);
                    title(['h-varcov(\phi(h)''z(s,t)) Lat ',num2str(lat),' Lon ',num2str(lon),' on ',datestr(obj.stem_datestamp.stamp(t))],'FontSize',14);
                    xlabel('h','FontSize',14);
                    ylabel('h','FontSize',14);
                    set(gca,'Ydir','normal','FontSize',14);
                    colormap(ax2,flipud(stem_misc.get_s_colormap()))
                    colorbar
                    axis square
                    ax2.Position(2:4)=ax1.Position(2:4);
                    ax2.PlotBoxAspectRatio=ax1.PlotBoxAspectRatio;
                end
                
            else
                
                basis_beta = obj.stem_fda.spline_basis_beta;
                b_range =getbasisrange(basis_beta);
                if b_range(1)<h(1)||b_range(2)>h(end)
                    warning(['Basis functions have range [',num2str(b_range(1)),',',num2str(b_range(2)),'] while covariates provided to profile_plot have range [',num2str(h(1)),',',num2str(h(end)),']. Extrapolation is performed']);
                end
                basis_range = b_range(1):0.1:b_range(2);

                grid_size=obj.stem_grid.grid_size;
                [~,idx]=min(abs(obj.stem_grid.coordinate(:,1)-lat)+abs(obj.stem_grid.coordinate(:,2)-lon));
                [lon_idx,lat_idx]=ind2sub(grid_size,idx);

             
                z=obj.zk_s;
                var_z=obj.diag_Pk_s; 
                basis_z = obj.stem_fda.spline_basis_z;
                b_z=eval_basis(basis_range,basis_z);

                profile=b_z*squeeze(z(lon_idx,lat_idx,t,:));
                v_coeff=diag(squeeze(var_z(lon_idx,lat_idx,t,:)));


                X_beta_inter = interp1(h,X_beta,basis_range,'linear','extrap');

                if obj.stem_fda.flag_beta_spline==1
                    b = size(X_beta,2);
                    k=getnbasis(obj.stem_fda.spline_basis_beta);
                    X_beta_tmp = nan(length(basis_range),b*k);
                    b_beta=full(eval_basis(basis_range,basis_beta));
                    for i = 1:b
                        X_beta_tmp(:,(i-1)*k+(1:k)) = X_beta_inter(:,i).*b_beta;
                    end
                    X_beta = X_beta_tmp;
                end

                profile1 = profile+X_beta*obj.stem_par.beta;

                v_spline=b_z*v_coeff*b_z';
                std_spline=sqrt(diag(v_spline));

                if nargout==1
                    y_hat = profile1;
                end

                if nargout==2
                    y_hat = profile1;
                    diag_Var_y_hat = std_spline;
                end
                quant = [0.005, 0.025, 0.05];
                
                figure;
                ax1=subplot(1,2,1);
                if vertical
                    plot(profile1,basis_range);
                    for j=1:3
                        coef = -norminv(quant(j),0,1);
                        profile1_up = profile1+coef*std_spline;
                        profile1_low = profile1-coef*std_spline;  
                        patch([profile1_up' fliplr(profile1_low')],[basis_range fliplr(basis_range)],...
                            'r','FaceAlpha',0.15,'EdgeAlpha',0)
                    end
                    ylim([b_range(1), b_range(2)]);
                    xlabel(obj.variable_name,'FontSize',14);
                    ylabel('h','FontSize',14);
                    set(gca,'Ydir','reverse','FontSize',14); 
                    axis square
                    title(['Lat ',num2str(lat),' Lon ',num2str(lon),' on ',...
                        datestr(obj.stem_datestamp.stamp(t))],'FontSize',14);

                    ax2=subplot(1,2,2);
                    imagesc(basis_range,basis_range,v_spline);
                    title(['h-varcov(',obj.variable_name,') Lat ',num2str(lat),' Lon ',num2str(lon),' on ',...
                        datestr(obj.stem_datestamp.stamp(t))],'FontSize',14);
                    xlabel('h','FontSize',14);
                    ylabel('h','FontSize',14);
                    set(gca,'Xdir','reverse','FontSize',14);
                    axis square 
                   
                    colormap(ax2,flipud(stem_misc.get_s_colormap()))
                    colorbar
                    ax2.Position(2:4)=ax1.Position(2:4);
                    ax2.PlotBoxAspectRatio=ax1.PlotBoxAspectRatio;
                else
                    plot(basis_range,profile1,'k');
                    for j=1:3
                        coef = -norminv(quant(j),0,1);
                        profile1_up = profile1+coef*std_spline;
                        profile1_low = profile1-coef*std_spline;  
                        patch([basis_range fliplr(basis_range)],[profile1_up' fliplr(profile1_low')],'r','FaceAlpha',0.15,'EdgeAlpha',0)
                    end
                    xlim([b_range(1), b_range(2)]);
                    ylabel(obj.variable_name,'FontSize',14);
                    xlabel('h','FontSize',14);
                    set(gca,'FontSize',14); 
                    axis square
                    title(['Lat ',num2str(lat),' Lon ',num2str(lon),' on ',...
                        datestr(obj.stem_datestamp.stamp(t))],'FontSize',14);

                    ax2=subplot(1,2,2);
                    imagesc(basis_range,basis_range,v_spline);
                    title(['h-varcov(',obj.variable_name,') Lat ',num2str(lat),' Lon ',num2str(lon),' on ',...
                        datestr(obj.stem_datestamp.stamp(t))],'FontSize',14);
                    xlabel('h','FontSize',14);
                    ylabel('h','FontSize',14);
                    set(gca,'Ydir','normal','FontSize',14);
                    axis square
                    colormap(ax2,flipud(stem_misc.get_s_colormap()))
                    colorbar
                    ax2.Position(2:4)=ax1.Position(2:4);
                    ax2.PlotBoxAspectRatio=ax1.PlotBoxAspectRatio;
                end
            end
        end 
    end
end