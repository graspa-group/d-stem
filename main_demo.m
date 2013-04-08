%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D-STEM - Distributed Space Time Expecation Maximization      %
%                                                              %
% Author: Francesco Finazzi                                    %
% E-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo - Dept. of Engineering    %
% Author website: http://www.unibg.it/pers/?francesco.finazzi  %
% Code website: https://code.google.com/p/d-stem/              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all


flag_bivariate=1;       %1: NO2 and PM2.5 data; 0: only NO2 data
flag_remote_data=1;     %1: remote sensing data are considered; 0: only ground level data are considered
flag_remote_corr=0;     %1: remote data are cross-correlated (if they are multivariate); 0: remote data are not cross-correlated

flag_beta_ground=1;     %1: the model includes the term X_beta(s,t)*beta; 0: X_beta(s,t)*beta is not included
flag_time_ground=1;     %1: the model includes the term X_z(s,t)*z(t); 0: X_z(s,t)*z(t) is not included
flag_w_ground=1;        %1: the model includes the term X_g(s,t)*w(s,t);  0: X_g(s,t)*w(s,t) is not included

flag_crossval=0;        %1: cross-validation is enabled; 0: cross-validation is not enabled
flag_tapering=0;        %1: tapering is enabled; 0: tapering is not enabled 

flag_estimate=1;        %1: model estimation is enabled; 0:model estimation is not enabled
flag_kriging=1;         %1: kriging is enabled; 0: kriging is not enabled


if flag_estimate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      Data  building     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %load the no2 ground level observations
    load ../Demo/Data/no2_ground/no2_ground_background.mat
    sd_g.Y{1}=no2_ground.data;
    sd_g.Y_name{1}='no2 ground';
    n1=size(sd_g.Y{1},1);
    T=size(sd_g.Y{1},2);
    
    if flag_bivariate
        %load the pm25 ground level observations
        load ../Demo/Data/pm25_ground/pm25_ground_background.mat
        sd_g.Y{2}=pm25_ground.data;
        sd_g.Y_name{2}='pm2.5 ground';
        n2=size(sd_g.Y{2},1);
    end
    
    %X_rg
    if flag_remote_data
        sd_g.X_rg{1}=ones(n1,1);
        sd_g.X_rg_name{1}={'constant'};
        if flag_bivariate
            sd_g.X_rg{2}=ones(n2,1);
            sd_g.X_rg_name{2}={'constant'};
        end
    else
        sd_g.X_rg=[];
        sd_g.X_rg_name{1}=[];
    end

    if flag_beta_ground
        %X_beta
        %load the covariates over all the NO2 monitoring stations
        load ../Demo/Data/no2_ground/no2_background_meteo_point.mat
        load ../Demo/Data/no2_ground/no2_background_elevation_point.mat
        load ../Demo/Data/no2_ground/no2_background_emission_point.mat
        load ../Demo/Data/no2_ground/no2_background_population_point.mat
        X_elevation_point=repmat(X_elevation_point,[1,1,T]);
        X_emission_point=repmat(X_emission_point,[1,1,T]);
        X_population_point(X_population_point==0)=0.02;
        X_population_point=log(X_population_point);
        X_population_point=repmat(X_population_point,[1,1,T]);
        %dummy variables for saturday and sunday
        X_dummy1=zeros(size(X_population_point,1),1,365);
        X_dummy2=zeros(size(X_population_point,1),1,365);
        X_dummy1(:,1,3:7:365)=1;
        X_dummy2(:,1,4:7:365)=1;
        X=cat(2,X_meteo_point,X_elevation_point,X_emission_point,X_population_point,X_dummy1,X_dummy2);
        X=double(X);        
        sd_g.X_beta{1}=X;
        sd_g.X_beta_name{1}={'wind speed','pressure','temperature','elevation','emission','population','saturday','sunday'};
        
        clear X_meteo_point
        clear X_emission_point
        clear X_elevation_point
        clear X_population_point
        
        if flag_bivariate
            load ../Demo/Data/pm25_ground/pm25_background_meteo_point.mat
            load ../Demo/Data/pm25_ground/pm25_background_elevation_point.mat
            load ../Demo/Data/pm25_ground/pm10_background_emission_point.mat
            load ../Demo/Data/pm25_ground/pm25_background_population_point.mat
            X_elevation_point=repmat(X_elevation_point,[1,1,T]);
            X_emission_point=repmat(X_emission_point,[1,1,T]);
            X_population_point(X_population_point==0)=0.02;
            X_population_point=log(X_population_point);
            X_population_point=repmat(X_population_point,[1,1,T]);
            %dummy variables for saturday and sunday
            X_dummy1=zeros(size(X_population_point,1),1,365);
            X_dummy2=zeros(size(X_population_point,1),1,365);
            X_dummy1(:,1,3:7:365)=1;
            X_dummy2(:,1,4:7:365)=1;
            X=cat(2,X_meteo_point,X_elevation_point,X_emission_point,X_population_point,X_dummy1,X_dummy2);
            X=double(X);
            sd_g.X_beta{2}=X;
            sd_g.X_beta_name{2}={'wind speed','pressure','temperature','elevation','emission','population','saturday','sunday'};
            
            clear X_meteo_point
            clear X_emission_point
            clear X_elevation_point
            clear X_population_point
            clear X
            clear X_dummy1
            clear X_dummy2
        end
    else
        sd_g.X_beta=[];
        sd_g.X_beta_name=[];
    end
    
    %X_z
    if flag_time_ground
        sd_g.X_z{1}=ones(n1,1,1);
        sd_g.X_z_name{1}={'constant'};
        if flag_bivariate
            sd_g.X_z{2}=ones(n2,1,1);
            sd_g.X_z_name{2}={'constant'};
        end
    else
        sd_g.X_z=[];
        sd_g.X_z_name=[];
    end
    
    if flag_w_ground
        x1_temp=sd_g.X_beta{1}(:,4,1);
        x2_temp=sd_g.X_beta{1}(:,5,1);
        x3_temp=sd_g.X_beta{1}(:,6,1);
        sd_g.X_g{1}=cat(4,x1_temp,x2_temp,x3_temp);
        sd_g.X_g_name{1}={'elevation','emission','population'};
        if flag_bivariate
            x1_temp=sd_g.X_beta{2}(:,4,1);
            x2_temp=sd_g.X_beta{2}(:,5,1);
            x3_temp=sd_g.X_beta{2}(:,6,1);
            sd_g.X_g{2}=cat(4,x1_temp,x2_temp,x3_temp);
            sd_g.X_g_name{2}={'elevation','emission','population'};
        end
    else
        sd_g.X_g=[];
        sd_g.X_g_name=[];
    end

    st_varset_g=stem_varset(sd_g.Y,sd_g.Y_name,sd_g.X_rg,sd_g.X_rg_name,sd_g.X_beta,sd_g.X_beta_name,sd_g.X_z,sd_g.X_z_name,sd_g.X_g,sd_g.X_g_name);
    
    if flag_tapering
        tapering_g=100; %km
        tapering_r=100; %km
    else
        tapering_g=[];
        tapering_r=[];
    end
    st_gridlist_g=stem_gridlist(tapering_g);
    st_grid=stem_grid([no2_ground.lat,no2_ground.lon],'deg','sparse','point');
    st_gridlist_g.add(st_grid);
    clear no2_ground
    if flag_bivariate
        st_grid=stem_grid([pm25_ground.lat,pm25_ground.lon],'deg','sparse','point');
        st_gridlist_g.add(st_grid);
        clear pm25_ground
    end
    
    % remote sensing data
    if flag_remote_data
        %Y
        load ../Demo/Data/no2_remote/no2_remote_025.mat
        
        no2_remote.data(no2_remote.data<1e+14)=NaN;
        no2_remote.data=reshape(no2_remote.data,size(no2_remote.data,1)*size(no2_remote.data,2),size(no2_remote.data,3));
        sd_r.Y{1}=no2_remote.data;
        sd_r.Y_name{1}='no2 remote';
        m1=size(sd_r.Y{1},1);
        
        if flag_bivariate
            load ../Demo/Data/aot_remote/aot_remote_025.mat
            aot_remote.data=squeeze(aot_remote.data(:,:,1,:));
            aot_remote.data=reshape(aot_remote.data,size(aot_remote.data,1)*size(aot_remote.data,2),size(aot_remote.data,3));
            sd_r.Y{2}=aot_remote.data;
            sd_r.Y_name{2}='aot remote';
            m2=size(sd_r.Y{2},1);
        end
        
        %X_rg
        sd_r.X_rg{1}=ones(m1,1,1);
        sd_r.X_rg_name{1}={'constant'};
        if flag_bivariate
            sd_r.X_rg{2}=ones(m2,1,1);
            sd_r.X_rg_name{2}={'constant'};
        end
        
        st_varset_r=stem_varset(sd_r.Y,sd_r.Y_name,sd_r.X_rg,sd_r.X_rg_name);
        
        st_gridlist_r=stem_gridlist(tapering_r);
        st_grid=stem_grid([no2_remote.lat(:),no2_remote.lon(:)],'deg','regular','pixel',size(no2_remote.lat),'square',0.25,0.25);
        st_gridlist_r.add(st_grid);
        clear no2_remote
        if flag_bivariate
            st_grid=stem_grid([aot_remote.lat(:),aot_remote.lon(:)],'deg','regular','pixel',size(aot_remote.lat),'square',0.25,0.25);
            st_gridlist_r.add(st_grid);
            clear aot_remote
        end
    else
        st_varset_r=[];
        st_gridlist_r=[];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      Model building     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    st_datestamp=stem_datestamp('01-01-2009','31-12-2009',T);
    %crossval
    if flag_crossval
        idx_step=2;
        st_crossval=stem_crossval('ground','no2',idx_step);
    else
        st_crossval=[];
    end
    
    if flag_remote_data
        if flag_remote_corr
            pixel_correlated=1;
        else
            pixel_correlated=0;
        end
    else
        pixel_correlated=[];
    end
    
    %data
    st_data=stem_data(st_varset_g,st_gridlist_g,st_varset_r,st_gridlist_r,st_datestamp,[],[],st_crossval,pixel_correlated);
    %par
    time_diagonal=1;
    st_par=stem_par(st_data,'exponential',pixel_correlated,time_diagonal);
    
    st_model=stem_model(st_data,st_par);
    st_model.set_system_size();
    clear sd_g
    clear sd_r
    
    %Data modification
    st_model.stem_data.log_transform;
    st_model.stem_data.standardize;
 
    %st_par initialization
    if flag_remote_data
        if flag_bivariate
            st_par.alpha_rg=[0.4 0.4 0.8 0.8]';
            if pixel_correlated
                st_par.theta_r=150;
                st_par.v_r=[1 0.6;0.6 1];
            else
                st_par.theta_r=[150 150]';
                st_par.v_r=eye(2);
            end
        else
            st_par.alpha_rg=[0.6 0.6]';
            st_par.theta_r=349;
            st_par.v_r=1;
        end
    end
    if flag_beta_ground
        st_par.beta=st_model.get_beta0();
    end
    if flag_w_ground
        if flag_bivariate
            st_par.alpha_g=[0.4 0.4;0.4 0.4;0.4 0.4]';
        else
            st_par.alpha_g=[0.6 0.6 0.6];
        end
        st_par.theta_g=[100 100 100]';
        for i=1:3
            if flag_bivariate
                v_g(:,:,i)=[1 0.6;0.6 1];
            else
                v_g(:,:,i)=1;
            end
        end
        st_par.v_g=v_g;
    end
    
    if flag_time_ground
        if flag_bivariate
            st_par.sigma_eta=diag(repmat(0.2,2,1));
            st_par.G=diag(repmat(0.8,2,1));
        else
            st_par.sigma_eta=diag(repmat(0.2,1,1));
            st_par.G=diag(repmat(0.8,1,1));
        end
    end
    
    if flag_remote_data
        if flag_bivariate
            st_par.sigma_eps=diag([0.3 0.3 0.3 0.3]);
        else
            st_par.sigma_eps=diag([0.3 0.3]);
        end
    else
        if flag_bivariate
            st_par.sigma_eps=diag([0.3 0.3]);
        else
            st_par.sigma_eps=0.3;
        end
    end
    st_model.set_initial_values(st_par);
    
    %Model estimation
    exit_toll=0.001;
    max_iterations=1;
    st_EM_options=stem_EM_options(exit_toll,max_iterations,'single');
    st_model.EM_estimate(st_EM_options);
    %st_model.set_varcov;
    %st_model.set_logL;    
end

if flag_kriging
    load ../Demo/Data/kriging/krig_elevation_005;
    krig_coordinates=[krig_elevation.lat(:),krig_elevation.lon(:)];
    krig_mask=krig_elevation.data_mask(:);

    st_krig=stem_krig(st_model);
    st_krig_grid=stem_grid(krig_coordinates,'deg','regular','pixel',[80,170],'square',0.05,0.05);
    back_transform=1;
    no_varcov=1;
    block_size=1000;
    X_krig='../Demo/Data/kriging/blocks';
    st_krig_result=st_krig.kriging('no2 ground',st_krig_grid,block_size,krig_mask,X_krig,back_transform,no_varcov);    
end

