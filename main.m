%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
% ground level data
flag_parallel=0;
flag_remote_data=1;

flag_time_ground=1;
flag_time_remote=0;
flag_beta_ground=1;
flag_beta_remote=0;
flag_w_ground=1;

flag_crossval=0;
flag_tapering=0;

flag_estimate=0;
flag_kriging=1;

pathparallel='/opt/matNfs/';
if flag_estimate
    %load the no2 ground level observations
    load ../Data/no2_ground/no2_ground_background.mat
    sd_g.Y{1}=no2_ground.data;
    sd_g.Y_name{1}='no2 ground';
    n1=size(sd_g.Y{1},1);
    T=size(sd_g.Y{1},2);
    
    %load the pm25 ground level observations
    load ../Data/pm25_ground/pm25_ground_background.mat
    sd_g.Y{2}=pm25_ground.data;
    sd_g.Y_name{2}='pm2.5 ground';
    n2=size(sd_g.Y{2},1);
    
    %X_rg
    if flag_remote_data
        sd_g.X_rg{1}=ones(n1,1);
        sd_g.X_rg_name{1}={'constant'};
        sd_g.X_rg{2}=ones(n2,1);
        sd_g.X_rg_name{2}={'constant'};        
    else
        sd_g.X_rg=[];
        sd_g.X_rg_name{1}=[];
    end

    if flag_beta_ground
        %X_beta
        %load the covariates over all the NO2 monitoring stations
        load ../Data/no2_ground/no2_background_meteo_point.mat
        load ../Data/no2_ground/no2_background_elevation_point.mat
        load ../Data/no2_ground/no2_background_emission_point.mat
        load ../Data/no2_ground/no2_background_population_point.mat
        X_elevation_point(isnan(X_elevation_point))=0;
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
        sd_g.X_beta_name{1}={'wind speed','temperature','pressure','elevation','emission','population','saturday','sunday'};
        
        clear X_meteo_point
        clear X_emission_point
        clear X_elevation_point
        clear X_population_point
        
        load ../Data/pm25_ground/pm25_background_meteo_point.mat
        load ../Data/pm25_ground/pm25_background_elevation_point.mat
        load ../Data/pm25_ground/pm10_background_emission_point.mat
        load ../Data/pm25_ground/pm25_background_population_point.mat
        X_elevation_point(isnan(X_elevation_point))=0;
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
        sd_g.X_beta_name{2}={'wind speed','temperature','pressure','elevation','emission','population','saturday','sunday'};
        
        clear X_meteo_point
        clear X_emission_point
        clear X_elevation_point
        clear X_population_point   
        clear X
    else
        sd_g.X_beta=[];
        sd_g.X_beta_name=[];
    end
    
    %X_time
    if flag_time_ground
        sd_g.X_time{1}=ones(n1,1,1);
        sd_g.X_time_name{1}={'constant'};
        sd_g.X_time{2}=ones(n2,1,1);
        sd_g.X_time_name{2}={'constant'};        
    else
        sd_g.X_time=[];
        sd_g.X_time_name=[];
    end
    
    if flag_w_ground
        x1_temp=sd_g.X_beta{1}(:,4,1);
        x2_temp=sd_g.X_beta{1}(:,5,1);
        x3_temp=sd_g.X_beta{1}(:,6,1);
        sd_g.X_g{1}=cat(4,x1_temp,x2_temp,x3_temp);
        sd_g.X_g_name{1}={'elevation','emission','population'};
        x1_temp=sd_g.X_beta{2}(:,4,1);
        x2_temp=sd_g.X_beta{2}(:,5,1);
        x3_temp=sd_g.X_beta{2}(:,6,1);
        sd_g.X_g{2}=cat(4,x1_temp,x2_temp,x3_temp);
        sd_g.X_g_name{2}={'elevation','emission','population'};        
%         sd_g.X_g{1}=ones(n1,1,1,1);
%         sd_g.X_g_name{1}={'constant'};
%         sd_g.X_g{2}=ones(n2,1,1,1);
%         sd_g.X_g_name{2}={'constant'};        
    else
        sd_g.X_g=[];
        sd_g.X_g_name=[];
    end

    st_varset_g=stem_varset(sd_g.Y,sd_g.Y_name,sd_g.X_rg,sd_g.X_rg_name,sd_g.X_beta,sd_g.X_beta_name,sd_g.X_time,sd_g.X_time_name,sd_g.X_g,sd_g.X_g_name);
    
    if flag_tapering
        tapering_g=400; %km
        tapering_r=250; %km
    else
        tapering_g=[];
        tapering_r=[];
    end
    st_gridlist_g=stem_gridlist(tapering_g);
    st_grid=stem_grid([no2_ground.lat,no2_ground.lon],'deg','sparse','point');
    st_gridlist_g.add(st_grid);
    clear no2_ground
    st_grid=stem_grid([pm25_ground.lat,pm25_ground.lon],'deg','sparse','point');
    st_gridlist_g.add(st_grid);
    clear pm25_ground    
    
    % remote sensing data
    if flag_remote_data
        %Y
        load ../Data/no2_remote/no2_remote_025.mat
        if 1
           lat=no2_remote.lat;
           lon=no2_remote.lon;
           Llat=(lat>36)&(lat<47.5);
           Llon=(lon>-10)&(lon<19);
           L=Llat&Llon;
           H=sum(L,1);
           K=sum(L,2);
           idx1=find(H);
           idx2=find(K);
           no2_remote.lat=no2_remote.lat(idx2,idx1);
           no2_remote.lon=no2_remote.lon(idx2,idx1);
           no2_remote.data=no2_remote.data(idx2,idx1,:);
        end
        
        no2_remote.data(no2_remote.data<1e+14)=NaN;
        no2_remote.data=reshape(no2_remote.data,size(no2_remote.data,1)*size(no2_remote.data,2),size(no2_remote.data,3));
        sd_r.Y{1}=no2_remote.data;
        sd_r.Y_name{1}='no2 remote';
        m1=size(sd_r.Y{1},1);
        
        load ../Data/aot_remote/aot_remote_025.mat
        if 1
           lat=aot_remote.lat;
           lon=aot_remote.lon;
           Llat=(lat>36)&(lat<47.5);
           Llon=(lon>-10)&(lon<19);
           L=Llat&Llon;
           H=sum(L,1);
           K=sum(L,2);
           idx1=find(H);
           idx2=find(K);
           aot_remote.lat=aot_remote.lat(idx2,idx1);
           aot_remote.lon=aot_remote.lon(idx2,idx1);
           aot_remote.data=aot_remote.data(idx2,idx1,:,:);
        end        
        aot_remote.data=squeeze(aot_remote.data(:,:,1,:));
        aot_remote.data=reshape(aot_remote.data,size(aot_remote.data,1)*size(aot_remote.data,2),size(aot_remote.data,3));
        sd_r.Y{2}=aot_remote.data;
        sd_r.Y_name{2}='aot remote';
        m2=size(sd_r.Y{2},1);        
        
        %X_rg
        sd_r.X_rg{1}=ones(m1,1,1);
        sd_r.X_rg_name{1}={'constant'};
        sd_r.X_rg{2}=ones(m2,1,1);
        sd_r.X_rg_name{2}={'constant'};        
        
        st_varset_r=stem_varset(sd_r.Y,sd_r.Y_name,sd_r.X_rg,sd_r.X_rg_name);
        
        st_gridlist_r=stem_gridlist(tapering_r);
        st_grid=stem_grid([no2_remote.lat(:),no2_remote.lon(:)],'deg','regular','pixel',size(no2_remote.lat),'square',0.25,0.25);
        st_gridlist_r.add(st_grid);
        clear no2_remote
        st_grid=stem_grid([aot_remote.lat(:),aot_remote.lon(:)],'deg','regular','pixel',size(aot_remote.lat),'square',0.25,0.25);
        st_gridlist_r.add(st_grid);
        clear aot_remote        
    else
        st_varset_r=[];
        st_gridlist_r=[];
    end
    
    %Model building
    
    %datestamp
    st_datestamp=stem_datestamp('01-01-2009','31-12-2009',T);
    %crossval
    if flag_crossval
        idx_step=2;
        st_crossval=stem_crossval('ground','no2',idx_step);
    else
        st_crossval=[];
    end
    
    if flag_remote_data
        remote_correlated=0;
    else
        remote_correlated=[];
    end
    
    %data
    st_data=stem_data(st_varset_g,st_gridlist_g,st_varset_r,st_gridlist_r,st_datestamp,[],[],st_crossval,remote_correlated);
    %par
    time_diagonal=1;
    st_par=stem_par(st_data,'exponential',remote_correlated,time_diagonal);
    
    st_model=stem_model(st_data,st_par);
    st_model.set_system_size();
    clear sd_g
    clear sd_r
    
    % data modification
    %st_model.stem_data.space_crop([36,47.5,-10,19]); %da croppare perché è già croppato solo il remote
    st_model.stem_data.space_crop([44,47.5,6,14]); %Norther-Italy
    st_model.stem_data.time_crop(180:180+90);
    %st_model.stem_data.site_crop('ground','no2',idx);
    st_model.stem_data.log_transform;
    st_model.stem_data.standardize;

        
    % st_par initialization
    if flag_remote_data
        st_par.alpha_rg=[0.4 0.4 0.8 0.8]';
        st_par.theta_r=[150 150]';
        st_par.v_r=eye(2);
    end
    if flag_beta_ground
        st_par.beta=st_model.get_beta0();
    end
    if flag_w_ground
        st_par.alpha_g=[0.4 0.4;0.4 0.4;0.4 0.4]';
        st_par.theta_g=[200 200 200]';
        for i=1:3
            v_g(:,:,i)=[1 0.6;0.6 1];
        end
        st_par.v_g=v_g;
    end
    
    if flag_time_ground||flag_time_remote
        st_par.sigma_eta=diag(repmat(0.2,2,1));
        st_par.G=diag(repmat(0.8,2,1));
    end
    
    if flag_remote_data
        st_par.sigma_eps=diag([0.3 0.3 0.3 0.3]);
    else
        st_par.sigma_eps=diag([0.3 0.3]);
    end
    st_model.set_initial_values(st_par);
    % model estimation
    st_EM_options=stem_EM_options(0.001,100,'single',[],0,[]);
    if flag_parallel
        st_EM_options.pathparallel=pathparallel;
    end
    
    %st_model.stem_par=st_model.stem_par_initial;
    %st_sim=stem_sim(st_model);
    %st_sim.simulate;
    
    st_model.EM_estimate(st_EM_options);
    %st_model.set_varcov;
    %st_model.set_logL;
    
    save(['st_model_',datestr(now,'yyyymmdd_HHMMSS')],'st_model');
end

if flag_kriging
    load ../Data/output/model/st_model_20130204_030938.mat;
    load('C:\Francesco\My Dropbox\air_quality_code_and_data\Output\krig_all_005_southwest_europe.mat');
    krig_lat=krig_meteo_005.lat;
    krig_lon=krig_meteo_005.lon;
    krig_mask=krig_elevation_005.data_mask(:);
    clear krig_elevation_005
    clear krig_emission_005
    clear krig_meteo_005
    clear krig_population_005

    st_krig=stem_krig(st_model);
    st_krig_grid=stem_grid([krig_lat(:),krig_lon(:)],'deg','regular','pixel',[230,580],'square',0.05,0.05);
    back_transform=1;
    no_varcov=0;
    block_size=500;
    X_krig='../Data/blocks/';
    st_krig_result=st_krig.kriging('no2 ground',st_krig_grid,block_size,krig_mask,X_krig,back_transform,no_varcov);
    save st_krig_result st_krig_result -v7.3
end

