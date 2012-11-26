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
flag_remote_data=0;

flag_time_ground=1;
flag_time_remote=0;
flag_beta_ground=1;
flag_beta_remote=0;
flag_w_ground=1;
flag_w_remote=0;

flag_crossval=1;
flag_tapering=0;

flag_estimate=1;
flag_kriging=0;

pathparallel='/opt/matNfs/';
if flag_estimate
    %Y
    %load the full dataset in order to have the information of the lat-lon
    %of all the NO2 monitoring stations
    load ../Data/no2_ground.mat
    lat_temp=no2_ground.lat;
    lon_temp=no2_ground.lon;
    
    %load the redeced dataset
    load ../Data/no2_ground_background_suburban.mat
    sd_g.Y{1}=no2_ground.data;
    sd_g.Y_name={'no2'};
    N=size(sd_g.Y{1},1);
    T=size(sd_g.Y{1},2);
    
    %X_rg
    if flag_w_remote
        sd_g.X_rg{1}=ones(size(sd_g.Y{1},1),1);
        sd_g.X_rg_name{1}={'constant'};
    else
        sd_g.X_rg=[];
        sd_g.X_rg_name{1}=[];
    end
    
    %X_beta
    %load the covariates over all the NO2 monitoring stations
    load ../Data/no2_meteo_point.mat
    load ../Data/no2_elevation_point.mat
    load ../Data/no2_emission_point.mat
    load ../Data/no2_population_point.mat
    load ../Data/no2_no2year_point.mat
    X_elevation_point(isnan(X_elevation_point))=0;
    X_elevation_point=repmat(X_elevation_point,[1,1,T]);
    X_emission_point=repmat(X_emission_point,[1,1,T]);
    X_no2year_point=repmat(no2_y,[1,1,T]);
    
    %extract the indices of the reduced dataset
    idx=[];
    for i=1:length(no2_ground.lat)
        finded=0;
        for j=1:length(lat_temp)
            if (no2_ground.lat(i)==lat_temp(j))&&(no2_ground.lon(i)==lon_temp(j))&&(finded==0)
                idx=[idx;j];
                finded=1;
            end
        end
    end
    
    X_dummy1=zeros(size(X_population_point,1),1,365);
    X_dummy2=zeros(size(X_population_point,1),1,365);
    X_dummy1(:,1,3:7:365)=1;
    X_dummy2(:,1,4:7:365)=1;
    %take the log of the population count
    X_population_point(X_population_point==0)=0.02;
    X_population_point=log(X_population_point);
    X_population_point=repmat(X_population_point,[1,1,T]);
    X_no2_remote_point=repmat(no2_y,[1,1,T]);
    X=cat(2,X_no2_remote_point,X_meteo_point,X_elevation_point,X_emission_point,X_population_point,X_dummy1,X_dummy2);
    X=double(X);
    %obtain the covariates over the reduced dataset
    X=X(idx,:,:);
    
    %exclude the sites above 1000 meters
    X_elevation_point=X_elevation_point(idx,1,1);
    idx=find(X_elevation_point>=1000);
    
    clear X_emission_point
    clear X_meteo_point
    clear X_elevation_point
    clear X_population_point
    clear X_no2_remote_point
    
    if flag_beta_ground
        sd_g.X_beta{1}=X(:,[5,6,8,9,10,11],:);
        sd_g.X_beta_name{1}={'temperature','wind speed','emission','population','dummy1','dummy2'};%{'pressure','temperature','wind speed','elevation','emission','population','dummy'};
        %sd_g.X_beta{1}=ones(size(X,1),1,1);
        %sd_g.X_beta_name{1}={'constant'};
    else
        sd_g.X_beta=[];
        sd_g.X_beta_name=[];
    end
    
    %X_time
    if flag_time_ground
        %sd_g.X_time{1}=X(:,[6],:);
        %sd_g.X_time_name{1}={'wind speed'};
        sd_g.X_time{1}=ones(size(X,1),1,1);
        sd_g.X_time_name{1}={'constant'};
    else
        sd_g.X_time=[];
        sd_g.X_time_name=[];
    end
    
    if flag_w_ground
        for t=1:size(X,3)
            temp=X(:,:,t);
            temp=reshape(temp,size(temp,1),1,1,size(temp,2));
            X_new(:,1,t,:)=temp;
        end
        %sd_g.X_g{1}=X_new(:,1,:,5:6);
        %sd_g.X_g_name{1}={'temperature','wind speed'};
        sd_g.X_g{1}=ones(size(X,1),1,1,1);
        sd_g.X_g_name{1}={'constant'};
        clear X_new
        clear temp
    else
        sd_g.X_g=[];
        sd_g.X_g_name=[];
    end
    clear X
    
    st_varset_g=stem_varset(sd_g.Y,sd_g.Y_name,sd_g.X_rg,sd_g.X_rg_name,sd_g.X_beta,sd_g.X_beta_name,sd_g.X_time,sd_g.X_time_name,sd_g.X_g,sd_g.X_g_name);
    
    if flag_tapering
        tapering_g=350; %km
        tapering_r=100; %km
    else
        tapering_g=[];
        tapering_r=[];
    end
    st_gridlist_g=stem_gridlist(tapering_g);
    st_grid=stem_grid([no2_ground.lat,no2_ground.lon],'deg','sparse','point');
    st_gridlist_g.add(st_grid);
    clear no2_ground
    
    % remote sensing data
    if flag_remote_data
        %Y
        load ../Data/no2_remote_025.mat
        no2_remote.data(no2_remote.data<1e+14)=NaN;
        no2_remote.data=reshape(no2_remote.data,size(no2_remote.data,1)*size(no2_remote.data,2),size(no2_remote.data,3));
        sd_r.Y{1}=no2_remote.data;
        sd_r.Y_name={'no2'};
        N=size(sd_r.Y{1},1);
        
        %X_rg
        if flag_w_remote
            sd_r.X_rg{1}=ones(N,1,1);
            sd_r.X_rg_name{1}={'constant'};
        else
            sd_r.X_rg=[];
            sd_r.X_rg_name{1}=[];
        end
        
        %X_beta
        if flag_beta_remote
            load ../Data/no2_meteo_pixel.mat
            load ../Data/no2_emission_pixel.mat
            X_emission_pixel=repmat(X_emission_pixel,[1,1,T]);
            X=cat(2,X_meteo_pixel,X_emission_pixel);
            sd_r.X_beta{1}=ones(N,1,1);%X;
            sd_r.X_beta_name{1}={'constant'};
        else
            sd_r.X_beta=[];
            sd_r.X_beta_name=[];
        end
        
        %X_time
        if flag_time_remote
            sd_r.X_time{1}=ones(N,1,1);
            sd_r.X_time_name{1}={'constant'};
        else
            sd_r.X_time=[];
            sd_r.X_time_name=[];
        end
        
        st_varset_r=stem_varset(sd_r.Y,sd_r.Y_name,sd_r.X_rg,sd_r.X_rg_name,sd_r.X_beta,sd_r.X_beta_name,sd_r.X_time,sd_r.X_time_name);
        
        st_gridlist_r=stem_gridlist(tapering_r);
        st_grid=stem_grid([no2_remote.lat(:),no2_remote.lon(:)],'deg','regular','pixel',[128 184],'square',0.25,0.25);
        st_gridlist_r.add(st_grid);
        clear no2_remote
    else
        st_varset_r=[];
        st_gridlist_r=[];
    end
    
    % model building
    %datestamp
    st_datestamp=stem_datestamp('01-01-2009','31-12-2009',T);
    %crossval
    if flag_crossval
        idx_step=2;
        st_crossval=stem_crossval('ground','no2',idx_step);
    else
        st_crossval=[];
    end
    %data
    st_data=stem_data(st_varset_g,st_gridlist_g,st_varset_r,st_gridlist_r,st_datestamp,[],[],st_crossval);
    %par
    if flag_remote_data
        remote_correlated=0;
    else
        remote_correlated=[];
    end
    time_diagonal=1;
    st_par=stem_par(st_data,'exponential',remote_correlated,time_diagonal);
    
    st_model=stem_model(st_data,st_par);
    st_model.set_system_size();
    clear sd_g
    clear sd_r
    
    % data modification
    %st_model.stem_data.space_crop([44,47,6,14]);
    %st_model.stem_data.time_crop(1:365);
    st_model.stem_data.site_crop('ground','no2',idx);
    st_model.stem_data.log_transform;
    st_model.stem_data.standardize;
    
    % st_par initialization
    if flag_remote_data
        st_par.alpha_rg=[0.4 0.7]';
        st_par.theta_r=400;
        st_par.v_r=1;
    end
    if flag_beta_ground
        st_par.beta=st_model.get_beta0();
    end
    if flag_w_ground
        st_par.alpha_g=[0.4];
        st_par.theta_g=[400]';
        for i=1:1
            v_g(:,:,i)=1;
        end
        st_par.v_g=v_g;
    end
    
    if flag_time_ground||flag_time_remote
        st_par.sigma_eta=diag(repmat(0.047,1,1));
        st_par.G=diag(repmat(0.8,1,1));
    end
    
    st_par.sigma_eps=diag([0.3]);
    st_model.set_initial_values(st_par);
    
    % model estimation
    st_EM_options=stem_EM_options(0.001,50,'single',[],0,[]);
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
    load ../Data/no2_krig_coordinates_005.mat;
    load ../Data/no2_krigmask_005.mat
    st_krig=stem_krig(st_model);
    st_krig_grid=stem_grid([krig_lat(:),krig_lon(:)],'deg','regular','pixel',[640,920],'square',0.05,0.05);
    back_transform=1;
    no_varcov=0;
    block_size=10000;
    mask=krig_mask;
    X_krig='../Data/blocks/';
    st_krig_result=st_krig.kriging('no2',st_krig_grid,block_size,mask,X_krig,back_transform,no_varcov);
    save st_krig_result st_krig_result -v7.3
end

