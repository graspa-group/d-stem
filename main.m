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

%% ground level data
flag_parallel=0;
flag_remote_data=0;

flag_time_ground=0;
flag_time_remote=0;
flag_beta_ground=0;
flag_beta_remote=0;
flag_w_ground=1;
flag_w_remote=0;

flag_crossval=0;
flag_tapering=1;
flag_kriging=0;
flag_residuals=1;
pathparallel='/opt/matNfs/';

if 1
    %Y
    if not(flag_residuals)
        load ../Data/no2_ground.mat
        sd_g.Y{1}=no2_ground.data;
    else
        load ../Data/no2_ground.mat
        load ../Data/res_ground.mat
        sd_g.Y{1}=res_ground;
        clear res_ground
    end
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
    load ../Data/no2_meteo_point.mat
    load ../Data/no2_elevation_point.mat
    load ../Data/no2_emission_point.mat
    load ../Data/no2_population_point.mat
    X_elevation_point(isnan(X_elevation_point))=0;
    X_elevation_point=repmat(X_elevation_point,[1,1,T]);
    X_emission_point=repmat(X_emission_point,[1,1,T]);
    %take the log of the population count
    X_population_point(X_population_point==0)=0.02;
    X_population_point=log(X_population_point);
    X_population_point=repmat(X_population_point,[1,1,T]);
    X_lat_point=repmat(no2_ground.lat,[1,1,T]);
    X_lon_point=repmat(no2_ground.lon,[1,1,T]);
    X=cat(2,X_meteo_point,X_elevation_point,X_emission_point,X_population_point,X_lat_point,X_lon_point);
    X=double(X);
    
    if flag_residuals
        X([2328,2592],:,:)=[];
    end
    
    clear X_emission_point
    clear X_meteo_point
    clear X_elevation_point
    clear X_population_point
    
    if flag_beta_ground
        sd_g.X_beta_name{1}={'pressure','temperature','wind speed','elevation','emission','population','lat','lon'};
        sd_g.X_beta{1}=X(:,3:end,:);
    else
        sd_g.X_beta_name=[];
        sd_g.X_beta=[];
    end
    
    %X_time
    if flag_time_ground
        sd_g.X_time{1}=X(:,6:end,:);
        sd_g.X_time_name{1}={'elevation','emission','population','lat','lon'};
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
        sd_g.X_g{1}=X_new(:,1,:,3:8);
        sd_g.X_g_name{1}={'pressure','temperature','wind speed','elevation','emission','population'}; %'pressure','temperature','wind speed','elevation','emission','population'
        clear X_new
        clear temp
    else
        sd_g.X_g=[];
        sd_g.X_g_name=[];
    end
    clear X

    st_varset_g=stem_varset(sd_g.Y,sd_g.Y_name,sd_g.X_rg,sd_g.X_rg_name,sd_g.X_beta,sd_g.X_beta_name,sd_g.X_time,sd_g.X_time_name,sd_g.X_g,sd_g.X_g_name);
    
    if flag_tapering
        tapering_g=300; %km
        tapering_r=100; %km
    else
        tapering_g=[];
        tapering_r=[];
    end
    st_gridlist_g=stem_gridlist(tapering_g);
    if flag_residuals
        no2_ground.lat([2328,2592])=[];
        no2_ground.lon([2328,2592])=[];
    end
    st_grid=stem_grid([no2_ground.lat,no2_ground.lon],'deg','sparse','point');
    st_gridlist_g.add(st_grid);
    clear no2_ground
    
    
    %% remote sensing data
    
    if flag_remote_data
        %Y
        if not(flag_residuals)
            load ../Data/no2_remote_025.mat
            no2_remote.data(no2_remote.data<1e+14)=NaN;
            no2_remote.data=reshape(no2_remote.data,size(no2_remote.data,1)*size(no2_remote.data,2),size(no2_remote.data,3));
            sd_r.Y{1}=no2_remote.data;
        else
            load ../Data/no2_remote_025.mat
            load ../Data/res_remote
            sd_r.Y{1}=res_remote;
            clear res_remote
        end
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
    
    %% model building
    %datestamp
    st_datestamp=stem_datestamp('01-01-2009','31-12-2009',T);
    %crossval
    if flag_crossval
        idx_step=3;
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
    time_diagonal=0;
    st_par=stem_par(st_data,'exponential',remote_correlated,time_diagonal);
    
    st_model=stem_model(st_data,st_par);
    st_model.set_system_size();
    clear sd_g
    clear sd_r
    
    %% data modification
    st_model.stem_data.space_crop([44,47,7,14]);
    st_model.stem_data.time_crop(1:10);
    %st_model.stem_data.log_transform;
    st_model.stem_data.standardize;
    
    %% st_par initialization
    
    if flag_remote_data
        st_par.alpha_rg=[0.4 0.7]';
        st_par.theta_r=400;
        st_par.v_r=1;
    end
    if flag_beta_ground
        st_par.beta=st_model.get_beta0();
    end
    if flag_w_ground
        st_par.alpha_g=[0.15 0.15 0.15 0.15 0.15 0.15];
        st_par.theta_g=[250 250 250 250 250 250]';
        for i=1:6
            v_g(:,:,i)=1;
        end
        st_par.v_g=v_g;
    end
    
    if flag_time_ground||flag_time_remote
        st_par.sigma_eta=diag(repmat(0.1,6,1));
        st_par.G=diag(repmat(0.7,6,1));
    end
    
    st_par.sigma_eps=diag([0.3 0.15]);
    st_model.stem_par_initial=st_par;
else
    load ../Data/st_model_small_area_residuals.mat
end

%% model estimation
st_EM_options=stem_EM_options(0.0001,100,'single',[],0,[]);
if flag_parallel
    st_EM_options.pathparallel=pathparallel;
end

st_model.stem_par=st_model.stem_par_initial;
st_sim=stem_sim(st_model);
st_sim.simulate;

st_model.EM_estimate(st_EM_options);
%st_model.set_Hessian;
%st_model.set_logL;

save(['st_model_',datestr(now,'yyyymmdd_HHMMSS')],'st_model');

if flag_kriging
    X_krig.name={'constant','elevation','emission','population'};
    X_krig.date_stamp=st_model.stem_data.stem_datestamp;
    back_transform=1;
    no_varcov=1;
    st_krig_result=st_krig.kriging(block_size,'no2',st_grid_krig,mask_krig,X_krig,back_transform,no_varcov);
    save st_krig_result st_krig_result -v7.3
end

