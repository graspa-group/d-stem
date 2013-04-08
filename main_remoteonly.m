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
% ground level data
flag_parallel=0;
flag_remote_data=0;

flag_time_ground=0;
flag_time_remote=0;
flag_beta_ground=0;
flag_beta_remote=0;
flag_w_ground=1;

flag_crossval=0;
flag_tapering=0;

flag_estimate=1;
flag_kriging=0;

pathparallel='/opt/matNfs/';
if flag_estimate    
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
    sd_g.Y{1}=no2_remote.data;
    sd_g.Y_name{1}='no2 remote';
    n1=size(sd_g.Y{1},1);
    
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
    sd_g.Y{2}=aot_remote.data;
    sd_g.Y_name{2}='aot remote';
    n2=size(sd_g.Y{2},1);
    
    %X_rg
    if flag_remote_data
      
    else
        sd_g.X_rg=[];
        sd_g.X_rg_name{1}=[];
    end

    if flag_beta_ground

    else
        sd_g.X_beta=[];
        sd_g.X_beta_name=[];
    end
    
    %X_z
    if flag_time_ground
       
    else
        sd_g.X_z=[];
        sd_g.X_z_name=[];
    end
    
    if flag_w_ground
        sd_g.X_g{1}=ones(n1,1,1,1);
        sd_g.X_g_name{1}={'constant'};
        sd_g.X_g{2}=ones(n2,1,1,1);
        sd_g.X_g_name{2}={'constant'};        
    else
        sd_g.X_g=[];
        sd_g.X_g_name=[];
    end

    st_varset_g=stem_varset(sd_g.Y,sd_g.Y_name,sd_g.X_rg,sd_g.X_rg_name,sd_g.X_beta,sd_g.X_beta_name,sd_g.X_z,sd_g.X_z_name,sd_g.X_g,sd_g.X_g_name);
    
    if flag_tapering

    else
        tapering_g=[];
        tapering_r=[];
    end
    
    st_gridlist_g=stem_gridlist(tapering_g);
    st_grid=stem_grid([no2_remote.lat(:),no2_remote.lon(:)],'deg','sparse','point');
    st_gridlist_g.add(st_grid);
    clear no2_remote
    st_grid=stem_grid([aot_remote.lat(:),aot_remote.lon(:)],'deg','sparse','point');
    st_gridlist_g.add(st_grid);
    clear aot_remote
   
    % remote sensing data
    if flag_remote_data

    else
        st_varset_r=[];
        st_gridlist_r=[];
    end
    
    %Model building
    
    %datestamp
    st_datestamp=stem_datestamp('01-01-2009','31-12-2009',365);
    %crossval
    if flag_crossval

    else
        st_crossval=[];
    end
    
    if flag_remote_data

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

    end
    if flag_beta_ground

    end
    if flag_w_ground
        st_par.alpha_g=[0.8 0.8]';
        st_par.theta_g=[230]';
        for i=1:1
            v_g(:,:,i)=[1 0.0;0.0 1];
        end
        st_par.v_g=v_g;
    end
    
    if flag_time_ground||flag_time_remote

    end
    
    if flag_remote_data

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

