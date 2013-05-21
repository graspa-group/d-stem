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

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Data  building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the no2 ground level observations
load ../Demo/Data/no2_ground/no2_ground_background.mat
sd_p.Y{1}=no2_ground.data;
sd_p.Y_name{1}='no2 ground';
n1=size(sd_p.Y{1},1);
T=size(sd_p.Y{1},2);

%load the pm25 ground level observations
load ../Demo/Data/pm25_ground/pm25_ground_background.mat
sd_p.Y{2}=pm25_ground.data;
sd_p.Y_name{2}='pm2.5 ground';
n2=size(sd_p.Y{2},1);

%No downscaler
sd_p.X_bp=[];
sd_p.X_bp_name=[];

%X_beta
%load the covariates for the NO2 monitoring stations
load ../Demo/Data/no2_ground/no2_ground_covariates.mat
sd_p.X_beta{1}=X;
sd_p.X_beta_name{1}={'wind speed','pressure','temperature','elevation','emission','population','saturday','sunday'};

%load the covariates for the PM2.5 monitoring stations
load ../Demo/Data/pm25_ground/pm25_ground_covariates.mat
sd_p.X_beta{2}=X;
sd_p.X_beta_name{2}={'wind speed','pressure','temperature','elevation','emission','population','saturday','sunday'};

%X_z
sd_p.X_z{1}=ones(n1,1,1);
sd_p.X_z_name{1}={'constant'};
sd_p.X_z{2}=ones(n2,1,1);
sd_p.X_z_name{2}={'constant'};

x1_temp=sd_p.X_beta{1}(:,4,1);
x2_temp=sd_p.X_beta{1}(:,6,1);
X=cat(4,x1_temp,x2_temp);
sd_p.X_p{1}=X;
sd_p.X_p_name{1}={'elevation','population'};
x1_temp=sd_p.X_beta{2}(:,4,1);
x2_temp=sd_p.X_beta{2}(:,6,1);
X=cat(4,x1_temp,x2_temp);
sd_p.X_p{2}=X;
sd_p.X_p_name{2}={'elevation','population'};

st_varset_p=stem_varset(sd_p.Y,sd_p.Y_name,sd_p.X_bp,sd_p.X_bp_name,sd_p.X_beta,sd_p.X_beta_name,sd_p.X_z,sd_p.X_z_name,sd_p.X_p,sd_p.X_p_name);

st_gridlist_p=stem_gridlist();

sd_p.coordinates=[no2_ground.lat,no2_ground.lon];
st_grid=stem_grid(sd_p.coordinates,'deg','sparse','point');
st_gridlist_p.add(st_grid);
sd_p.coordinates=[pm25_ground.lat,pm25_ground.lon];
st_grid=stem_grid(sd_p.coordinates,'deg','sparse','point');
st_gridlist_p.add(st_grid);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Model building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

st_datestamp=stem_datestamp('01-01-2009','31-12-2009',T);

%stem_data object creation
st_data=stem_data(st_varset_p,st_gridlist_p,[],[],st_datestamp);
%stem_par object creation
flag_time_diagonal=0;
st_par=stem_par(st_data,'exponential',[],flag_time_diagonal);
%stem_model object creation
st_model=stem_model(st_data,st_par);
clear sd_p

%Data transform
st_model.stem_data.log_transform;
st_model.stem_data.standardize;

%st_par object initialization

st_par.beta=st_model.get_beta0();
st_par.alpha_p=[0.4 0.4;0.4 0.4]';
st_par.theta_p=[100 100]';
for i=1:length(st_par.theta_p)
    v_p(:,:,i)=[1 0.6;0.6 1];
end
st_par.v_p=v_p;
st_par.sigma_eta=diag([0.2 0.2]);
st_par.G=diag([0.8 0.8]);
st_par.sigma_eps=diag([0.3 0.3]);

st_model.set_initial_values(st_par);

%Model estimation
exit_toll=0.002;
max_iterations=100;
st_EM_options=stem_EM_options(exit_toll,max_iterations);
st_model.EM_estimate(st_EM_options);
st_model.set_varcov;
st_model.set_logL;

load ../Demo/Data/kriging/krig_elevation_005;
krig_coordinates=[krig_elevation.lat(:),krig_elevation.lon(:)];
krig_mask=krig_elevation.data_mask(:);
%kriging
st_krig=stem_krig(st_model);
st_krig_grid=stem_grid(krig_coordinates,'deg','regular','pixel',[80,170],'square',0.05,0.05);
back_transform=1;
no_varcov=0;
block_size=1000;
X_krig='../Demo/Data/kriging/blocks';
st_krig_result=st_krig.kriging('no2 ground',st_krig_grid,block_size,krig_mask,X_krig,back_transform,no_varcov);
        
