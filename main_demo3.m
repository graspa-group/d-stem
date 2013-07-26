%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% D-STEM - Distributed Space Time Expecation Maximization              %
%%%                                                                      %
%%% Author: Francesco Finazzi                                            %
%%% E-mail: francesco.finazzi@unibg.it                                   %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Economics and Quantitative Methods %
%%% Author website: http://www.unibg.it/pers/?francesco.finazzi          %
%%% Code website: https://code.google.com/p/d-stem/                      %
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

clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Ground level data  building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the no2 ground level observations
load ../Demo/Data/no2_ground.mat
ground.Y{1}=no2_ground.data;
ground.Y_name{1}='no2 ground';
n1=size(ground.Y{1},1);
T=size(ground.Y{1},2);

%load the pm25 ground level observations
load ../Demo/Data/pm25_ground.mat
ground.Y{2}=pm25_ground.data;
ground.Y_name{2}='pm2.5 ground';
n2=size(ground.Y{2},1);

%downscaler loading coefficients
ground.X_bp{1}=ones(n1,1);
ground.X_bp_name{1}={'constant'};
ground.X_bp{2}=ones(n2,1);
ground.X_bp_name{2}={'constant'};

%X_beta
%load the covariates for the NO2 monitoring stations
load ../Demo/Data/no2_ground_covariates.mat
ground.X_beta{1}=X;
ground.X_beta_name{1}={'wind speed','elevation','sunday'};

%load the covariates for the PM2.5 monitoring stations
load ../Demo/Data/pm25_ground_covariates.mat
ground.X_beta{2}=X;
ground.X_beta_name{2}={'wind speed','elevation','sunday'};

%X_z
ground.X_z{1}=ones(n1,1,1);
ground.X_z_name{1}={'constant'};
ground.X_z{2}=ones(n2,1,1);
ground.X_z_name{2}={'constant'};

ground.X_p{1}=ground.X_beta{1}(:,2,1);
ground.X_p_name{1}={'elevation'};
ground.X_p{2}=ground.X_beta{2}(:,2,1);
ground.X_p_name{2}={'elevation'};

obj_stem_varset_p=stem_varset(ground.Y,ground.Y_name,ground.X_bp,ground.X_bp_name,ground.X_beta,ground.X_beta_name,ground.X_z,ground.X_z_name,ground.X_p,ground.X_p_name);

tapering_p=[];
obj_stem_gridliobj_stem_p=stem_gridlist(tapering_p);

ground.coordinates=[no2_ground.lat,no2_ground.lon];
obj_stem_grid=stem_grid(ground.coordinates,'deg','sparse','point');
obj_stem_gridliobj_stem_p.add(obj_stem_grid);
ground.coordinates=[pm25_ground.lat,pm25_ground.lon];
obj_stem_grid=stem_grid(ground.coordinates,'deg','sparse','point');
obj_stem_gridliobj_stem_p.add(obj_stem_grid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Remote sensing data  building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Y
load ../Demo/Data/no2_remote_025.mat
remote.Y{1}=no2_remote.data;
remote.Y_name{1}='no2 remote';
m1=size(remote.Y{1},1);

load ../Demo/Data/aot_remote_025.mat
remote.Y{2}=aot_remote.data;
remote.Y_name{2}='aot remote';
m2=size(remote.Y{2},1);
    
%X_bp
remote.X_bp{1}=ones(m1,1);
remote.X_bp_name{1}={'constant'};
remote.X_bp{2}=ones(m2,1);
remote.X_bp_name{2}={'constant'};

obj_stem_varset_b=stem_varset(remote.Y,remote.Y_name,remote.X_bp,remote.X_bp_name);

obj_stem_gridliobj_stem_b=stem_gridlist();
remote.coordinates=[no2_remote.lat(:),no2_remote.lon(:)];
obj_stem_grid=stem_grid(remote.coordinates,'deg','regular','pixel',size(no2_remote.lat),'square',0.25,0.25);
obj_stem_gridliobj_stem_b.add(obj_stem_grid);
remote.coordinates=[aot_remote.lat(:),aot_remote.lon(:)];
obj_stem_grid=stem_grid(remote.coordinates,'deg','regular','pixel',size(aot_remote.lat),'square',0.25,0.25);
obj_stem_gridliobj_stem_b.add(obj_stem_grid);
clear aot_remote
clear no2_remote

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Model building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj_stem_datestamp=stem_datestamp('01-01-2009 00:00','31-12-2009 00:00',T);

flag_pixel_correlated=0;
flag_time_diagonal=0;

%stem_data object creation
shape=shaperead('../Demo/Maps/worldmap');
obj_stem_data=stem_data(obj_stem_varset_p,obj_stem_gridliobj_stem_p,obj_stem_varset_b,obj_stem_gridliobj_stem_b,obj_stem_datestamp,shape,[],[],flag_pixel_correlated);
%stem_par object creation
obj_stem_par=stem_par(obj_stem_data,'exponential',flag_pixel_correlated,flag_time_diagonal);
%stem_model object creation
obj_stem_model=stem_model(obj_stem_data,obj_stem_par);
clear ground
clear remote

%Data transform
obj_stem_model.stem_data.log_transform;
obj_stem_model.stem_data.standardize;

obj_stem_par.alpha_bp=[0.4 0.4 0.8 0.8]';
if flag_pixel_correlated
    obj_stem_par.theta_b=100;
    obj_stem_par.v_b=[1 0.6;0.6 1];
else
    obj_stem_par.theta_b=[100 100];
    obj_stem_par.v_b=eye(2);
end

obj_stem_par.beta=obj_stem_model.get_beta0();
obj_stem_par.alpha_p=[0.6 0.6]';
obj_stem_par.theta_p=100;
obj_stem_par.v_p=[1 0.6;0.6 1];
obj_stem_par.sigma_eta=diag([0.2 0.2]);
obj_stem_par.G=diag([0.8 0.8]);
obj_stem_par.sigma_eps=diag([0.3 0.3 0.3 0.3]);

obj_stem_model.set_initial_values(obj_stem_par);

%Model estimation
exit_toll=0.002;
max_iterations=100;
obj_stem_EM_options=stem_EM_options(exit_toll,max_iterations);
obj_stem_model.EM_estimate(obj_stem_EM_options);
obj_stem_model.set_varcov;
obj_stem_model.set_logL;

obj_stem_model.print;        

subplot(2,1,1);
obj_stem_model.stem_data.plot('no2 remote','pixel',25);

size=obj_stem_model.stem_data.stem_gridlist_b.grid{1}.grid_size;
E_wb_y1=obj_stem_model.stem_EM_result.E_wb_y1(1:size(1)*size(2),25);
coordinate=obj_stem_model.stem_data.stem_gridlist_b.grid{1}.coordinate;
lat=coordinate(:,1);
lon=coordinate(:,2);
lat=reshape(lat,size);
lon=reshape(lon,size);
E_wb_y1=reshape(E_wb_y1,size);
subplot(2,1,2);
stem_misc.plot_map(lat,lon,E_wb_y1,obj_stem_model.stem_data.shape,'no2 remote estimated on 25-Jan-2009','longitude','latitude');