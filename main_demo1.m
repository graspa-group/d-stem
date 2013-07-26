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

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Data  building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the no2 ground level observations
load ../Demo/Data/no2_ground.mat
ground.Y{1}=no2_ground.data;
ground.Y_name{1}='no2 ground';
n1=size(ground.Y{1},1);
T=size(ground.Y{1},2);

%X_beta
%load the covariates for the NO2 monitoring stations
load ../Demo/Data/no2_ground_covariates.mat
ground.X_beta{1}=X;
ground.X_beta_name{1}={'wind speed','elevation','sunday'};

ground.X_z{1}=ones(n1,1);
ground.X_z_name{1}={'constant'};

ground.X_p{1}=ground.X_beta{1}(:,2,1);
ground.X_p_name{1}={'elevation'};

obj_stem_varset_p=stem_varset(ground.Y,ground.Y_name,[],[],ground.X_beta,ground.X_beta_name,ground.X_z,ground.X_z_name,ground.X_p,ground.X_p_name);

%coordinates
obj_stem_gridlist_p=stem_gridlist();
ground.coordinates=[no2_ground.lat,no2_ground.lon];
obj_stem_grid=stem_grid(ground.coordinates,'deg','sparse','point');
obj_stem_gridlist_p.add(obj_stem_grid);
clear no2_ground

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Model building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj_stem_datestamp=stem_datestamp('01-01-2009 00:00','31-12-2009 00:00',T);

%stem_data object creation
shape=shaperead('../Demo/Maps/worldmap');
obj_stem_data=stem_data(obj_stem_varset_p,obj_stem_gridlist_p,[],[],obj_stem_datestamp,shape);
%stem_par object creation
obj_stem_par=stem_par(obj_stem_data,'exponential');
%stem_model object creation
obj_stem_model=stem_model(obj_stem_data,obj_stem_par);
clear ground

%Data transform
obj_stem_model.stem_data.log_transform;
obj_stem_model.stem_data.standardize;

%obj_stem_par object initialization
obj_stem_par.beta=obj_stem_model.get_beta0();
obj_stem_par.alpha_p=0.6;
obj_stem_par.theta_p=100; %km
obj_stem_par.v_p=1;
obj_stem_par.sigma_eta=0.2;
obj_stem_par.G=0.8;
obj_stem_par.sigma_eps=0.3;
 
obj_stem_model.set_initial_values(obj_stem_par);

%Model estimation
exit_toll=0.002;
max_iterations=100;
obj_stem_EM_options=stem_EM_options(exit_toll,max_iterations);
obj_stem_model.EM_estimate(obj_stem_EM_options);
obj_stem_model.set_varcov;
obj_stem_model.set_logL;

load ../Demo/Data/kriging/krig_elevation_005;
krig_coordinates=[krig_elevation.lat(:),krig_elevation.lon(:)];
krig_mask=krig_elevation.data_mask(:);
%kriging
obj_stem_krig=stem_krig(obj_stem_model);
obj_stem_krig_grid=stem_grid(krig_coordinates,'deg','regular','pixel',[80,170],'square',0.05,0.05);
back_transform=1;
no_varcov=0;
block_size=1000;
X_krig='../Demo/Data/kriging/blocks';
obj_stem_krig_result=obj_stem_krig.kriging('no2 ground',obj_stem_krig_grid,block_size,krig_mask,X_krig,back_transform,no_varcov);    

obj_stem_model.print;    
obj_stem_model.stem_EM_result.stem_kalmansmoother_result.plot;
obj_stem_krig_result.plot(100); %April 10, 2009
