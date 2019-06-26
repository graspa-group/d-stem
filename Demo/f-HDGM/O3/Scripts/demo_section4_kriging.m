clc
clearvars

addpath('../../../../../D-STEM_dev/Src/'); %D-STEM
addpath('../../../../../D-STEM_dev/Src/fda/'); %fda library

%% Data loading

load ../Output/ozone_model
load ../Data/kriging/krig_coordinates;

%% Objects creation

[lon_mat,lat_mat] = meshgrid(krig_coordinates.lon,krig_coordinates.lat);
krig_coordinates_2D = [lat_mat(:) lon_mat(:)];

o_krig_grid = stem_grid(krig_coordinates_2D,'deg','regular','pixel',size(lat_mat),'square',0.05,0.05);
o_krig_data = stem_krig_data(o_krig_grid);

o_krig = stem_krig(ozone_model,o_krig_data);
o_krig_options = stem_krig_options();
o_krig_options.back_transform = 0;
o_krig_options.no_varcov = 0;
o_krig_options.block_size = 30; 

%% Kriging

o_krig_result = o_krig.kriging(o_krig_options);

save('../Output/kriging/ozone_kriging','o_krig_result','-v7.3');

%% Figure 4 of the paper
load ../Data/kriging/X_beta_t_100;

t=100; %day
h=10.5; %time within the day in hours
o_krig_result.surface_plot(h,t,X_beta_t_100);

%% Figure 5 of the paper
load ../Data/kriging/X_beta_h

h = 0:24;
lon = 116.25;
lat = 40.45;
t = 880;

o_krig_result.profile_plot(h,lon,lat,t,X_beta_h)