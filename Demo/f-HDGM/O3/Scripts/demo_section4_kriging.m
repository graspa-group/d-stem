clc
clearvars

addpath(genpath('../../../../Src/')); %path of D-STEM source code

%% Data loading

%load the previously estimated model
if exist('../Output/ozone_model','file')==2
    error('Run demo_section4_model_estimate.m first');
else
    load ../Output/ozone_model
end

%% Kriging

%kriging grid
step_deg = 0.05; %grid step with respect to latitude and longitude
lat = 39.4:step_deg:41.1;
lon = 115.4:step_deg:117.5;
[lon_mat,lat_mat] = meshgrid(lon,lat);
krig_coordinates = [lat_mat(:) lon_mat(:)];

o_krig_grid = stem_grid(krig_coordinates,'deg','regular','pixel',size(lat_mat),'square',0.05,0.05);
o_krig_data = stem_krig_data(o_krig_grid);

%Beijing shapefile. Will be used for plotting maps
o_model.stem_data.shape = shaperead('../Maps/Beijing_adm1.shp');
o_krig = stem_krig(o_model,o_krig_data);

o_krig_options = stem_krig_options();
o_krig_options.block_size = 30; %number of sites per block
o_krig_options.back_transform = 0; %y is not back-transformed to the original scale

o_krig_result = o_krig.kriging(o_krig_options);

%% Plots

%Figure 4 of the paper

%Covariates loading. Covariates are specific to t=100 and h=10.5
load ../Data/kriging/X_beta_t_100;

t = 100; %day
h = 10.5; %time within the day in hours
o_krig_result.surface_plot(h,t,X_beta_t_100);

%Figure 5 of the paper

%Covariates loading. Covariates are specific to h=0:24, lon=116.25, lat=40.45 and t=880
load ../Data/kriging/X_beta_h

h = 0:24; %covariates domain
lon = 116.25; %longitude
lat = 40.45; %latitude
t = 880; %time point

o_krig_result.profile_plot(h,lon,lat,t,X_beta_h)

%% Kriging result saving

save('../Output/ozone_kriging','o_krig_result','-v7.3');