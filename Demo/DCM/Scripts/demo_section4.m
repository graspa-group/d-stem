clc
clearvars

addpath('../../../../D-STEM_dev/Src/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Data  building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the no2 ground level observations
load ../Data/no2_ground.mat
ground.Y{1} = no2_ground.data;
ground.Y_name{1} = 'no2 ground';
n1 = size(ground.Y{1}, 1);
T = size(ground.Y{1}, 2);

%load covariates for the NO2 monitoring stations
load ../Data/no2_ground_covariates.mat
ground.X_beta{1} = X;
ground.X_beta_name{1} = {'wind speed', 'elevation', 'sunday'};

ground.X_z{1} = ones(n1, 1);
ground.X_z_name{1} = {'constant'};

ground.X_p{1} = ground.X_beta{1}(:, 2, 1);
ground.X_p_name{1} = {'elevation'};

obj_stem_varset_p = stem_varset(ground.Y, ground.Y_name, [], [], ...
                                ground.X_beta, ground.X_beta_name, ...
                                ground.X_z, ground.X_z_name, ...
                                ground.X_p,ground.X_p_name);

%Coordinates
obj_stem_gridlist_p = stem_gridlist();
ground.coordinates{1} = [no2_ground.lat, no2_ground.lon];
obj_stem_grid = stem_grid(ground.coordinates{1}, 'deg', 'sparse', 'point');
obj_stem_gridlist_p.add(obj_stem_grid);
clear no2_ground

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Model building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj_stem_datestamp = stem_datestamp('01-01-2009 00:00', '31-12-2009 00:00', T);

%stem_data object creation
shape = shaperead('../Maps/worldmap');

obj_stem_modeltype = stem_modeltype('DCM');
obj_stem_data = stem_data(obj_stem_varset_p, obj_stem_gridlist_p, ...
                          [], [], obj_stem_datestamp, [], obj_stem_modeltype, shape);
%stem_par object creation
obj_stem_par_constraints=stem_par_constraints();
obj_stem_par = stem_par(obj_stem_data, 'exponential',obj_stem_par_constraints);
%stem_model object creation
obj_stem_model = stem_model(obj_stem_data, obj_stem_par);
clear ground

%Data transform
obj_stem_model.stem_data.log_transform;
obj_stem_model.stem_data.standardize;

%Starting values
obj_stem_par.beta = obj_stem_model.get_beta0();
obj_stem_par.theta_p = 100; %km
obj_stem_par.v_p = 1;
obj_stem_par.sigma_eta = 0.2;
obj_stem_par.G = 0.8;
obj_stem_par.sigma_eps = 0.3;
 
obj_stem_model.set_initial_values(obj_stem_par);

%Model estimation
exit_toll = 0.001;
max_iterations = 100;
obj_stem_EM_options = stem_EM_options(exit_toll, max_iterations);
obj_stem_EM_options.max_iterations=100;
obj_stem_model.EM_estimate(obj_stem_EM_options);
obj_stem_model.set_varcov;
obj_stem_model.set_logL;

%Kriging
load ../Data/kriging/krig_elevation_005;
load ../Data/kriging/krig_covariates
[LON,LAT] = meshgrid(krig_covariates.lon,krig_covariates.lat);
krig_coordinates = [LAT(:) LON(:)];
krig_mask = krig_elevation.data_mask;

obj_stem_krig_grid = stem_grid(krig_coordinates, 'deg', 'regular','pixel',size(LAT),'square',0.75,0.75);

obj_stem_krig_data = stem_krig_data(obj_stem_krig_grid,krig_covariates.data,krig_covariates.names,krig_mask);
obj_stem_krig = stem_krig(obj_stem_model,obj_stem_krig_data);

obj_stem_krig_options = stem_krig_options();
obj_stem_krig_options.block_size = 500;

obj_stem_krig_result = obj_stem_krig.kriging(obj_stem_krig_options);

obj_stem_krig_result{1}.plot(1)