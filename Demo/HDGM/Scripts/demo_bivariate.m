clc
clearvars

addpath('../../../Src/'); %D-STEM

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Data  building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the no2 ground level observations

load ../Data/no2_ground.mat
ground.Y{1} = no2_ground.data;
ground.Y_name{1} = 'no2 ground';
n1 = size(ground.Y{1}, 1);
T = size(ground.Y{1}, 2);

%load the pm25 ground level observations
load ../Data/pm25_ground.mat
ground.Y{2} = pm25_ground.data;
ground.Y_name{2} = 'pm2.5 ground';
n2 = size(ground.Y{2}, 1);

%load covariates for the NO2 monitoring stations
load ../Data/no2_ground_covariates.mat
ground.X_beta{1} = X;
ground.X_beta_name{1} = {'wind speed', 'elevation', 'sunday'};

%load covariates for the PM2.5 monitoring stations

load ../Data/pm25_ground_covariates.mat
ground.X_beta{2} = X;
ground.X_beta_name{2} = {'wind speed', 'elevation', 'sunday'};


%X_z
ground.X_z{1} = ones(n1, 1);
ground.X_z_name{1} = {'constant'};
ground.X_z{2} = ones(n2, 1);
ground.X_z_name{2} = {'constant'};

obj_stem_varset_p = stem_varset(ground.Y, ground.Y_name, [], [], ...
                                ground.X_beta, ground.X_beta_name, ... 
                                ground.X_z, ground.X_z_name);
                            
obj_stem_gridlist_p = stem_gridlist();
ground.coordinates{1} = [no2_ground.lat, no2_ground.lon];
ground.coordinates{2} = [pm25_ground.lat, pm25_ground.lon];
obj_stem_grid1 = stem_grid(ground.coordinates{1}, 'deg', 'sparse', 'point');
obj_stem_grid2 = stem_grid(ground.coordinates{2}, 'deg', 'sparse', 'point');
obj_stem_gridlist_p.add(obj_stem_grid1);
obj_stem_gridlist_p.add(obj_stem_grid2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Model building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj_stem_datestamp = stem_datestamp('01-01-2009 00:00','31-12-2009 00:00',T);

%stem_data object creation
shape = [];
%obj_stem_validation=[];
S_val1=1:5:n1;
S_val2=1:5:n2;
obj_stem_validation=stem_validation({'no2 ground','pm2.5 ground'},{S_val1,S_val2},0,{'point','point'});

obj_stem_modeltype = stem_modeltype('HDGM');
obj_stem_data = stem_data(obj_stem_varset_p, obj_stem_gridlist_p, ...
                          [], [], obj_stem_datestamp, obj_stem_validation, obj_stem_modeltype, shape);

%stem_par object creation
obj_stem_par_constraints=stem_par_constraints();
obj_stem_par_constraints.time_diagonal=0;
obj_stem_par = stem_par(obj_stem_data, 'exponential',obj_stem_par_constraints);
%stem_model object creation
obj_stem_model = stem_model(obj_stem_data, obj_stem_par);
clear ground

%Data transform
obj_stem_model.stem_data.log_transform;
obj_stem_model.stem_data.standardize;

%Starting values
obj_stem_par.beta = obj_stem_model.get_beta0();
obj_stem_par.theta_z = 100;
obj_stem_par.v_z = [1 0.6;0.6 1];
obj_stem_par.sigma_eta = diag([0.2 0.2]);
obj_stem_par.G = diag([0.8 0.8]);
obj_stem_par.sigma_eps = diag([0.3 0.3]);

obj_stem_model.set_initial_values(obj_stem_par);

%Model estimation
exit_toll = 0.001;
max_iterations = 200;
obj_stem_EM_options = stem_EM_options(exit_toll, max_iterations);
obj_stem_model.EM_estimate(obj_stem_EM_options);
obj_stem_model.set_varcov;
obj_stem_model.set_logL;


obj_stem_model.print;      
% validation result for no2
obj_stem_model.plot_val(0,'no2 ground');
% validation result for pm2.5
obj_stem_model.plot_val(0,'pm2.5 ground');