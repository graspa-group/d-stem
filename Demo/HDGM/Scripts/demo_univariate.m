clc
clearvars

%%

addpath('../../../../D-STEM_dev/Src/'); %D-STEM
addpath('../../../../D-STEM_dev/Src/fda/'); %fda library
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Data  building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

load ../Data/BJ_Data.mat

load ../Data/no2_ground.mat
ground.Y{1} = no2_ground.data;
ground.Y_name{1} = 'no2 ground';

load ../Data/no2_ground_covariates.mat
ground.X_beta{1} = X;
ground.X_beta_name{1} = {'wind speed', 'elevation', 'sunday'};

n1 = size(ground.Y{1}, 1);
ground.X_z{1} = ones(n1, 1);
ground.X_z_name{1} = {'constant'};

obj_stem_varset_p = stem_varset(ground.Y, ground.Y_name, [], [], ...
                                ground.X_beta, ground.X_beta_name, ... 
                                ground.X_z, ground.X_z_name);

 
ground.coordinates{1} = [no2_ground.lat, no2_ground.lon];                            
obj_stem_gridlist_p = stem_gridlist();
obj_stem_grid = stem_grid(ground.coordinates{1}, 'deg', 'sparse', 'point');
obj_stem_gridlist_p.add(obj_stem_grid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Model building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = size(ground.Y{1}, 2);
obj_stem_datestamp = stem_datestamp('01-01-2009 00:00','31-12-2009 00:00',T);

%stem_data object creation
S_xval=1:6:n1;
obj_stem_crossval=stem_crossval({'no2 ground'},{S_xval},0,{'point'});
%obj_stem_crossval=[];
shape = [];
obj_stem_modeltype = stem_modeltype('HDGM');
obj_stem_data = stem_data(obj_stem_varset_p, obj_stem_gridlist_p, ...
                          [], [], obj_stem_datestamp, obj_stem_crossval, obj_stem_modeltype, shape);

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
obj_stem_par.v_z = 20;
obj_stem_par.sigma_eta = 0.2;
obj_stem_par.G = 0.8;
obj_stem_par.sigma_eps = 0.3;

obj_stem_model.set_initial_values(obj_stem_par);

%Model estimation
exit_toll = 0.001;
max_iterations = 200;
obj_stem_EM_options = stem_EM_options(exit_toll, max_iterations);
obj_stem_model.EM_estimate(obj_stem_EM_options);
obj_stem_model.set_varcov;
obj_stem_model.set_logL;

%%
obj_stem_model.print
obj_stem_model.plot_xval

obj_stem_model_HDGM = obj_stem_model;
save('../Output/o_HDGM_uni.mat', 'obj_stem_model_HDGM', '-v7.3');