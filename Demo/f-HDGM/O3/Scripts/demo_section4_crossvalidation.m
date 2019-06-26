clc
clearvars

addpath('../../../../../D-STEM_dev/Src/'); %D-STEM
addpath('../../../../../D-STEM_dev/Src/fda/'); %fda library

%% Data loading

%Ozone data set
load ../Data/Beijing_O3

%Beijing shapefile
shape = shaperead('../Maps/Beijing_adm1.shp');

%% Objects creation

%model type
o_modeltype = stem_modeltype('f-HDGM');

tmp=datenum('01-01-2015 00:00','dd-mm-yyyy HH:MM')+Beijing_O3.Time_step-1;
Beijing_O3.Time=datetime(datestr(tmp,'dd-mm-yyyy HH:MM'),'InputFormat','dd-MM-yyyy HH:mm');

[o_varset,o_gridlist,o_datestamp] = stem_misc.data_formatter(Beijing_O3);

%yaqiong
%date stamps
%T=o_varset.T;
%o_datestamp = stem_datestamp('01-01-2015 00:00','31-12-2017 00:00',T);

%cross-validation
S_xval = [1,7,10];
o_crossval = stem_crossval('O3',S_xval);

%basis functions
input_fda.spline_type = 'Fourier';
input_fda.spline_range = [0 24];
input_fda.spline_nbasis_z = 7;
input_fda.spline_nbasis_beta = 5;
input_fda.spline_nbasis_sigma = 5;

o_fda = stem_fda(input_fda);

input_data.stem_modeltype=o_modeltype;
input_data.stem_varset_p=o_varset;
input_data.stem_gridlist_p=o_gridlist;
input_data.stem_datestamp=o_datestamp;
input_data.stem_crossval=o_crossval;
input_data.stem_fda=o_fda;
input_data.shape=shape;

o_data = stem_data(input_data);

%spatial correlation function
sp_corr = 'exponential';

o_par = stem_par(o_data,sp_corr);

o_model = stem_model(o_data,o_par);

%Initial values for model parameters
n_basis=o_fda.get_basis_number;

o_par.beta = o_model.get_beta0();
o_par.sigma_eps = o_model.get_coe_log_sigma_eps0();
o_par.theta_z = ones(1,n_basis.z)*20;
o_par.G = diag(ones(n_basis.z,1)*0.5);
o_par.v_z = eye(n_basis.z)*10;

o_model.set_initial_values(o_par);

%% Model estimation

%EM parameters
exit_toll = 0.0001;
max_iterations = 200;
o_EM_options = stem_EM_options(exit_toll,max_iterations);

%EM estimation
o_model.EM_estimate(o_EM_options);

ozone_model_xval = o_model;
save('../Output/ozone_model_xval.mat', 'ozone_model_xval', '-v7.3');