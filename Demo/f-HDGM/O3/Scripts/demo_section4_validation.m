clc
clearvars

addpath(genpath('../../../../Src/')); %path of D-STEM source code

%% Data loading

%The ozone profiles data set is loaded. The data set is a MATLAB table
load ../Data/Beijing_O3

%% Objects creation

%model type
o_modeltype = stem_modeltype('f-HDGM');

%A MATLAB structure is used to store spline information for all model terms
input_fda.spline_type = 'Fourier';
input_fda.spline_range = [0 24];
input_fda.spline_nbasis_z = 7;
input_fda.spline_nbasis_beta = 5;
input_fda.spline_nbasis_sigma = 5;
%The structure is passed as input of stem_fda constructor
o_fda = stem_fda(input_fda);

%A MATLAB structure is used to store the minimal information needed to
%create an object of class stem_data
input_data.data_table = Beijing_O3;
input_data.stem_modeltype = o_modeltype;
input_data.stem_fda = o_fda;
input_data.shape = shape;
%The structure is passed as input of stem_data constructor
o_data = stem_data(input_data);

%validation
S_val = [1,7,10]; %indexes of the validation sites
o_data.stem_validation = stem_validation('O3',S_val);

%An object of class stem_par is created. The stem_data object and the 
%spatial correlation type are needed to shape the stem_par object internal 
%structure
o_par = stem_par(o_data,'exponential');

%Data and model parameters define a stem_model object
o_model = stem_model(o_data,o_par);

%Initial values for model parameters are provided
n_basis = o_fda.get_basis_number;

o_par.beta = o_model.get_beta0();
o_par.sigma_eps = o_model.get_coe_log_sigma_eps0();
o_par.theta_z = ones(1,n_basis.z)*0.18;
o_par.G = diag(ones(n_basis.z,1)*0.5);
o_par.v_z = eye(n_basis.z)*10;

o_model.set_initial_values(o_par);

%% Model estimation

%EM parameters
o_EM_options = stem_EM_options();
o_EM_options.exit_tol_par = 0.0001; %exit tolerance on parameter vector changes
o_EM_options.exit_tol_loglike = 0.0001; %exit tolerance on log-likelihood changes
o_EM_options.max_iterations = 200; %maximum number of EM iterations

%EM estimation
o_model.EM_estimate(o_EM_options);

%% Model saving

if ~exist('../Output/','dir')
    mkdir('../Output/')
end
save('../Output/ozone_model_val.mat', 'o_model', '-v7.3');