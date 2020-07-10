clc
clearvars

addpath(genpath('../../../../Src/')); %path of D-STEM source code

%% Data loading

%The temperature profiles data set is loaded. The data set is a MATLAB table
load ../Data/RAOB.mat

%% Objects creation

%model type
o_modeltype = stem_modeltype('f-HDGM');

%basis functions
spline_order = 2;
rng_spline = [50,925]; %range of the vertical temperature profiles
knots_number = 5; %number of spline knots
knots = linspace(rng_spline(1),rng_spline(2),knots_number); %knot locations (equispaced)

%A MATLAB structure is used to store spline information for all model terms
input_fda.spline_type = 'Bspline';
input_fda.spline_order = spline_order;
input_fda.spline_knots = knots;
input_fda.spline_range = rng_spline;
%The structure is passed as input of stem_fda constructor
o_fda = stem_fda(input_fda);

%A MATLAB structure is used to store the minimal information needed to
%create an object of class stem_data
input_data.data_table = RAOB;
input_data.stem_modeltype = o_modeltype;
input_data.stem_fda = o_fda;
%The structure is passed as input of stem_data constructor
o_data=stem_data(input_data);

%validation information
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

n = o_data.stem_varset_p.dim(1); %the number of spatial sites
S_val = sort(randperm(n,round(n*0.3))); %30% of the sites are used for validation
o_validation = stem_validation('Temperature',S_val);

%validation object is passed to o_data
o_data.stem_validation=o_validation;

%An object of class stem_par is created. The stem_data object and the 
%spatial correlation type are needed to shape the stem_par object internal 
%structure
o_par = stem_par(o_data,'exponential');

%Data and model parameters define a stem_model object
o_model = stem_model(o_data,o_par);

%Partitioning is adopted to speed up model estimation
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1000));
k = 5; %number of partitions
trials = 100; %number of trials of the modified k-means algorithm
lambda = 5000; %penalization parameter. lambda high favours partitions of similar sizes
partitions = o_data.kmeans_partitioning(k,trials,lambda); %the output is the dimension of partitions

%Initial values for model parameters are provided
n_basis = o_fda.get_basis_number;

o_par.beta = o_model.get_beta0();
o_par.theta_z = ones(1,n_basis.z)*3.6;
o_par.v_z = eye(n_basis.z);
o_par.G = diag(ones(n_basis.z,1)*0.8);
o_par.sigma_eps = ones(n_basis.sigma,1);

o_model.set_initial_values(o_par);

%Figure 6 of the paper
lat = 0;
lon = 0;
o_model.plot_profile(lat,lon);

%% Model estimation

%EM parameters
o_EM_options = stem_EM_options();
o_EM_options.exit_tol_par = 0.0001; %exit tolerance on parameter vector changes
o_EM_options.exit_tol_loglike = 0.0001; %exit tolerance on log-likelihood changes
o_EM_options.max_iterations = 100; %maximum number of EM iterations
o_EM_options.workers = 2; %number of workers (CPU-cores) for parallel computing
o_EM_options.partitions = partitions;

%EM estimation
o_model.EM_estimate(o_EM_options);

%Variance-covariance matrix evaluation
delta = 0.05; %parameter for approximating the variance-covariance matrix
o_model.set_varcov(delta);

%% Print and plots

%model estimation results
o_model.print;

%Figures 7 and 8 of the paper
vertical = 1;
o_model.plot_validation(vertical);

%% Model saving
if ~exist('../Output/','dir')
    mkdir('../Output/')
end
save('../Output/roab_model.mat', 'o_model', '-v7.3');

%% Kriging
           
%kriging grid
step_deg = 3; %grid step with respect to latitude and longitude
lat = 84-step_deg/2:-step_deg:-84+step_deg/2;
lon = -180+step_deg/2:step_deg:180-step_deg/2;
[lon_mat,lat_mat] = meshgrid(lon,lat);
krig_coordinates = [lat_mat(:) lon_mat(:)];

o_krig_grid = stem_grid(krig_coordinates, 'deg', 'regular','pixel',size(lat_mat),'square',step_deg,step_deg);
o_krig_data = stem_krig_data(o_krig_grid);
o_krig = stem_krig(o_model,o_krig_data);

o_krig_options = stem_krig_options();
o_krig_options.block_size = 150; %number of sites per block
o_krig_options.nn_size = 10; %number of nearest-neighbour sites
o_krig_options.workers = 2; %number of workers (CPU-cores) for parallel computing

o_krig_result = o_krig.kriging(o_krig_options);

%Figure 9 and 10 of the paper
h = 875.3; %hPa
t = 12; %time point
o_krig_result.surface_plot(h,t);

%% Kriging result saving

save('../Output/roab_kriging_result.mat', 'o_krig_result', '-v7.3');