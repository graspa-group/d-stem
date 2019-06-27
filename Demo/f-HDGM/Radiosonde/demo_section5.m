clc
clearvars

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

addpath(genpath('../../../../Src/')); %D-STEM

%% Data loading

%Temperature profiles data set
load Data/RAOB.mat

%% Objects creation

%model type
o_modeltype=stem_modeltype('f-HDGM');

[o_varset,o_gridlist,o_datestamp]=stem_misc.data_formatter(RAOB);

%cross-validation
n=o_varset.dim(1);
S_xval=sort(randperm(n,round(n*0.3))); % 30% of the sites are used for cross-validation
o_crossval=stem_crossval('Temperature',S_xval);

%basis functions
spline_order=2;
rng_spline=[50,925];
knots_number=5;
knots=linspace(rng_spline(1),rng_spline(2),knots_number);

input_fda.spline_type='Bspline';
input_fda.spline_order_beta=spline_order;
input_fda.spline_order_z=spline_order;
input_fda.spline_order_sigma=spline_order;
input_fda.spline_knots_beta=knots;
input_fda.spline_knots_z=knots;
input_fda.spline_knots_sigma=knots;
input_fda.spline_range=rng_spline;
o_fda=stem_fda(input_fda);

input_data.stem_modeltype=o_modeltype;
input_data.stem_varset_p=o_varset;
input_data.stem_gridlist_p=o_gridlist;
input_data.stem_datestamp=o_datestamp;
input_data.stem_crossval=o_crossval;
input_data.stem_fda=o_fda;

o_data=stem_data(input_data);

%spatial correlation function type
sp_corr = 'exponential';

o_par=stem_par(o_data,sp_corr);

o_model=stem_model(o_data,o_par);

%k-means for block tapering
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1000));
k=5; 
trials=100;
lambda=5000;
block_size=o_data.kmeans_globe(k,trials,lambda);

%Initial values for model parameters
n_basis=o_fda.get_basis_number;

o_par.beta = o_model.get_beta0();
o_par.theta_z=ones(1,n_basis.z)*400;
o_par.v_z=eye(n_basis.z);
o_par.G=diag(ones(n_basis.z,1)*0.8);
o_par.sigma_eps=ones(n_basis.sigma,1);

o_model.set_initial_values(o_par);

%Figure 6 of the paper

lat=0;
lon=0;
o_model.plot_profile(lat,lon);

%% Model estimation

%EM parameters
exit_toll=10^-4;
max_iterations=100;
o_EM_options=stem_EM_options(exit_toll,max_iterations);
o_EM_options.workers=4;
o_EM_options.block_tapering_block_size=block_size;

%EM estimation
o_model.EM_estimate(o_EM_options);

%Variance-covariance matrix evaluation
delta=0.05;
o_model.set_varcov(delta);

%Figures 7 and 8 of the paper
vertical=1;
o_model.plot_xval(vertical);

%% Kriging
           
%kriging grid
step_deg=6; %degrees
lat=84-step_deg/2:-step_deg:-84+step_deg/2;
lon=-180+step_deg/2:step_deg:180-step_deg/2;
[lon_mat,lat_mat]=meshgrid(lon,lat);

krig_coordinates=[lat_mat(:) lon_mat(:)];
o_krig_grid = stem_grid(krig_coordinates, 'deg', 'regular','pixel',size(lat_mat),'square',step_deg,step_deg);
o_krig_data = stem_krig_data(o_krig_grid);
o_krig = stem_krig(o_model,o_krig_data);

o_krig_options = stem_krig_options();
o_krig_options.block_size = 150; %sites per block
o_krig_options.nn_size = 10; %nearest-neighbour sites
o_krig_options.workers = 4;

o_krig_result = o_krig.kriging(o_krig_options);

%Figure 9 and 10 of the paper
h=875.3; %hPa
t=12;
o_krig_result.surface_plot(h,t);