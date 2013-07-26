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

flag_beta_ground=1;     %1: the model includes the term X_beta(s,t)*beta; 0: X_beta(s,t)*beta is not included
flag_time_ground=0;     %1: the model includes the term X_z(s,t)*z(t); 0: X_z(s,t)*z(t) is not included
flag_w_ground=1;        %1: the model includes the term X_g(s,t)*w(s,t);  0: X_g(s,t)*w(s,t) is not included

flag_time_diagonal=1;   %1: G and sigma_eta matrice are diagonal; 0: matrices are full
flag_tapering=1;        %1: tapering is enabled; 0: tapering is not enabled
flag_kriging=0;         %1: kriging is enabled; 0: kriging is not enabled
flag_parallel=0;

pathparallel='/home/finazzi/matNfs/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Data  building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the be ground level observations
load ../Data2013/be_ground/be_ground_background.mat
be_ground=out;
clear out
sd_g.Y{1}=be_ground.data;
sd_g.Y_name{1}='be ground';
n(1)=size(sd_g.Y{1},1);
T=size(sd_g.Y{1},2);
%load the co ground level observations
load ../Data2013/co_ground/co_ground_background.mat
co_ground=out;
clear out
sd_g.Y{2}=co_ground.data;
sd_g.Y_name{2}='co ground';
n(2)=size(sd_g.Y{2},1);
%load the no2 ground level observations
load ../Data2013/no2_ground/no2_ground_background.mat
no2_ground=out;
clear out
sd_g.Y{3}=no2_ground.data;
sd_g.Y_name{3}='no2 ground';
n(3)=size(sd_g.Y{3},1);
%load the o3 ground level observations
load ../Data2013/o3_ground/o3_ground_background.mat
o3_ground=out;
clear out
sd_g.Y{4}=o3_ground.data;
sd_g.Y_name{4}='o3 ground';
n(4)=size(sd_g.Y{4},1);
%load the pm10 ground level observations
load ../Data2013/pm10_ground/pm10_ground_background.mat
pm10_ground=out;
clear out
sd_g.Y{5}=pm10_ground.data;
sd_g.Y_name{5}='pm10 ground';
n(5)=size(sd_g.Y{5},1);
%load the pm2.5 ground level observations
load ../Data2013/pm25_ground/pm25_ground_background.mat
pm25_ground=out;
clear out
sd_g.Y{6}=pm25_ground.data;
sd_g.Y_name{6}='pm2.5 ground';
n(6)=size(sd_g.Y{6},1);
%load the so2 ground level observations
load ../Data2013/so2_ground/so2_ground_background.mat
so2_ground=out;
clear out
sd_g.Y{7}=so2_ground.data;
sd_g.Y_name{7}='so2 ground';
n(7)=size(sd_g.Y{7},1);

%no remote sensing data
sd_g.X_rg=[];
sd_g.X_rg_name=[];


if flag_beta_ground
    %X_beta
    load ../Data2013/be_ground/be_background_meteo_point.mat
    load ../Data2013/be_ground/be_background_elevation_point.mat
    load ../Data2013/be_ground/be_background_population_point.mat
    X_elevation_point=repmat(X_elevation_point,[1,1,T]);
    X_population_point(X_population_point==0)=0.02;
    X_population_point=log(X_population_point);
    X_population_point=repmat(X_population_point,[1,1,T]);
    %dummy variables for saturday and sunday
    X_dummy1=zeros(size(X_population_point,1),1,T);
    X_dummy2=zeros(size(X_population_point,1),1,T);
    X_dummy1(:,1,3:7:T)=1;
    X_dummy2(:,1,4:7:T)=1;
    X=cat(2,X_meteo_point,X_elevation_point,X_population_point,X_dummy1,X_dummy2);
    sd_g.X_beta{1}=X;
    sd_g.X_beta_name{1}={'wind speed','pressure','temperature','elevation','population','saturday','sunday'};

    load ../Data2013/co_ground/co_background_meteo_point.mat
    load ../Data2013/co_ground/co_background_elevation_point.mat
    load ../Data2013/co_ground/co_background_population_point.mat
    X_elevation_point=repmat(X_elevation_point,[1,1,T]);
    X_population_point(X_population_point==0)=0.02;
    X_population_point=log(X_population_point);
    X_population_point=repmat(X_population_point,[1,1,T]);
    %dummy variables for saturday and sunday
    X_dummy1=zeros(size(X_population_point,1),1,T);
    X_dummy2=zeros(size(X_population_point,1),1,T);
    X_dummy1(:,1,3:7:T)=1;
    X_dummy2(:,1,4:7:T)=1;
    X=cat(2,X_meteo_point,X_elevation_point,X_population_point,X_dummy1,X_dummy2);
    sd_g.X_beta{2}=X;
    sd_g.X_beta_name{2}={'wind speed','pressure','temperature','elevation','population','saturday','sunday'};
    
    load ../Data2013/no2_ground/no2_background_meteo_point.mat
    load ../Data2013/no2_ground/no2_background_elevation_point.mat
    load ../Data2013/no2_ground/no2_background_population_point.mat
    X_elevation_point=repmat(X_elevation_point,[1,1,T]);
    X_population_point(X_population_point==0)=0.02;
    X_population_point=log(X_population_point);
    X_population_point=repmat(X_population_point,[1,1,T]);
    %dummy variables for saturday and sunday
    X_dummy1=zeros(size(X_population_point,1),1,T);
    X_dummy2=zeros(size(X_population_point,1),1,T);
    X_dummy1(:,1,3:7:T)=1;
    X_dummy2(:,1,4:7:T)=1;
    X=cat(2,X_meteo_point,X_elevation_point,X_population_point,X_dummy1,X_dummy2);
    sd_g.X_beta{3}=X;
    sd_g.X_beta_name{3}={'wind speed','pressure','temperature','elevation','population','saturday','sunday'};
    
    load ../Data2013/o3_ground/o3_background_meteo_point.mat
    load ../Data2013/o3_ground/o3_background_elevation_point.mat
    load ../Data2013/o3_ground/o3_background_population_point.mat
    X_elevation_point=repmat(X_elevation_point,[1,1,T]);
    X_population_point(X_population_point==0)=0.02;
    X_population_point=log(X_population_point);
    X_population_point=repmat(X_population_point,[1,1,T]);
    %dummy variables for saturday and sunday
    X_dummy1=zeros(size(X_population_point,1),1,T);
    X_dummy2=zeros(size(X_population_point,1),1,T);
    X_dummy1(:,1,3:7:T)=1;
    X_dummy2(:,1,4:7:T)=1;
    X=cat(2,X_meteo_point,X_elevation_point,X_population_point,X_dummy1,X_dummy2);
    sd_g.X_beta{4}=X;
    sd_g.X_beta_name{4}={'wind speed','pressure','temperature','elevation','population','saturday','sunday'};
    
    load ../Data2013/pm10_ground/pm10_background_meteo_point.mat
    load ../Data2013/pm10_ground/pm10_background_elevation_point.mat
    load ../Data2013/pm10_ground/pm10_background_population_point.mat
    X_elevation_point=repmat(X_elevation_point,[1,1,T]);
    X_population_point(X_population_point==0)=0.02;
    X_population_point=log(X_population_point);
    X_population_point=repmat(X_population_point,[1,1,T]);
    %dummy variables for saturday and sunday
    X_dummy1=zeros(size(X_population_point,1),1,T);
    X_dummy2=zeros(size(X_population_point,1),1,T);
    X_dummy1(:,1,3:7:T)=1;
    X_dummy2(:,1,4:7:T)=1;
    X=cat(2,X_meteo_point,X_elevation_point,X_population_point,X_dummy1,X_dummy2);
    sd_g.X_beta{5}=X;
    sd_g.X_beta_name{5}={'wind speed','pressure','temperature','elevation','population','saturday','sunday'};
    
    load ../Data2013/pm25_ground/pm25_background_meteo_point.mat
    load ../Data2013/pm25_ground/pm25_background_elevation_point.mat
    load ../Data2013/pm25_ground/pm25_background_population_point.mat
    X_elevation_point=repmat(X_elevation_point,[1,1,T]);
    X_population_point(X_population_point==0)=0.02;
    X_population_point=log(X_population_point);
    X_population_point=repmat(X_population_point,[1,1,T]);
    %dummy variables for saturday and sunday
    X_dummy1=zeros(size(X_population_point,1),1,T);
    X_dummy2=zeros(size(X_population_point,1),1,T);
    X_dummy1(:,1,3:7:T)=1;
    X_dummy2(:,1,4:7:T)=1;
    X=cat(2,X_meteo_point,X_elevation_point,X_population_point,X_dummy1,X_dummy2);
    sd_g.X_beta{6}=X;
    sd_g.X_beta_name{6}={'wind speed','pressure','temperature','elevation','population','saturday','sunday'};
    
    load ../Data2013/so2_ground/so2_background_meteo_point.mat
    load ../Data2013/so2_ground/so2_background_elevation_point.mat
    load ../Data2013/so2_ground/so2_background_population_point.mat
    X_elevation_point=repmat(X_elevation_point,[1,1,T]);
    X_population_point(X_population_point==0)=0.02;
    X_population_point=log(X_population_point);
    X_population_point=repmat(X_population_point,[1,1,T]);
    %dummy variables for saturday and sunday
    X_dummy1=zeros(size(X_population_point,1),1,T);
    X_dummy2=zeros(size(X_population_point,1),1,T);
    X_dummy1(:,1,3:7:T)=1;
    X_dummy2(:,1,4:7:T)=1;
    X=cat(2,X_meteo_point,X_elevation_point,X_population_point,X_dummy1,X_dummy2);
    sd_g.X_beta{7}=X;
    sd_g.X_beta_name{7}={'wind speed','pressure','temperature','elevation','population','saturday','sunday'};
    
    clear X
    clear X_dummy1
    clear X_dummy2
    clear X_elevation_point
    clear X_meteo_point
    clear X_population_point
else
    %no covariates
    sd_g.X_beta=[];
    sd_g.X_beta_name=[];
end

%X_z
if flag_time_ground
    for i=1:length(sd_g.Y)
        sd_g.X_z{i}=ones(n(i),1,1);
        sd_g.X_z_name{i}={'constant'};
    end
else
    sd_g.X_z=[];
    sd_g.X_z_name=[];
end

if flag_w_ground
    for i=1:length(sd_g.Y)
        x1_temp=sd_g.X_beta{i}(:,4,1);
        x2_temp=sd_g.X_beta{i}(:,5,1);
        sd_g.X_g{i}=cat(4,x1_temp,x2_temp);
        sd_g.X_g_name{i}={'elevation','population'};
    end
    clear x1_temp
    clear x2_temp
else
    sd_g.X_g=[];
    sd_g.X_g_name=[];
end

st_varset_g=stem_varset(sd_g.Y,sd_g.Y_name,sd_g.X_rg,sd_g.X_rg_name,sd_g.X_beta,sd_g.X_beta_name,sd_g.X_z,sd_g.X_z_name,sd_g.X_g,sd_g.X_g_name);
clear sd_g
clear sd_r

if flag_tapering
    tapering_g=200; %km
    tapering_r=200; %km
else
    tapering_g=[];
    tapering_r=[];
end

st_gridlist_g=stem_gridlist(tapering_g);
st_grid=stem_grid([be_ground.lat,be_ground.lon],'deg','sparse','point');
st_gridlist_g.add(st_grid);
st_grid=stem_grid([co_ground.lat,co_ground.lon],'deg','sparse','point');
st_gridlist_g.add(st_grid);
st_grid=stem_grid([no2_ground.lat,no2_ground.lon],'deg','sparse','point');
st_gridlist_g.add(st_grid);
st_grid=stem_grid([o3_ground.lat,o3_ground.lon],'deg','sparse','point');
st_gridlist_g.add(st_grid);
st_grid=stem_grid([pm10_ground.lat,pm10_ground.lon],'deg','sparse','point');
st_gridlist_g.add(st_grid);
st_grid=stem_grid([pm25_ground.lat,pm25_ground.lon],'deg','sparse','point');
st_gridlist_g.add(st_grid);
st_grid=stem_grid([so2_ground.lat,so2_ground.lon],'deg','sparse','point');
st_gridlist_g.add(st_grid);

clear be_ground
clear co_ground
clear no2_ground
clear o3_ground
clear pm10_ground
clear pm25_ground
clear so2_ground

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Model building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

st_datestamp=stem_datestamp('01-01-2009','31-12-2011',T);
%crossval
% if flag_crossval
%     indices=1:2:n1;
%     st_crossval=stem_crossval('no2 ground','point',indices);
% else
%     st_crossval=[];
% end

%stem_data object creation
st_data=stem_data(st_varset_g,st_gridlist_g,[],[],st_datestamp);
%stem_par object creation
st_par=stem_par(st_data,'exponential',[],flag_time_diagonal);
%stem_model object creation
st_model=stem_model(st_data,st_par);
st_model.set_system_size();

st_model.stem_data.time_crop(1:2);
%Data transform
st_model.stem_data.log_transform;
st_model.stem_data.standardize;

%st_par object initialization
if flag_beta_ground
    %st_par.beta=st_model.get_beta0();
    st_par.beta=[-0.22469   0.0039199    -0.36778   -0.093424    0.075568   -0.012365   -0.035818    -0.22922    -0.03354    -0.37148   -0.035384     0.10002    -0.02021   -0.042864    -0.24189   -0.048992    -0.25627   -0.080688     0.32053   -0.073611     -0.1423     0.23284    0.019088     0.45239     0.14007    -0.11393     0.03966    0.066022    -0.26265    0.078944    -0.14133     0.01239     0.12167   -0.023468   -0.056904    -0.29643     0.12848    -0.22467     0.10089     0.13876    -0.00933   -0.031729    -0.11379   -0.064869    -0.14249    0.026245     0.12173   -0.015284   -0.025182]';
end
if flag_w_ground
    st_par.alpha_p=[0.39913     0.43455     0.46606      0.5465     0.53752     0.54563     0.47474;     0.34398     0.35515     0.44259     0.36287     0.38852     0.36926     0.40006]';
    st_par.theta_p=[119.8098      103.5975]';
    v_p(:,:,1)=[1.0000    0.0649    0.0916   -0.0958    0.0911    0.1088    0.0450;
                0.0649    1.0000    0.1305   -0.1142    0.1740    0.1563    0.1095;
                0.0916    0.1305    1.0000   -0.2124    0.1987    0.2088    0.0936;
               -0.0958   -0.1142   -0.2124    1.0000   -0.1109   -0.1237   -0.0196;
                0.0911    0.1740    0.1987   -0.1109    1.0000    0.3633    0.1694;
                0.1088    0.1563    0.2088   -0.1237    0.3633    1.0000    0.1289;
                0.0450    0.1095    0.0936   -0.0196    0.1694    0.1289    1.0000];

    v_p(:,:,2)=[1.0000    0.0658    0.0700   -0.0737    0.0659    0.0675    0.0507;
                0.0658    1.0000    0.0919   -0.0881    0.1222    0.1038    0.0872;
                0.0700    0.0919    1.0000   -0.1618    0.1548    0.1640    0.1270;
               -0.0737   -0.0881   -0.1618    1.0000   -0.0909   -0.0990   -0.0317;
                0.0659    0.1222    0.1548   -0.0909    1.0000    0.2640    0.1362;
                0.0675    0.1038    0.1640   -0.0990    0.2640    1.0000    0.1027;
                0.0507    0.0872    0.1270   -0.0317    0.1362    0.1027    1.0000];

    st_par.v_p=v_p;
end

if flag_time_ground
    st_par.sigma_eta=diag(repmat(0.2,7,1));
    st_par.G=diag(repmat(0.8,7,1));
end

st_par.sigma_eps=diag([0.61768     0.60308     0.32316     0.34159     0.45583     0.45869      0.6316]);
st_model.set_initial_values(st_par);

%Model estimation
exit_toll=0.001;
max_iterations=100;
st_EM_options=stem_EM_options(exit_toll,max_iterations,'single');
if flag_parallel
 st_EM_options.pathparallel=pathparallel;
end
st_model.EM_estimate(st_EM_options);


if flag_kriging
    load ../Demo/Data/kriging/krig_elevation_005;
    krig_coordinates=[krig_elevation.lat(:),krig_elevation.lon(:)];
    krig_mask=krig_elevation.data_mask(:);
    %kriging
    st_krig=stem_krig(st_model);
    st_krig_grid=stem_grid(krig_coordinates,'deg','regular','pixel',[80,170],'square',0.05,0.05);
    back_transform=1;
    no_varcov=1;
    block_size=1000;
    X_krig='../Demo/Data/kriging/blocks';
    st_krig_result=st_krig.kriging('no2 ground',st_krig_grid,block_size,krig_mask,X_krig,back_transform,no_varcov);    
end

