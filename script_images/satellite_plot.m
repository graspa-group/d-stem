clc
clear all
s=shaperead('../Maps/CNTR_BN_03M_2010');
load ../Data/st_model_20120603_135334.mat
load ../Data/grid_pixel.mat
lat=reshape(grid_pixel(:,1),128,184);
lon=reshape(grid_pixel(:,2),128,184);

data=st_model.stem_data.stem_varset_r.Y{1};
data=reshape(data,[128 184 35]);
temp=data(:,:,11);
temp(temp<-5)=NaN;

mapshow(s,'Color','black');
hold on
mapshow(lon,lat,temp,'DisplayType','texturemap');