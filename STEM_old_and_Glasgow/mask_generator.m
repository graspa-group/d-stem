%script for generating mask from shapefile
clc
clear all
load ../Data/krig_lat_scotland2009
load ../Data/krig_lon_scotland2009
shape = shaperead('../Maps/Scotland/scotland_only');
shape.X(shape.Y>56)=[];
shape.Y(shape.Y>56)=[];

%shape.X=shape.X(120433:229355); %scotland inland
%shape.Y=shape.Y(120433:229355);

[lat,lon]=meshgrid(krig_lat,krig_lon);
L=lat<56;
lat=lat(L);
lon=lon(L);
IN = inpolygon(lon,lat,shape.X,shape.Y);

IN=reshape(IN,961,180);
IN=double(IN);
IN(IN==0)=NaN;

load ../Data/krig_mask_scotland2009
mask(:,1:180)=mask(:,1:180).*IN;