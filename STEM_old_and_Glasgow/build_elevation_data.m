function elevation_data = build_elevation_data

elevation=load('-ascii','../Data/Elevation/elevation.ascii');
long_X=reshape(elevation(:,1),1239,426);
lat_X=reshape(elevation(:,2),1239,426);
data_X=reshape(elevation(:,3),1239,426);
data_X(data_X==-9999)=0;


data_X=reshape(data_X(3:3:1239,3:3:426),413,142);
data_X=fliplr(data_X);
data_X=rot90(data_X);
long_X=reshape(long_X(3:3:1239,3:3:426),413*142,1);
lat_X=reshape(lat_X(3:3:1239,3:3:426),413*142,1);

elevation_data.elevation=data_X;
elevation_data.lat_max=max(lat_X);
elevation_data.lat_min=min(lat_X);
elevation_data.long_max=max(long_X);
elevation_data.long_min=min(long_X);
elevation_data.step=long_X(2)-long_X(1);

end