s=shaperead('../Maps/CNTR_BN_03M_2010');

load ../Data/st_krig_result2.mat

grid_pixel=st_krig_result.stem_grid.coordinate;
lat=reshape(grid_pixel(:,1),128,184);
lon=reshape(grid_pixel(:,2),128,184);

data=st_krig_result.var_y_hat(:,:,30);

mapshow(s,'Color','black');
hold on
mapshow(lon,lat,data,'DisplayType','texturemap');