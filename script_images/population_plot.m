s=shaperead('../../Maps/CNTR_BN_03M_2010');
load ../../Data/population.mat
pop2(pop2==0)=NaN;
load ../../Data/grid_pixel.mat
lat=reshape(grid_pixel(:,1),128,184);
lon=reshape(grid_pixel(:,2),128,184);
mapshow(s,'Color','black');
hold on
mapshow(lon,lat,log10(pop2),'DisplayType','texturemap');
