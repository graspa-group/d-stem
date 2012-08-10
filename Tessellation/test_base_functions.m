%Test bisquare
clc
clear all
load coast
lat_eu=lat;
long_eu=long;

st_gridlist=stem_gridlist();

%first data grid
lat1=[30:0.25:75];
lon1=[-25:0.25:45];
[LON1,LAT1]=meshgrid(lon1,lat1);
st_grid=stem_grid([LAT1(:),LON1(:)],'sparse','point',[]);
st_gridlist.add(st_grid);

%second data grid
lat2=[30:0.5:75];
lon2=[-25:0.5:45];
[LON2,LAT2]=meshgrid(lon2,lat2);
st_grid=stem_grid([LAT2(:),LON2(:)],'sparse','point',[]);
st_gridlist.add(st_grid);

%knots grid generation
resolution=[3,4,5];
box=[30,-25,75,45];
rotation=[0 0;3,3;2,2];
[st_gridlist_knot,dim,side_length] = stem_tessellation.get_knots(resolution,box,rotation);

%K matrix generation
bisquare_par=side_length*1.5;
K = stem_base_functions.get_weight_matrix(st_gridlist,st_gridlist_knot,'bisquare',bisquare_par);


%Following the recommendation of Cressie and Johannesson (2008) to set the 
%range w of the bisquare functions equal to 1.5 times the distance 
%of two adjacent centers at the same resolution.

%plot
k=K(1:181*281,655);
k=reshape(k,181,281);

ax = worldmap('Europe');
geoshow(ax, lat_eu, long_eu,'DisplayType', 'polygon', 'FaceColor', [1 1 1],'LineWidth',2);
geoshow(ax, st_gridlist_knot.grid{1}.coordinate(:,1),st_gridlist_knot.grid{1}.coordinate(:,2),'DisplayType','multipoint','Marker','o');
geoshow(ax, st_gridlist_knot.grid{2}.coordinate(:,1),st_gridlist_knot.grid{2}.coordinate(:,2),'DisplayType','multipoint','Marker','+');
geoshow(ax, st_gridlist_knot.grid{3}.coordinate(:,1),st_gridlist_knot.grid{3}.coordinate(:,2),'DisplayType','multipoint','Marker','.');
geoshow(LAT1,LON1,k,'DisplayType','texturemap')



