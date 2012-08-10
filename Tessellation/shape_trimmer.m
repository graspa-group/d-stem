%refine a shape file given a bounding box
clc
clear all
shape=shaperead('c:\W.shp');
latmin_ref=54.5;
latmax_ref=62.5;
lonmin_ref=-9;
lonmax_ref=0;
toremove=[];
for i=1:length(shape)
    box=shape(i).BoundingBox;
    latmin=box(1,2);
    latmax=box(2,2);
    lonmin=box(1,1);
    lonmax=box(2,1);
    if (latmax<latmin_ref)||(latmin>latmax_ref)||(lonmax<lonmin_ref)||(lonmin>lonmax_ref)||(abs(latmin-latmax)<0.05)||(abs(lonmin-lonmax)<0.05)
        toremove=[toremove i];
    end
    i
end
shape(toremove)=[];