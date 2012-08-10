classdef stem_tessellation < handle
    methods
        function obj = stem_tessellation()
            if nargin>0
                error('No input arguments are necessary.');
            end
        end
    end
    
    methods (Static)
        function [st_gridlist,dim,side_length] = get_knots(levels,box,rot)
            %INPUT
            %levels: Rx1 vector of tessellation resolutions
            %box: 4x1 vector [minlat, minlon, maxlat, maxlon]
            %rot: Rx2 matrix. Each 2x1 row contains the grid rotation along latitude and longitude respectively (with respect to the centre of the sphere)
            %
            %OUTPUT
            %st_gridlist: stem_gridlist object. Each grid is a knots grid.
            
            if nargin<1
                error('The first input arguments must be provided');
            end
            for i=1:length(levels)
                if levels(i)<=0
                    error('The elements of levels must be > 0');
                end
            end
            if (nargin<2)||(isempty(box))
                box=[-90,-180,90,180];
            end
            if nargin<3
                rot=zeros(length(levels),2);
            end
            if not(length(box)==4)
                error('The input box must be a 4 elements vector');
            end
            if not((size(rot,2)==2)&&(size(rot,1)==length(levels)))
                error('The input rot must be a levels x 2 matrix');
            end
            disp('Knots computation...');
            st_gridlist=stem_gridlist();
            for i=1:length(levels)
                fv = sphere_tri('ico',levels(i),1);
                [lon,lat] = cart2sph(fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3));
                lat=rad2deg(lat)+rot(i,1);
                lon=rad2deg(lon)+rot(i,2);
                v1(1)=lat(fv.faces(1,2));
                v1(2)=lon(fv.faces(1,2));
                v2(1)=lat(fv.faces(1,3));
                v2(2)=lon(fv.faces(1,3));
                side_length(i)=distdim(distance(v1,v2), 'deg', 'km');
                lat(lat>90)=lat(lat>90)-180;
                lat(lat<-90)=lat(lat<-90)+180;
                lon(lon>180)=lon(lon>180)-360;
                lon(lon<-180)=lon(lon<-180)+360;
                L=(lat>=box(1))&(lon>=box(2))&(lat<=box(3))&(lon<=box(4));
                lat=lat(L);
                lon=lon(L);
                st_grid=stem_grid([lat lon],'sparse','point',[]);
                st_gridlist.add(st_grid);
                dim(i)=length(lat);
            end
            disp('Knots computationd ended.');
        end
    end
    
end