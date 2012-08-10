classdef stem_base_functions < handle
    methods
        function obj = stem_base_functions()
            if nargin>0
                error('No input arguments are necessary.');
            end
        end
    end
    
    methods (Static)
        function K = get_weight_matrix(st_gridlist,st_gridlist_knot,type,func_par)
            if nargin<3
                error('All the input arguments are required');
            end
            if not(strcmp(type,'bisquare'))
                error('Only the bisquare function type is supported');
            end
            if strcmp(type,'bisquare')
                if not(length(func_par)==length(st_gridlist_knot.grid))
                    error('func_par must be a vector with length equal to the number of knot grids');
                end
            end
            dim_r=[];
            for i=1:length(st_gridlist.grid)
                dim_r(i)=size(st_gridlist.grid{i}.coordinate,1);
            end
            dim_c=[];
            for i=1:length(st_gridlist_knot.grid)
                dim_c(i)=size(st_gridlist_knot.grid{i}.coordinate,1);
            end
            K=zeros(sum(dim_r),sum(dim_c));
            blocks_c=[0 cumsum(dim_c)];
            d = stem_gridlist.get_jointcoordinates(st_gridlist);
            for i=1:length(st_gridlist_knot.grid)
                for k=1:length(st_gridlist_knot.grid{i}.coordinate)
                    D=distdim(distance(d,st_gridlist_knot.grid{i}.coordinate(k,:)), 'deg', 'km');
                    L=D<func_par(i);
                    D=(1-(D./func_par(i)).^2).^2.*L;
                    K(:,k+blocks_c(i))=D;
                end
            end
        end
    end
        
end