%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef stem_gridlist < handle
    properties
        grid={};
        tap=[];
    end
    
    properties (SetAccess=private)
        box=[];
    end
    
    methods
        function obj = stem_gridlist(tap)
            if nargin>=1
                if not(isempty(tap))
                    obj.tap=tap;
                end
            end
        end
        
        function coordinates = get_jointcoordinates(obj)
            coordinates=[];
            for i=1:length(obj.grid)
                coordinates=[coordinates;obj.grid{i}.coordinate];
            end
            if isempty(coordinates)
                warning('The stem_gridlist object does not contain any stem_grid');
            end
        end        
        
        function DistMat = get_distance_matrix(obj)
            % return the distance matrix given a cell vector of grids
            d = obj.get_jointcoordinates();
            if isempty(obj.tap)
                DistMat=zeros(size(d,1));
                for z=1:length(d)
                    DistMat(z,z+1:end)=distdim(distance(d(z,:),d(z+1:end,:)), obj.grid{1}.unit, 'km');
                end
                DistMat=DistMat+DistMat';
            else
                idx_r=[];
                idx_c=[];
                elements=[];
                for z=1:length(d)
                    %evaluate the distance between the z-th coordinate and
                    %the vector ahead (z included)
                    dist_vec=distdim(distance(d(z,:),d(z:end,:)), obj.grid{1}.unit, 'km');
                    %IMPORTANT! the distances equal to zero are setted to
                    %eps so they are not confused with the zero generated
                    %by tapering
                    dist_vec(dist_vec==0)=eps;
                    
                    L=dist_vec<=obj.tap;
                    idx_r=[idx_r;ones(sum(L),1)*z];
                    idx_c=[idx_c;find(L)+z-1];
                    elements=[elements;dist_vec(L)];
                    %traspose
                    idx_c=[idx_c;ones(sum(L),1)*z];
                    idx_r=[idx_r;find(L)+z-1];
                    elements=[elements;dist_vec(L)];
                end
                DistMat=sparse(idx_r,idx_c,elements);
            end
        end
        
        function add(obj,stem_grid)
            if not(isa(stem_grid,'stem_grid'))
                error('The argument must be of class stem_grid');
            end
            if length(obj.grid)>1
                if not(strcmp(obj.grid{1}.unit,stem_grid.unit))
                    error('The grid unit differs from the unit of the grids already in stem_gridlist');
                end
            end 
            counter=length(obj.grid)+1;
            obj.grid{counter}=stem_grid;
            obj.updatebox();
        end
        
        function remove(obj,indices)
            if min(indices)<1
                error('The minimum index must be higher than zero');
            end
            if max(indices)>1
                error('The maximum index cannot be higher than the maximum number of grids');
            end
            obj.grid(indices)=[];
            obj.updatebox();
        end
        
        function set.grid(obj,grid)
            obj.grid=grid;
            obj.updatebox();
        end
        
        function set.tap(obj,tap)
            if tap<=0
                error('The tapering parameter must be > 0');
            end
            obj.tap=tap;
        end
    end
    
    methods (Access=private)
        function updatebox(obj)
            temp=[];
            for i=1:length(obj.grid)
                temp=[temp obj.grid{i}.box(1)];
            end
            obj.box(1)=min(temp);
            temp=[];
            for i=1:length(obj.grid)
                temp=[temp obj.grid{i}.box(2)];
            end
            obj.box(2)=max(temp);
            temp=[];
            for i=1:length(obj.grid)
                temp=[temp obj.grid{i}.box(3)];
            end
            obj.box(3)=min(temp);
            temp=[];
            for i=1:length(obj.grid)
                temp=[temp obj.grid{i}.box(4)];
            end
            obj.box(4)=max(temp);
        end
    end
end