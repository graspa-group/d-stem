%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% D-STEM - Distributed Space Time Expecation Maximization              %
%%%                                                                      %
%%% Author: Francesco Finazzi                                            %
%%% E-mail: francesco.finazzi@unibg.it                                   %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Economics and Quantitative Methods %
%%% Author website: http://www.unibg.it/pers/?francesco.finazzi          %
%%% Code website: https://code.google.com/p/d-stem/                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This file is part of D-STEM.
% 
% D-STEM is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 2 of the License, or
% (at your option) any later version.
% 
% D-STEM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with D-STEM. If not, see <http://www.gnu.org/licenses/>.

classdef stem_gridlist < handle
    properties
        grid={};    %[stem_grid object] {qx1} cell-array of stem_grid objects (one for each variable)
        tap=[];     %[double >0]        (1x1) the tapering parameter. It is the maximum distance after which the spatial correlation is zero
    end
    
    properties (SetAccess=private)
        box=[];     %[double]           (4x1) the bounding box of the geographic area covered by all the grids [lat_min,lat_max,lon_min,lon_max]
    end
    
    methods
        function obj = stem_gridlist(tap)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %<tap>        - [double >0]             (1x1) the tapering parameter. It is the maximum distance after which the spatial correlation is zero. If it is empty the full distance matrix is evaluated
            %
            %OUTPUT
            %obj          - [stem_gridlist object]  (1x1)            
            if nargin>=1
                if not(isempty(tap))
                    obj.tap=tap;
                end
            end
        end
        
        function coordinates = get_jointcoordinates(obj)
            %DESCRIPTION: get the coordinates of all the spatial locations of all the variables (pixel or points)
            %
            %INPUT
            %obj            - [stem_gridlist object]    (1x1)
            %
            %OUTPUT
            %coordinates    - [double]                  (N_r|N_gx2)              
            coordinates=[];
            for i=1:length(obj.grid)
                coordinates=[coordinates;obj.grid{i}.coordinate];
            end
            if isempty(coordinates)
                warning('The stem_gridlist object does not contain any stem_grid');
            end
        end        
        
        function DistMat = get_distance_matrix(obj,type)
            %DESCRIPTION: get the distance matrix of all the variables
            %
            %INPUT
            %obj        - [stem_gridlist object]    (1x1)
            %<type>     - [boolean]                 (1x1) 1: also the cross-distances are evaluated (distances between different variables); 0: the distance matrix is block-diagonal with respect to the variables
            %
            %OUTPUT
            %DistMat    - [double]                  (N_rxN_r|N_gxN_g)  The distance matrix
            if nargin<2
                type=1;
            end
            
            if isempty(obj.tap)
                d = obj.get_jointcoordinates();
                DistMat=zeros(size(d,1));
                for z=1:length(d)
                    DistMat(z,z+1:end)=distdim(distance(d(z,:),d(z+1:end,:)), obj.grid{1}.unit, 'km');
                end
                DistMat=DistMat+DistMat';
            else
                if (type==1||length(obj.grid)==1)
                    d = obj.get_jointcoordinates();
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
                else
                    for i=1:length(obj.grid)
                        d = obj.grid{i}.coordinate;
                        evaluate=1;
                        if i>1
                            try
                                s=d-obj.grid{1}.coordinate;
                                s=s(:);
                                s=sum(abs(s));
                                if s==0
                                    evaluate=0;
                                end
                            catch
                                %if the difference cannot be evaluated as
                                %the size is different
                                evaluate=1;
                            end
                        end
                        if evaluate
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
                            Dist{i}=sparse(idx_r,idx_c,elements);
                        else
                            Dist{i}=Dist{1};
                        end
                    end
                    DistMat=blkdiag(Dist{1},Dist{2});
                    for i=3:length(obj.grid)
                        DistMat=blkdiag(DistMat,Dist{i});
                    end
                end
            end
        end
        
        function add(obj,stem_grid)
            %DESCRIPTION: add a stem_grid object to the stem_gridlist object
            %
            %INPUT
            %obj        - [stem_gridlist object]    (1x1)
            %stem_grid  - [stem_grid object]        (1x1) the stem_grid object to add
            %
            %OUTPUT
            %none: the stem_grid object is added to the list          
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
            %DESCRIPTION: remove a stem_grid object from the stem_gridlist object
            %
            %INPUT
            %obj        - [stem_gridlist object]    (1x1)
            %indices    - [integer >0]              (1x1) the index of the stem_grid object to remove
            %
            %OUTPUT
            %none: the stem_grid object is removed from the list               
            if min(indices)<1
                error('The minimum index must be higher than zero');
            end
            if max(indices)>1
                error('The maximum index cannot be higher than the maximum number of grids');
            end
            obj.grid(indices)=[];
            obj.updatebox();
        end
        
        %Class set methods
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
            %DESCRIPTION: update the box property of the object
            %
            %INPUT
            %obj        - [stem_gridlist object] (1x1)
            %
            %OUTPUT
            %none: the box property is updated            
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