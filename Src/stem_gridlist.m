%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% D-STEM - Distributed Space Time Expecation Maximization              %
%%%                                                                      %
%%% Author: Francesco Finazzi                                            %
%%% E-mail: francesco.finazzi@unibg.it                                   %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Information and                    %
%%%              Production Engineering                                  %
%%% Author website: http://www.unibg.it/pers/?francesco.finazzi          %
%%%                                                                      %
%%% Author: Yaqiong Wang                                                 %
%%% E-mail: yaqiongwang@pku.edu.cn                                       %
%%% Affiliation: Peking University,                                      %
%%%              Guanghua school of management,                          %
%%%              Business Statistics and Econometrics                    %
%%%                                                                      %
%%% Author: Alessandro Fassò                                             %
%%% E-mail: alessandro.fasso@unibg.it                                    %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Information and                    %
%%%              Production Engineering                                  %
%%% Author website: http://www.unibg.it/pers/?alessandro.fasso           %
%%%                                                                      %
%%% Code website: https://github.com/graspa-group/d-stem                 %
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
    
    %PROPERTIES
    %Each class property or method property is defined as follows
    %
    %"Name"="Default value";    %["type"]    "dimension"     "description" 
    %
    %DIMENSION NOTATION
    %(1 x 1) is a scalar
    %(N x 1) is a Nx1 vector
    %(N x T) is a NxT matrix
    %(N x B x T) is a NxBxT array
    %{q} is a cell array of length q
    %{q}{p} is a cell array of length q, each cell is a cell array of length p
    %{q}(NxT) is a cell array of length q, each cell is a NxT matrix
    %
    %CONSTANTS
    %q    - the number of variables
    %N_p  - total number of point sites
    %N_b  - total number of pixel sites
    
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
            %tap          - [double >0]             (1x1) the tapering parameter. It is the maximum distance after which the spatial correlation is zero. If it is empty the full distance matrix is evaluated
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
                coordinates=cat(1,coordinates,obj.grid{i}.coordinate);
            end
            if isempty(coordinates)
                warning('The stem_gridlist object does not contain any stem_grid');
            end
        end        
        
        function DistMat = get_distance_matrix(obj,stem_modeltype,correlation_type,cross_type,idx_var)
            %DESCRIPTION: get the distance matrix of all the variables. If the model is of the emulator type, a distance matrix is given for each dimension of the inputs
            %
            %INPUT
            %obj                - [stem_gridlist object]    (1x1)
            %stem_modeltype     - [stem_modeltype object]   (1x1) object of class stem_modeltype
            %correlation_type   - [string]                  (1x1) the name of the correlation function (see the correlation_function method of the stem_misc class for valid names)
            %cross_type         - [boolean]                 (1x1) 1: also the cross-distances are evaluated (distances between different variables); 0: the distance matrix is block-diagonal with respect to the variables
            %idx_var            - [integer>0]               (1x1) the index of the variable in order to get the distance matrix of that variable only
            %
            %OUTPUT
            %DistMat            - [double]                  (N_pxN_p|N_bxN_b)|{d}  The distance matrix
            
            if nargin<3
                error('stem_modeltype and correlation_type must be provided');
            end
            if nargin<4
                cross_type=1;
            end
            if nargin<5
                idx_var=[];
            end
            if nargin>=3
                if not(isempty(idx_var))
                    if idx_var<1||idx_var>length(obj.grid)
                        error(['idx_var must be between 1 and ',num2str(length(obj.grid))]);
                    end
                end
            end
            
            if strcmp(correlation_type,'expsphere') && not(strcmp(obj.grid{1}.unit,'deg'))
                error('Coordiantes must be given in lat/lon when the correlation function is expsphere');
            end

            if isempty(idx_var)
                d = obj.get_jointcoordinates();
            else
                d = obj.grid{idx_var}.coordinate;
            end
            if isempty(obj.tap)
                if not(strcmp(correlation_type,'expsphere'))
                    if strcmp(obj.grid{1}.unit,'deg')
                        if strcmp(stem_modeltype,'Emulator')
                            error('The deg unit of measure for grids is not allowed when model_name is ''Emulator''');
                        end
                        DistMat=zeros(size(d,1));
                        for z=1:length(d)
                            DistMat(z,z+1:end)=distance(d(z,:),d(z+1:end,:));
                        end
                        DistMat=DistMat+DistMat';
                    else
                        if not(strcmp(stem_modeltype,'Emulator'))
                            DistMat=squareform(pdist(d));
                        else
                            DistMat=cell(size(d,2),1);
                            for i=1:size(d,2)
                                DistMat{i}=squareform(pdist(d(:,i)));
                            end
                        end
                    end
                else
                    DistMat=cell(2,1);
                    DistMat{1}=zeros(size(d,1));
                    for z=1:length(d)
                        DistMat{1}(z,z+1:end)=distance(d(z,:),d(z+1:end,:));
                    end
                    DistMat{1}=DistMat{1}+DistMat{1}';
                    DistMat{2}=squareform(pdist(d(:,1)));
                end
            else
                if not(strcmp(correlation_type,'expsphere'))
                    if not(strcmp(stem_modeltype,'Emulator'))
                        if (cross_type==1||length(obj.grid)==1||not(isempty(idx_var)))
                            idx_r=[];
                            idx_c=[];
                            elements=[];
                            
                            for z=1:length(d)
                                %evaluate the distance between the z-th coordinate and
                                %the vector ahead (z included)
                                if strcmp(obj.grid{1}.unit,'deg')
                                    dist_vec=distance(d(z,:),d(z:end,:));
                                else
                                    temp=d(z,:);
                                    temp2=d(z:end,:);
                                    temp=repmat(temp,[size(temp2,1),1]);
                                    dist_vec=sqrt(sum((temp-temp2).^2,2));
                                end
                                %IMPORTANT! the distances equal to zero are setted to
                                %eps so they are not confused with the zero generated
                                %by tapering
                                dist_vec(dist_vec==0)=eps;
                                
                                L=dist_vec<=obj.tap;
                                idx_r=cat(1,idx_r,ones(sum(L),1)*z);
                                idx_c=cat(1,idx_c,find(L)+z-1);
                                elements=cat(1,elements,dist_vec(L));
                                %traspose
                                idx_c=cat(1,idx_c,ones(sum(L),1)*z);
                                idx_r=cat(1,idx_r,find(L)+z-1);
                                elements=cat(1,elements,dist_vec(L));
                            end
                            DistMat=sparse(idx_r,idx_c,elements);
                        else
                            Dist=cell(length(obj.grid),1);
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
                                        if strcmp(obj.grid{1}.unit,'deg')
                                            dist_vec=distance(d(z,:),d(z:end,:));
                                        else
                                            temp=d(z,:);
                                            temp2=d(z:end,:);
                                            temp=repmat(temp,[size(temp2,1),1]);
                                            dist_vec=sqrt(sum((temp-temp2).^2,2));
                                        end
                                        
                                        %IMPORTANT! the distances equal to zero are setted to
                                        %eps so they are not confused with the zero generated
                                        %by tapering
                                        dist_vec(dist_vec==0)=eps;
                                        
                                        L=dist_vec<=obj.tap;
                                        idx_r=cat(1,idx_r,ones(sum(L),1)*z);
                                        idx_c=cat(1,idx_c,find(L)+z-1);
                                        elements=cat(1,elements,dist_vec(L));
                                        %traspose
                                        idx_c=cat(1,idx_c,ones(sum(L),1)*z);
                                        idx_r=cat(1,idx_r,find(L)+z-1);
                                        elements=cat(1,elements,dist_vec(L));
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
                    else
                        %emulator case
                        DistMat=cell(size(d,2),1);
                        for i=1:size(d,2)
                            idx_r=[];
                            idx_c=[];
                            elements=[];
                            
                            for z=1:length(d)
                                %evaluate the distance between the z-th coordinate and
                                %the vector ahead (z included)
                                temp=d(z,i);
                                temp2=d(z:end,i);
                                dist_vec=abs(temp-temp2);
                                
                                %IMPORTANT! the distances equal to zero are setted to
                                %eps so they are not confused with the zero generated
                                %by tapering
                                dist_vec(dist_vec==0)=eps;
                                
                                L=dist_vec<=obj.tap;
                                idx_r=cat(1,idx_r,ones(sum(L),1)*z);
                                idx_c=cat(1,idx_c,find(L)+z-1);
                                elements=cat(1,elements,dist_vec(L));
                                %traspose
                                idx_c=cat(1,idx_c,ones(sum(L),1)*z);
                                idx_r=cat(1,idx_r,find(L)+z-1);
                                elements=cat(1,elements,dist_vec(L));
                            end
                            DistMat{i}=sparse(idx_r,idx_c,elements);
                        end
                    end
                else
                    DistMat=cell(2,1);
                    if (cross_type==1||length(obj.grid)==1||not(isempty(idx_var)))
                        idx_r=[];
                        idx_c=[];
                        elements=[];
                        
                        for z=1:length(d)
                            %evaluate the distance between the z-th coordinate and
                            %the vector ahead (z included)
                            if strcmp(obj.grid{1}.unit,'deg')
                                dist_vec=distance(d(z,:),d(z:end,:));
                            else
                                temp=d(z,:);
                                temp2=d(z:end,:);
                                temp=repmat(temp,[size(temp2,1),1]);
                                dist_vec=sqrt(sum((temp-temp2).^2,2));
                            end
                            %IMPORTANT! the distances equal to zero are setted to
                            %eps so they are not confused with the zero generated
                            %by tapering
                            dist_vec(dist_vec==0)=eps;
                            
                            L=dist_vec<=obj.tap;
                            idx_r=cat(1,idx_r,ones(sum(L),1)*z);
                            idx_c=cat(1,idx_c,find(L)+z-1);
                            elements=cat(1,elements,dist_vec(L));
                            %traspose
                            idx_c=cat(1,idx_c,ones(sum(L),1)*z);
                            idx_r=cat(1,idx_r,find(L)+z-1);
                            elements=cat(1,elements,dist_vec(L));
                        end
                        DistMat{1}=sparse(idx_r,idx_c,elements);
                    else
                        Dist=cell(length(obj.grid),1);
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
                                    if strcmp(obj.grid{1}.unit,'deg')
                                        dist_vec=distance(d(z,:),d(z:end,:));
                                    else
                                        temp=d(z,:);
                                        temp2=d(z:end,:);
                                        temp=repmat(temp,[size(temp2,1),1]);
                                        dist_vec=sqrt(sum((temp-temp2).^2,2));
                                    end
                                    
                                    %IMPORTANT! the distances equal to zero are setted to
                                    %eps so they are not confused with the zero generated
                                    %by tapering
                                    dist_vec(dist_vec==0)=eps;
                                    
                                    L=dist_vec<=obj.tap;
                                    idx_r=cat(1,idx_r,ones(sum(L),1)*z);
                                    idx_c=cat(1,idx_c,find(L)+z-1);
                                    elements=cat(1,elements,dist_vec(L));
                                    %traspose
                                    idx_c=cat(1,idx_c,ones(sum(L),1)*z);
                                    idx_r=cat(1,idx_r,find(L)+z-1);
                                    elements=cat(1,elements,dist_vec(L));
                                end
                                Dist{i}=sparse(idx_r,idx_c,elements);
                            else
                                Dist{i}=Dist{1};
                            end
                        end
                        DistMat{1}=blkdiag(Dist{1},Dist{2});
                        for i=3:length(obj.grid)
                            DistMat{1}=blkdiag(DistMat{1},Dist{i});
                        end
                    end
                    
                    idx_r=[];
                    idx_c=[];
                    elements=[];
                    
                    for z=1:length(d)
                        %evaluate the distance between the z-th coordinate and
                        %the vector ahead (z included)
                        temp=d(z,1);
                        temp2=d(z:end,1);
                        dist_vec=abs(temp-temp2);
                        
                    
                        %IMPORTANT! the distances equal to zero are setted to
                        %eps so they are not confused with the zero generated
                        %by tapering
                        dist_vec(dist_vec==0)=eps;
                        
                        L=dist_vec<=obj.tap;
                        idx_r=cat(1,idx_r,ones(sum(L),1)*z);
                        idx_c=cat(1,idx_c,find(L)+z-1);
                        elements=cat(1,elements,dist_vec(L));
                        %traspose
                        idx_c=cat(1,idx_c,ones(sum(L),1)*z);
                        idx_r=cat(1,idx_r,find(L)+z-1);
                        elements=cat(1,elements,dist_vec(L));
                    end
                    DistMat{2}=sparse(idx_r,idx_c,elements);
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
        
        function permute(obj,indices)
            %DESCRIPTION: permute a stem_grid object from the stem_gridlist object
            %
            %INPUT
            %obj        - [stem_gridlist object]    (1x1)
            %indices    - [integer >0]              (1x1) the indices of the stem_grid objects to permute
            %
            %OUTPUT
            %none: the stem_grid object is permuted 
            for i=1:length(obj.grid)
                if iscell(indices)
                    obj.grid{i}.permute(indices{i});
                else
                    obj.grid{i}.permute(indices);
                end
            end
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
            
            for j=1:length(obj.grid{1}.box)
                temp=[];
                for i=1:length(obj.grid)
                    temp=cat(2,temp,obj.grid{i}.box(j));
                end
                if mod(j,2)==1
                    obj.box(j)=min(temp);
                else
                    obj.box(j)=max(temp);
                end
            end
        end
    end
end