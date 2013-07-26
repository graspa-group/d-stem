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


classdef stem_data < handle
    
    %CONSTANTS
    %N_p = n1_p+...+nq_p - total number of point sites
    %N_b = n1_b+...+nq_b - total number of pixel sites
    %N   = N_p+N_b - total number of observation sites
    %N_b = n1_b+...+nq_b+n1_b+...+nq_b - total number of covariates
    %S   = 2 if both point and pixel data are considered. S = 1 if only point data are considered.
    %T   - number of temporal steps
    %TT  = T if the space-time varying coefficients are time-variant and TT=1 if they are time-invariant
    %p   - dimension of the latent temporal variable z
    
    properties
        stem_varset_p=[];       %[stem_varset object]   (1x1) stem_varset object for the point variables
        stem_varset_b=[];       %[stem_varset object]   (1x1) stem_varset object for the pixel variables
        stem_gridlist_p=[];     %[stem_gridlist object] (1x1) stem_gridlist object for the point variables        
        stem_gridlist_b=[];     %[stem_gridlist object] (1x1) stem_gridlist object for the pixel variables
        stem_datestamp=[]       %[stem_datestamp object](1x1) stem_datestamp object with information on time steps
        stem_crossval=[];       %[stem_crossval object] (1x1) stem_crossval object with information on crossvalidation
        
        shape=[];               %[struct]               (1x1) boundary of the geographic region loaded from a shape file
        simulated=0;            %[boolean]              (1x1) 1: the data have been simulated; 0: observed data
        pixel_correlated=0;     %[boolean]              (1x1) 1: the pixel data are cross-correlated
        X_z=[];                 %[double]               (NxpxTT) the full X_z matrix
    end
    
    properties (SetAccess = private) 
        Y=[];                   %[double]     (NxT) the full observation matrix
        X_bp=[];                %[double]     (Nx1xTT) the full X_bp matrix
        X_beta=[];              %[double]     (NxN_bxTT) the full X_beta matrix
        X_p=[];                 %[double]     (Nx1xTTxK) the full X_p matrix
        DistMat_p=[];           %[double]     (N_pxN_p) distance matrix of the point sites
        DistMat_b=[];           %[double]     (N_bxN_b) distance matrix of the pixel sites
        M=[];                   %[integer >1] (N_px1) vector of indices of the pixel mapped on the point sites
        %flags
        can_reset=0;            %[boolean]    (1x1) 1: data are saved on disk and can be reloaded; 0: data are only on RAM
        X_bp_tv=0;              %[boolean]    (1x1) 1: X_bp is time variant; 0: otherwise
        X_beta_tv=0;            %[boolean]    (1x1) 1: X_beta is time variant; 0: otherwise
        X_z_tv=0;               %[boolean]    (1x1) 1: X_z is time variant; 0: otherwise
        X_p_tv=0;               %[boolean]    (1x1) 1: X_p is time variant; 0: otherwise
        X_tv=0;                 %[boolean]    (1x1) 1: at least one between X_bp, X_beta, X_z and X_p is time variant; 0:otherwise
    end
    
    methods
        
        function obj = stem_data(stem_varset_p,stem_gridlist_p,stem_varset_b,stem_gridlist_b,stem_datestamp,shape,can_reset,stem_crossval,pixel_correlated)
            %DESCRIPTION: is the constructor of the class stem_data
            %
            %INPUT
            %
            %stem_varset_p      - [stem_varset object]    (1x1) stem_varset object for the point variables
            %stem_gridlist_p    - [stem_gridlist object]  (1x1) stem_gridlist object for the point variables  
            %<stem_varset_b>    - [stem_varset object]    (1x1) (default: []) stem_varset object for the pixel variables
            %<stem_gridlist_b>  - [stem_gridlist object]  (1x1) (default: []) stem_gridlist object for the pixel variables
            %<stem_datestamp>   - [stem_datestamp object] (1x1) (default: []) stem_datestamp object with information on time steps
            %<shape>            - [struct]                      (default: world boundaries) geographic data structure loaded from a shapefile with the boundary of the geographic region
            %<can_reset>        - [boolean]               (1x1) (default: 0) 1: the data are saved on disk and they can be reloaded using the method reset of this class after, for example, data transformation
            %<stem_crossval>    - [stem_crossval object]  (1x1) (default: []) stem_crossval object with information on crossvalidation
            %
            %OUTPUT
            %obj                - [stem_data object]      (1x1) the stem_data object
            
            if nargin<2
                error('Not enough input parameters');
            end
            obj.stem_varset_p=stem_varset_p;
            obj.stem_gridlist_p=stem_gridlist_p;
            if nargin==3
                error('stem_gridlist_b must be provided');
            end
            if nargin>2
                if not(isempty(stem_varset_b))
                    obj.stem_varset_b=stem_varset_b;
                    obj.stem_gridlist_b=stem_gridlist_b;
                end
            end
            if nargin>4
                obj.stem_datestamp=stem_datestamp;
            end

            if nargin>=6
                if not(isempty(shape))
                    obj.shape=shape;
                    obj.shape(1,1).Geometry='Line';
                else
                    obj.shape=shaperead('landareas.shp');
                end
            else
                obj.shape=shaperead('landareas.shp');
            end
            
            if nargin>=7
                if not(isempty(can_reset))
                    obj.can_reset=can_reset;
                else
                    obj.can_reset=0;
                end
            else
                obj.can_reset=0;
            end
            
            if nargin>=8
                if not(isempty(stem_crossval))
                    obj.stem_crossval=stem_crossval;
                    idx_var = obj.stem_varset_p.get_Y_index(obj.stem_crossval.variable_name);
                    if isempty(idx_var)
                        error('Cross-validation variable not found');
                    end
                end
            end
            
            if nargin>=9
                if not(isempty(pixel_correlated))
                    obj.pixel_correlated=pixel_correlated;
                end
            end
            
%             if obj.can_reset
%                 obj.set_original;
%             end
            
            obj.update_data();
            obj.update_distance();
            if not(isempty(stem_varset_b))
                obj.update_M();
            end
            obj.remove_duplicated_sites();
        end
        
        function update_data(obj)
            %DESCRIPTION: generates the matrices Y, X_bp, X_beta, X_z and X_p 
            %
            %INPUT
            %obj - [stem_data object] (1x1) the stem_data object
            %
            %OUTPUT
            %
            %none: the matrices listed above are updated
            
            disp('Generating data matrices...');
            %Y
            Y=[];
            for i=1:length(obj.stem_varset_p.Y)
                Y=cat(1,Y,obj.stem_varset_p.Y{i});
            end
            if not(isempty(obj.stem_varset_b))
                for i=1:length(obj.stem_varset_b.Y)
                    Y=cat(1,Y,obj.stem_varset_b.Y{i});
                end
            end
            obj.Y=Y;
            clear Y;
            %X_bp
            X_bp=[];
            if not(isempty(obj.stem_varset_b))
                if not(isempty(obj.stem_varset_b.X_bp))&&not(isempty(obj.stem_varset_p.X_bp))
                    done=0;
                    if size(obj.stem_varset_p.X_bp{1},3)==1 && size(obj.stem_varset_b.X_bp{1},3)==obj.T
                        for i=1:length(obj.stem_varset_p.X_bp)
                            X_bp=cat(1,X_bp,repmat(obj.stem_varset_p.X_bp{i},[1,1,obj.T]));
                        end
                        for i=1:length(obj.stem_varset_b.X_bp)
                            X_bp=cat(1,X_bp,obj.stem_varset_b.X_bp{i});
                        end
                        done=1;
                    end
                    if size(obj.stem_varset_p.X_bp{1},3)==obj.T && size(obj.stem_varset_b.X_bp{1},3)==1 && not(done)
                        for i=1:length(obj.stem_varset_p.X_bp)
                            X_bp=cat(1,X_bp,obj.stem_varset_p.X_bp{i});
                        end
                        for i=1:length(obj.stem_varset_b.X_bp)
                            X_bp=cat(1,X_bp,repmat(obj.stem_varset_b.X_bp{i},[1,1,obj.T]));
                        end
                        done=1;
                    end
                    if size(obj.stem_varset_p.X_bp{1},3)==size(obj.stem_varset_b.X_bp{1},3) && not(done)
                        for i=1:length(obj.stem_varset_p.X_bp)
                            X_bp=cat(1,X_bp,obj.stem_varset_p.X_bp{i});
                        end
                        for i=1:length(obj.stem_varset_b.X_bp)
                            X_bp=cat(1,X_bp,obj.stem_varset_b.X_bp{i});
                        end
                        done=1;
                    end
                end
                obj.X_bp=X_bp;
                if size(obj.X_bp,3)>1
                    obj.X_bp_tv=1;
                end
                clear X_bp;
            end
%             else
%                 for i=1:length(obj.stem_varset_p.X_bp)
%                     X_bp=cat(1,X_bp,obj.stem_varset_p.X_bp{i});
%                 end
%             end

            %X_beta
            X_beta=[];
            if not(isempty(obj.stem_varset_b))
                if not(isempty(obj.stem_varset_b.X_beta))&&not(isempty(obj.stem_varset_p.X_beta))
                    done=0;
                    if size(obj.stem_varset_p.X_beta{1},3)==1 && size(obj.stem_varset_b.X_beta{1},3)==obj.T
                        for t=1:obj.T
                            X_temp=[];
                            for i=1:length(obj.stem_varset_p.X_beta)
                                X_temp=blkdiag(X_temp,obj.stem_varset_p.X_beta{i});
                            end
                            for i=1:length(obj.stem_varset_b.X_beta)
                                X_temp=blkdiag(X_temp,obj.stem_varset_b.X_beta{i}(:,:,t));
                            end
                            X_beta(:,:,t)=X_temp;
                        end
                        done=1;
                    end
                    if size(obj.stem_varset_p.X_beta{1},3)==obj.T && size(obj.stem_varset_b.X_beta{1},3)==1 && not(done)
                        for t=1:obj.T
                            X_temp=[];
                            for i=1:length(obj.stem_varset_p.X_beta)
                                X_temp=blkdiag(X_temp,obj.stem_varset_p.X_beta{i}(:,:,t));
                            end
                            for i=1:length(obj.stem_varset_b.X_beta)
                                X_temp=blkdiag(X_temp,obj.stem_varset_b.X_beta{i});
                            end
                            X_beta(:,:,t)=X_temp;
                        end
                        done=1;
                    end
                    if (size(obj.stem_varset_p.X_beta{1},3)==obj.T)&&(size(obj.stem_varset_b.X_beta{1},3)==obj.T) && not(done)
                        for t=1:size(obj.stem_varset_p.X_beta{1},3)
                            X_temp=[];
                            for i=1:length(obj.stem_varset_p.X_beta)
                                X_temp=blkdiag(X_temp,obj.stem_varset_p.X_beta{i}(:,:,t));
                            end
                            for i=1:length(obj.stem_varset_b.X_beta)
                                X_temp=blkdiag(X_temp,obj.stem_varset_b.X_beta{i}(:,:,t));
                            end
                            X_beta(:,:,t)=X_temp;
                        end
                        done=1;
                    end
                    if (size(obj.stem_varset_p.X_beta{1},3)==1)&&(size(obj.stem_varset_b.X_beta{1},3)==1) && not(done)
                        X_temp=[];
                        for i=1:length(obj.stem_varset_p.X_beta)
                            X_temp=blkdiag(X_temp,obj.stem_varset_p.X_beta{i});
                        end
                        for i=1:length(obj.stem_varset_b.X_beta)
                            X_temp=blkdiag(X_temp,obj.stem_varset_b.X_beta{i});
                        end
                        X_beta=X_temp;
                        done=1;
                    end
                else
                    if not(isempty(obj.stem_varset_p.X_beta))
                        done=0;
                        if size(obj.stem_varset_p.X_beta{1},3)==obj.T
                            for t=1:size(obj.stem_varset_p.X_beta{1},3)
                                X_temp=[];
                                for i=1:length(obj.stem_varset_p.X_beta)
                                    X_temp=blkdiag(X_temp,obj.stem_varset_p.X_beta{i}(:,:,t));
                                end
                                %X_temp=cat(1,X_temp,zeros(obj.stem_varset_b.N,size(X_temp,2),size(X_temp,3)));
                                X_beta(:,:,t)=X_temp;
                            end
                            done=1;
                        end
                        if size(obj.stem_varset_p.X_beta{1},3)==1 && not(done)
                            X_temp=[];
                            for i=1:length(obj.stem_varset_p.X_beta)
                                X_temp=blkdiag(X_temp,obj.stem_varset_p.X_beta{i});
                            end
                            %X_temp=cat(1,X_temp,zeros(obj.stem_varset_b.N,size(X_temp,2),size(X_temp,3)));
                            X_beta=X_temp;
                            done=1;
                        end
                    end
                end
            else
                if not(isempty(obj.stem_varset_p.X_beta))
                    done=0;
                    if size(obj.stem_varset_p.X_beta{1},3)==obj.T
                        nbeta=0;
                        for i=1:length(obj.stem_varset_p.X_beta)
                            nbeta=nbeta+size(obj.stem_varset_p.X_beta{i},2);
                        end
                        X_beta=zeros(obj.N,nbeta,obj.T);
                        for t=1:size(obj.stem_varset_p.X_beta{1},3)
                            X_temp=[];
                            for i=1:length(obj.stem_varset_p.X_beta)
                                X_temp=blkdiag(X_temp,obj.stem_varset_p.X_beta{i}(:,:,t));
                            end
                            X_beta(:,:,t)=X_temp;
                        end
                        done=1;
                    end
                    if size(obj.stem_varset_p.X_beta{1},3)==1 && not(done)
                        X_temp=[];
                        for i=1:length(obj.stem_varset_p.X_beta)
                            X_temp=blkdiag(X_temp,obj.stem_varset_p.X_beta{i});
                        end
                        X_beta=X_temp;
                        done=1;
                    end
                end
            end
            obj.X_beta=X_beta;
            if size(obj.X_beta,3)>1
                obj.X_beta_tv=1;
            end
            clear X_beta;

            %X_z
            X_z=[];
            if not(isempty(obj.stem_varset_b))
                if not(isempty(obj.stem_varset_b.X_z))&&not(isempty(obj.stem_varset_p.X_z))
                    done=0;
                    if size(obj.stem_varset_p.X_z{1},3)==1 && size(obj.stem_varset_b.X_z{1},3)==obj.T
                        for t=1:obj.T
                            X_temp=[];
                            for i=1:length(obj.stem_varset_p.X_z)
                                X_temp=blkdiag(X_temp,obj.stem_varset_p.X_z{i});
                            end
                            for i=1:length(obj.stem_varset_b.X_z)
                                X_temp=blkdiag(X_temp,obj.stem_varset_b.X_z{i}(:,:,t));
                            end
                            X_z(:,:,t)=X_temp;
                        end
                        done=1;
                    end
                    if size(obj.stem_varset_p.X_z{1},3)==obj.T && size(obj.stem_varset_b.X_z{1},3)==1 && not(done)
                        for t=1:obj.T
                            X_temp=[];
                            for i=1:length(obj.stem_varset_p.X_z)
                                X_temp=blkdiag(X_temp,obj.stem_varset_p.X_z{i}(:,:,t));
                            end
                            for i=1:length(obj.stem_varset_b.X_z)
                                X_temp=blkdiag(X_temp,obj.stem_varset_b.X_z{i});
                            end
                            X_z(:,:,t)=X_temp;
                        end
                        done=1;
                    end
                    if (size(obj.stem_varset_p.X_z{1},3)==obj.T)&&(size(obj.stem_varset_b.X_z{1},3)==obj.T) && not(done)
                        for t=1:size(obj.stem_varset_p.X_z{1},3)
                            X_temp=[];
                            for i=1:length(obj.stem_varset_p.X_z)
                                X_temp=blkdiag(X_temp,obj.stem_varset_p.X_z{i}(:,:,t));
                            end
                            for i=1:length(obj.stem_varset_b.X_z)
                                X_temp=blkdiag(X_temp,obj.stem_varset_b.X_z{i}(:,:,t));
                            end
                            X_z(:,:,t)=X_temp;
                        end
                        done=1;
                    end
                    if (size(obj.stem_varset_p.X_z{1},3)==1)&&(size(obj.stem_varset_b.X_z{1},3)==1) && not(done)
                        X_temp=[];
                        for i=1:length(obj.stem_varset_p.X_z)
                            X_temp=blkdiag(X_temp,obj.stem_varset_p.X_z{i});
                        end
                        for i=1:length(obj.stem_varset_b.X_z)
                            X_temp=blkdiag(X_temp,obj.stem_varset_b.X_z{i});
                        end
                        X_z=X_temp;
                        done=1;
                    end
                else
                    if not(isempty(obj.stem_varset_p.X_z))
                        done=0;
                        if size(obj.stem_varset_p.X_z{1},3)==obj.T
                            for t=1:size(obj.stem_varset_p.X_z{1},3)
                                X_temp=[];
                                for i=1:length(obj.stem_varset_p.X_z)
                                    X_temp=blkdiag(X_temp,obj.stem_varset_p.X_z{i}(:,:,t));
                                end
                                X_z(:,:,t)=X_temp;
                            end
                            done=1;
                        end
                        if size(obj.stem_varset_p.X_z{1},3)==1 && not(done)
                            X_temp=[];
                            for i=1:length(obj.stem_varset_p.X_z)
                                X_temp=blkdiag(X_temp,obj.stem_varset_p.X_z{i});
                            end
                            X_z=X_temp;
                            done=1;
                        end
                        %X_z=cat(1,X_z,zeros(size(obj.Y,1)-size(X_z,1),size(X_z,2),size(X_z,3)));
                    end
                end
            else
                if not(isempty(obj.stem_varset_p.X_z))
                    done=0;
                    if size(obj.stem_varset_p.X_z{1},3)==obj.T
                        for t=1:size(obj.stem_varset_p.X_z{1},3)
                            X_temp=[];
                            for i=1:length(obj.stem_varset_p.X_z)
                                X_temp=blkdiag(X_temp,obj.stem_varset_p.X_z{i}(:,:,t));
                            end
                            X_z(:,:,t)=X_temp;
                        end
                        done=1;
                    end
                    if size(obj.stem_varset_p.X_z{1},3)==1 && not(done)
                        X_temp=[];
                        for i=1:length(obj.stem_varset_p.X_z)
                            X_temp=blkdiag(X_temp,obj.stem_varset_p.X_z{i});
                        end
                        X_z=X_temp;
                        done=1;
                    end
                end
            end
            if not(isempty(X_z))
                obj.X_z=X_z;
                if size(obj.X_z,3)>1
                    obj.X_z_tv=1;
                end
            end
            clear X_z

            %X_p
            if not(isempty(obj.stem_varset_p.X_p))
                X_p=[];
                for i=1:length(obj.stem_varset_p.X_p)
                    X_p=cat(1,X_p,obj.stem_varset_p.X_p{i});
                end
                obj.X_p=X_p;
                if size(obj.X_p,3)>1
                    obj.X_p_tv=1;
                end
                clear X_p;
            end
            obj.X_tv=obj.X_bp_tv | obj.X_z_tv | obj.X_p_tv;
            disp('Generation ended.');
        end
        
        function update_M(obj)
            %DESCRIPTION: generates the vector M
            %
            %INPUT
            %obj - [stem_data object] (1x1) the stem_data object
            %
            %OUTPUT         
            %
            %none: the vector M is generated ad updated
            
            disp('Generating M replication vector...');
            M=[];
            blocks=[0 cumsum(obj.stem_varset_b.dim)];
            for j=1:obj.stem_varset_p.nvar
                dmax=distdim(distance(0,0,obj.stem_gridlist_b.grid{j}.pixel_side_w,obj.stem_gridlist_b.grid{j}.pixel_side_w), obj.stem_gridlist_p.grid{1}.unit, 'km');
                for i=1:size(obj.stem_gridlist_p.grid{j}.coordinate,1)
                    d=distdim(distance(obj.stem_gridlist_p.grid{j}.coordinate(i,:),obj.stem_gridlist_b.grid{j}.coordinate), obj.stem_gridlist_p.grid{1}.unit, 'km');
                    [m,idx]=min(d);
                    if d>dmax
                        %warning(['Point ',num2str(i),' of point variable ',num2str(j),' does not belong to any pixel. The nearest pixel at ',num2str(m),' km is considered']);
                    end
                    M=[M;idx+blocks(j)];
                end
            end
            obj.M=M;
            disp('Generation ended.');
        end
        
        function update_distance(obj,type,force)
            %DESCRIPTION: generates the distance matrices
            %
            %INPUT
            %obj    - [stem_data object] (1x1) the stem_data object
            %<type> - [string]           (1x1) (Default: 'both') 'point': only the distance matrix for the point data is evaluated. 
            %                                                    'pixel': only the distance matrix for the pixel data is evaluated.
            %                                                    'both':  both the matrices are evaluated.
            %
            %OUTPUT         
            %
            %none: the DistMat_p and DistMat_b property are generated and updated

            if nargin<2
                type='both';
            end
            if nargin<3
                force=0;
            end
            cmp=strcmp(type,{'both','point','pixel'});
            if sum(cmp)==0
                error('type must be point, pixel or both');
            end
            
            if strcmp(type,'point')||strcmp(type,'both')
                if not(isempty(obj.stem_varset_p.X_p))||force
                    disp('Generating point distance matrices...');
                    obj.DistMat_p=obj.stem_gridlist_p.get_distance_matrix();
                    disp('Generation ended.');
                end
            end
            if strcmp(type,'pixel')||strcmp(type,'both')
                if not(isempty(obj.stem_gridlist_b))&&not(isempty(obj.stem_varset_b.X_bp))
                    disp('Generating pixel data distance matrices...');
                    obj.DistMat_b=obj.stem_gridlist_b.get_distance_matrix(obj.pixel_correlated);
                    disp('Generation ended.');
                end
            end
        end     
       
        %Data transform
        function detrend_Y(obj)
            %DESCRIPTION: remove the mean from each time series in Y
            %
            %INPUT
            %obj - [stem_data object] (1x1) the stem_data object
            %
            %OUTPUT
            %
            %none: the Y property is updated
            
            disp('Point level data detrend started...');
            obj.stem_varset_p.detrend;
            disp('Point level data detrend ended.');
            if not(isempty(obj.stem_varset_b))
                disp('Pixel data detrend started...');
                obj.stem_varset_b.detrend;
                disp('Pixel data detrend ended.');
            end
            disp('Updtaing data matrices after detrend...');
            obj.update_data;            
        end
        
        function standardize_Y(obj)
            %DESCRIPTION: each time series in Y is standardized
            %
            %INPUT
            %obj - [stem_data object] (1x1) the stem_data object
            %
            %OUTPUT
            %
            %none: the Y property is updated            
            
            disp('Point level data site by site standardization started...');
            obj.stem_varset_p.standardize_Y;
            disp('Point level data site by site standardization ended.');
            if not(isempty(obj.stem_varset_b))
                disp('Pixel data site by site standardization started...');
                obj.stem_varset_b.standardize_Y;
                disp('Pixel data site by site standardization ended.');
            end
            disp('Updtaing data matrices after site by site standardization...');
            obj.update_data;
        end
        
        function standardize(obj)
            %DESCRIPTION: standardize the matrices Y, X_bp, X_beta, X_z and X_p with respect to their overall mean and overall standard deviation
            %
            %INPUT
            %obj - [stem_data object] (1x1) the stem_data object
            %
            %OUTPUT
            %
            %none: the matrices listed above are updated
            
            disp('Point level data standardization started...');
            obj.stem_varset_p.standardize;
            disp('Point level data standardization ended.');
            if not(isempty(obj.stem_varset_b))
                disp('Pixel data standardization started...');
                obj.stem_varset_b.standardize;
                disp('Pixel data standardization ended.');
            end
            disp('Updtaing data matrices after standardization...');
            obj.update_data;
        end
        
        function log_transform(obj)
            %DESCRIPTION: log-transforms the matrix Y
            %
            %INPUT
            %obj - [stem_data object] (1x1) the stem_data object
            %
            %OUTPUT
            %
            %none: the matrix Y is updated            
            
            disp('Point level data log-transformation started...');
            obj.stem_varset_p.log_transform;
            disp('Point level data log-transformation ended.');
            if not(isempty(obj.stem_varset_b))
                disp('Pixel data log-transformation started...');
                obj.stem_varset_b.log_transform;
                disp('Pixel data log-transformation ended.');
            end
            disp('Updtaing data matrices after log-transformation...');
            obj.update_data;
        end
        
        function time_average(obj,n_steps)
            %DESCRIPTION: computes time averages of n_steps for the matrice the matrices Y, X_bp, X_beta, X_z and X_p
            %
            %INPUT
            %obj        - [stem_data object] (1x1) the stem_data object
            %n_steps    - [integer >0]       (1x1) the number of temporal steps to average
            %
            %OUTPUT
            %
            %none: the matrices listed above are updated
            
            if nargin<2
                error('The subsampling factor must be provided');
            end
            if n_steps<=1
                error('n_steps must be greater than 1');
            end
            if round(n_steps)~=n_steps
                error('n_steps must be an integer value');
            end
   
            indices=0:n_steps:obj.T;
            if indices(end)~=obj.T
                indices=[indices,obj.T];
            end
            obj.stem_datestamp.average_stamps(indices);
            
            disp('Time averaging started...');
            for i=1:length(obj.stem_varset_p.Y)
                for j=1:length(indices)-1
                    Y_temp{i}(:,j)=nanmean(obj.stem_varset_p.Y{i}(:,indices(j)+1:indices(j+1)),2);
                end
            end
            obj.stem_varset_p.Y=Y_temp;
            clear Y_temp

            if not(isempty(obj.stem_varset_p.X_bp))
                if obj.stem_varset_p.X_bp_tv
                    for i=1:length(obj.stem_varset_p.X_bp)
                        for j=1:length(indices)-1
                            X_bp_temp{i}(:,:,j)=nanmean(obj.stem_varset_p.X_bp{i}(:,:,indices(j)+1:indices(j+1)),3);
                        end
                    end
                    obj.stem_varset_p.X_bp=X_bp_temp;
                    clear X_bp_temp
                end
            end
            
            if not(isempty(obj.stem_varset_p.X_beta))
                if obj.stem_varset_p.X_beta_tv
                    for i=1:length(obj.stem_varset_p.X_beta)
                        for j=1:length(indices)-1
                            X_beta_temp{i}(:,:,j)=nanmean(obj.stem_varset_p.X_beta{i}(:,:,indices(j)+1:indices(j+1)),3);
                        end
                    end
                    obj.stem_varset_p.X_beta=X_beta_temp;
                    clear X_beta_temp
                end
            end
            
            if not(isempty(obj.stem_varset_p.X_z))
                if obj.stem_varset_p.X_z_tv
                    for i=1:length(obj.stem_varset_p.X_z)
                        for j=1:length(indices)-1
                            X_z_temp{i}(:,:,j)=nanmean(obj.stem_varset_p.X_z{i}(:,:,indices(j)+1:indices(j+1)),3);
                        end
                    end
                    obj.stem_varset_p.X_z=X_z_temp;
                    clear X_z_temp
                end
            end   
            
            if not(isempty(obj.stem_varset_p.X_p))
                if obj.stem_varset_p.X_p_tv
                    for i=1:length(obj.stem_varset_p.X_p)
                        for j=1:length(indices)-1
                            X_p_temp{i}(:,:,j,:)=nanmean(obj.stem_varset_p.X_p{i}(:,:,indices(j)+1:indices(j+1),:),3);
                        end
                    end
                    obj.stem_varset_p.X_p=X_p_temp;        
                    clear X_p_temp
                end
            end
            
            
            if not(isempty(obj.stem_varset_b))
                for i=1:length(obj.stem_varset_b.Y)
                    for j=1:length(indices)-1
                        Y_temp{i}(:,j)=nanmean(obj.stem_varset_b.Y{i}(:,indices(j)+1:indices(j+1)),2);
                    end
                end
                obj.stem_varset_b.Y=Y_temp;
                clear Y_temp
                
                if not(isempty(obj.stem_varset_p.X_bp))
                    if obj.stem_varset_b.X_bp_tv
                        for i=1:length(obj.stem_varset_b.X_bp)
                            for j=1:length(indices)-1
                                X_bp_temp{i}(:,:,j)=nanmean(obj.stem_varset_b.X_bp{i}(:,:,indices(j)+1:indices(j+1)),3);
                            end
                        end
                        obj.stem_varset_b.X_bp=X_bp_temp;
                        clear X_bp_temp
                    end
                end
                
                if not(isempty(obj.stem_varset_b.X_beta))
                    if obj.stem_varset_b.X_beta_tv
                        for i=1:length(obj.stem_varset_b.X_beta)
                            for j=1:length(indices)-1
                                X_beta_temp{i}(:,:,j)=nanmean(obj.stem_varset_b.X_beta{i}(:,:,indices(j)+1:indices(j+1)),3);
                            end
                        end
                        obj.stem_varset_b.X_beta=X_beta_temp;
                        clear X_beta_temp
                    end
                end
                
                if not(isempty(obj.stem_varset_b.X_z))
                    if obj.stem_varset_b.X_z_tv
                        for i=1:length(obj.stem_varset_b.X_z)
                            for j=1:length(indices)-1
                                X_z_temp{i}(:,:,j)=nanmean(obj.stem_varset_b.X_z{i}(:,:,indices(j)+1:indices(j+1)),3);
                            end
                        end
                        obj.stem_varset_b.X_z=X_z_temp;
                    end
                end
                
            end

            obj.update_data;    
            disp('Time averaging ended.');
        end
        
        function time_crop(obj,dates_or_indices)
            %DESCRIPTION: crop the matrices Y, X_bp, X_beta, X_z and X_p with respect to time
            %
            %INPUT
            %obj                 - [stem_data object]       (1x1) the stem_data object
            %dates_or_indices    - [string | integer >0]    {2x1}|(dTx1) a {2x1} cell vector of starting and ending date in the format dd-mm-yyyy or a dTx1 vector of (possibly non consecutive) temporal indices
            %
            %OUTPUT
            %
            %none: the matrices listed above are updated            

            if iscell(dates_or_indices)
                start_date=datenum(dates_or_indices{1},'dd-mm-yyyy');
                end_date=datenum(dates_or_indices{2},'dd-mm-yyyy');
                if start_date>end_date
                    error('Starting date must be lower or equal to ending date');
                end
                sd=find(obj.stem_datestamp.stamp==start_date);
                ed=find(obj.stem_datestamp.stamp==end_date);
                if isempty(sd)||isempty(ed)
                    error('Dates out of range');
                end
                indices=sd:ed;
            else
                if min(dates_or_indices)<0
                    error('Starting date must be higher than zero');
                end
                if max(dates_or_indices)>obj.T
                    error('Ending date must be lower or equal to the total number of temporal steps');
                end
                indices=dates_or_indices;
            end
            
            disp('Time crop started...');
            Y=obj.stem_varset_p.Y;
            for i=1:length(Y)
                Y{i}=Y{i}(:,indices);
            end
            obj.stem_varset_p.Y=Y;
            if not(isempty(obj.stem_varset_p.X_beta))
                if obj.stem_varset_p.X_beta_tv
                    X_beta=obj.stem_varset_p.X_beta;
                    for i=1:length(X_beta)
                        X_beta{i}=X_beta{i}(:,:,indices);
                    end
                    obj.stem_varset_p.X_beta=X_beta;
                end
            end
            if not(isempty(obj.stem_varset_p.X_bp))
                if obj.stem_varset_p.X_bp_tv
                    X_bp=obj.stem_varset_p.X_bp;
                    for i=1:length(X_bp)
                        X_bp{i}=X_bp{i}(:,:,indices);
                    end
                    obj.stem_varset_p.X_bp=X_bp;
                end
            end               
            if not(isempty(obj.stem_varset_p.X_z))
                if obj.stem_varset_p.X_z_tv
                    X_z=obj.stem_varset_p.X_z;
                    for i=1:length(X_z)
                        X_z{i}=X_z{i}(:,:,indices);
                    end
                    obj.stem_varset_p.X_z=X_z;
                end
            end   
            if not(isempty(obj.stem_varset_p.X_p))
                if obj.stem_varset_p.X_p_tv
                    X_p=obj.stem_varset_p.X_p;
                    for i=1:length(X_p)
                        X_p{i}=X_p{i}(:,:,indices,:);
                    end
                    obj.stem_varset_p.X_p=X_p;
                end
            end
            
            if not(isempty(obj.stem_varset_b))
                if not(isempty(obj.stem_varset_b.Y))
                    Y=obj.stem_varset_b.Y;
                    for i=1:length(Y)
                        Y{i}=Y{i}(:,indices);
                    end
                    obj.stem_varset_b.Y=Y;
                    if not(isempty(obj.stem_varset_b.X_beta))
                        if obj.stem_varset_b.X_beta_tv
                            X_beta=obj.stem_varset_b.X_beta;
                            for i=1:length(obj.stem_varset_b.X_beta)
                                X_beta{i}=X_beta{i}(:,:,indices);
                            end
                            obj.stem_varset_b.X_beta=X_beta;
                        end
                    end
                    if not(isempty(obj.stem_varset_b.X_bp))
                        if obj.stem_varset_b.X_bp_tv
                            X_bp=obj.stem_varset_b.X_bp;
                            for i=1:length(X_bp)
                                X_bp{i}=X_bp{i}(:,:,indices);
                            end
                            obj.stem_varset_b.X_z=X_bp;
                        end
                    end                    
                    if not(isempty(obj.stem_varset_b.X_z))
                        if obj.stem_varset_b.X_z_tv
                            X_z=obj.stem_varset_b.X_z;
                            for i=1:length(X_z)
                                X_z{i}=X_z{i}(:,:,indices);
                            end
                            obj.stem_varset_b.X_z=X_z;
                        end
                    end
                end
            end
            
            obj.stem_datestamp.subset_stamps(indices);
            %looks for line of all missing for the sparse grids of the point data
            changed=0;
            for i=1:obj.stem_varset_p.nvar
                indices=sum(isnan(obj.stem_varset_p.Y{i}),2)==size(obj.stem_varset_p.Y{i},2);
                if sum(indices)>0
                    obj.stem_varset_p.Y{i}(indices,:)=[];
                    if not(isempty(obj.stem_varset_p.X_beta))
                        obj.stem_varset_p.X_beta{i}(indices,:,:)=[];
                    end
                    if not(isempty(obj.stem_varset_p.X_z))
                        obj.stem_varset_p.X_z{i}(indices,:,:)=[];
                    end
                    if not(isempty(obj.stem_varset_p.X_p))
                        obj.stem_varset_p.X_p{i}(indices,:,:,:)=[];
                    end
                    if not(isempty(obj.stem_varset_p.X_bp))
                        obj.stem_varset_p.X_bp{i}(indices,:,:)=[];
                    end
                    obj.stem_gridlist_p.grid{i}.coordinate(indices,:)=[];
                    disp(['Deleted ',num2str(sum(indices)),' site(s) for the point variable ',obj.stem_varset_p.Y_name{i},' due to all missing.']);
                    changed=1;
                end
            end
            if changed
                disp('Updating point distance matrix after time crop...');
                obj.update_distance('point'); %only point because the pixel data are not deleted from the data matrix even if they are NaN for all th time steps
                disp('Update ended.');
                if not(isempty(obj.stem_varset_b))
                    disp('Updating M replication vector after time crop...');
                    obj.update_M;
                    disp('Update ended.');
                end
            end
            disp('Updating data matrix after time crop...');
            obj.update_data;
            disp('Update ended.');
            disp('Time crop ended.');
        end     
        
        function space_crop(obj,box)
            %DESCRIPTION: crop the matrices Y, X_bp, X_beta, X_z and X_p with respect to space
            %
            %INPUT
            %obj                 - [stem_data object]   (1x1) the stem_data object
            %box                 - [double]             (4x1) the geographic box [lat_min,lat_max,lon_min,lon_max]
            %
            %OUTPUT
            %
            %none: the matrices listed above are updated   
            
            lat_min=box(1);
            lat_max=box(2);
            lon_min=box(3);
            lon_max=box(4);
            if nargin<2
                error('The bounding-box must be provided');
            end
            if lat_min>lat_max
                error('The lat_min value must be lower than the lat_max value');
            end
            if lon_min>lon_max
                error('The lon_min value must be lower than the lon_max value');
            end
            disp('Point level data space crop started...');
            for i=1:obj.stem_varset_p.nvar
                GY=obj.stem_gridlist_p.grid{i}.coordinate;
                indices = GY(:,2) >= lon_min & ...
                    GY(:,2) <= lon_max & ...
                    GY(:,1) >= lat_min & ...
                    GY(:,1) <= lat_max;
                if sum(indices)>0
                    if strcmp(obj.stem_gridlist_p.grid{i}.grid_type,'regular')
                        grid_lat=obj.stem_gridlist_p.grid{i}.coordinate(:,1);
                        grid_lon=obj.stem_gridlist_p.grid{i}.coordinate(:,2);
                        grid_lat=reshape(grid_lat,obj.stem_gridlist_p.grid{i}.grid_size);
                        grid_lat=grid_lat(:,1);
                        grid_lon=reshape(grid_lon,obj.stem_gridlist_p.grid{i}.grid_size);
                        grid_lon=grid_lon(1,:);
                        grid_lat(grid_lat<lat_min)=[];
                        grid_lat(grid_lat>lat_max)=[];
                        grid_lon(grid_lon<lon_min)=[];
                        grid_lon(grid_lon>lon_max)=[];
                    end
                    obj.stem_gridlist_p.grid{i}.coordinate=obj.stem_gridlist_p.grid{i}.coordinate(indices,:);
                    
                    if strcmp(obj.stem_gridlist_p.grid{i}.grid_type,'regular')
                        obj.stem_gridlist_p.grid{i}.grid_size=[length(grid_lat),length(grid_lon)];
                    end                    

                    obj.stem_varset_p.Y{i}=obj.stem_varset_p.Y{i}(indices,:);
                    
                    if not(isempty(obj.stem_varset_p.X_bp))
                        obj.stem_varset_p.X_bp{i}=obj.stem_varset_p.X_bp{i}(indices,:,:);
                    end
                    
                    if not(isempty(obj.stem_varset_p.X_beta))
                        obj.stem_varset_p.X_beta{i}=obj.stem_varset_p.X_beta{i}(indices,:,:);
                    end
                    if not(isempty(obj.stem_varset_p.X_z))
                        obj.stem_varset_p.X_z{i}=obj.stem_varset_p.X_z{i}(indices,:,:);
                    end
                    if not(isempty(obj.stem_varset_p.X_p))
                        obj.stem_varset_p.X_p{i}=obj.stem_varset_p.X_p{i}(indices,:,:,:);
                    end
                else
                    error(['    Variable ',obj.stem_varset_p.Y_name{i},' does not have sites in the crop area.']);
                end
            end
            disp('Point level data space crop ended.');
            if not(isempty(obj.stem_varset_b))
                disp('Pixel data space crop started...');
                for i=1:obj.stem_varset_b.nvar
                    GY=obj.stem_gridlist_b.grid{i}.coordinate;
                    indices = GY(:,2) >= lon_min & ...
                        GY(:,2) <= lon_max & ...
                        GY(:,1) >= lat_min & ...
                        GY(:,1) <= lat_max;
                    
                    if sum(indices)>0
                        if strcmp(obj.stem_gridlist_b.grid{i}.grid_type,'regular')
                            grid_lat=obj.stem_gridlist_b.grid{i}.coordinate(:,1);
                            grid_lon=obj.stem_gridlist_b.grid{i}.coordinate(:,2);
                            grid_lat=reshape(grid_lat,obj.stem_gridlist_b.grid{i}.grid_size);
                            grid_lat=grid_lat(:,1);
                            grid_lon=reshape(grid_lon,obj.stem_gridlist_b.grid{i}.grid_size);
                            grid_lon=grid_lon(1,:);
                            grid_lat(grid_lat<lat_min)=[];
                            grid_lat(grid_lat>lat_max)=[];
                            grid_lon(grid_lon<lon_min)=[];
                            grid_lon(grid_lon>lon_max)=[];
                        end
                       
                        obj.stem_gridlist_b.grid{i}.coordinate=obj.stem_gridlist_b.grid{i}.coordinate(indices,:);
                        if strcmp(obj.stem_gridlist_b.grid{i}.grid_type,'regular')
                            obj.stem_gridlist_b.grid{i}.grid_size=[length(grid_lat),length(grid_lon)];
                        end
                        obj.stem_varset_b.Y{i}=obj.stem_varset_b.Y{i}(indices,:);
                        if not(isempty(obj.stem_varset_b.X_bp))
                            obj.stem_varset_b.X_bp{i}=obj.stem_varset_b.X_bp{i}(indices,:,:);
                        end
                        if not(isempty(obj.stem_varset_b.X_beta))
                            obj.stem_varset_b.X_beta{i}=obj.stem_varset_b.X_beta{i}(indices,:,:);
                        end
                        if not(isempty(obj.stem_varset_b.X_z))
                            obj.stem_varset_b.X_z{i}=obj.stem_varset_b.X_z{i}(indices,:,:);
                        end
                    else
                        error(['    Variable ',obj.stem_varset_b.Y_name{i},' does not have sites in the crop area.']);
                    end
                end
                disp('Pixel data space crop ended.');
            end
            disp('Updating data matrices after space crop...');
            obj.update_data;
            disp('Update ended.');
            disp('Updating distance matrices after space crop...');
            obj.update_distance;
            disp('Update ended.');

            if not(isempty(obj.stem_varset_b))
                disp('Updating M replication vector after space crop...');
                obj.update_M;
                disp('Update ended.');
            end

        end   
        
        function site_crop(obj,type,var_name,indices)
            %DESCRIPTION: remove specific sites from the dataset
            %
            %INPUT
            %obj                 - [stem_data object]   (1x1) the stem_data object
            %type                - [string]             (1x1) 'point': remove the sites from the point dataset; 'pixel': remove the sites from the pixel dataset
            %var_name            - [string]             (1x1) the name of the variable from which to remove the sites
            %indices             - [integer >0]         (dNx1) the indices of the sites to remove
            %
            %OUTPUT
            %
            %none: the matrices Y, X_bp, X_beta, X_z and X_p with are updated
            
            if sum(strcmp(type,{'point','pixel'}))==0
                error('Type must be either point or pixel');
            end
            if min(indices<1)
                error('The minimum value of indices cannot be lower than 1');
            end
            if strcmp(type,'point')
                idx_var=obj.stem_varset_p.get_Y_index(var_name);
                if isempty(idx_var)
                    error('Variable not found');
                end
                N=obj.stem_varset_p.N;
                if max(indices>N)
                    error(['The maximum value of indices cannot be greater than ',num2str(N)]);
                end
                obj.stem_gridlist_p.grid{idx_var}.coordinate(indices,:)=[];
                obj.stem_varset_p.Y{idx_var}(indices,:)=[];
                
                if not(isempty(obj.stem_varset_p.X_bp))
                    obj.stem_varset_p.X_bp{idx_var}(indices,:,:)=[];
                end
                if not(isempty(obj.stem_varset_p.X_beta))
                    obj.stem_varset_p.X_beta{idx_var}(indices,:,:)=[];
                end
                if not(isempty(obj.stem_varset_p.X_z))
                    obj.stem_varset_p.X_z{idx_var}(indices,:,:)=[];
                end
                if not(isempty(obj.stem_varset_p.X_p))
                    obj.stem_varset_p.X_p{idx_var}(indices,:,:,:)=[];
                end
            else
                if isempty(obj.stem_varset_b)
                    error('No pixel data');
                end
                idx_var=obj.stem_varset_b.get_Y_index(var_name);
                if isempty(idx_var)
                    error('Variable not found');
                end
                obj.stem_gridlist_b.grid{idx_var}.coordinate(indices,:)=[];
                obj.stem_varset_b.Y{idx_var}(indices,:)=[];
                if not(isempty(obj.stem_varset_b.X_beta))
                    obj.stem_varset_b.X_beta{idx_var}(indices,:,:)=[];
                end
                if not(isempty(obj.stem_varset_b.X_z))
                    obj.stem_varset_b.X_z{idx_var}(indices,:,:)=[];
                end
            end
            disp('Updating data matrices after site crop...');
            obj.update_data;
            disp('Update ended.');
            disp('Updating distance matrices after site crop...');
            obj.update_distance(type);
            disp('Update ended.');
            if not(isempty(obj.stem_varset_b))
                disp('Updating M replication vector after site crop...');
                obj.update_M;
                disp('Update ended.');
            end
        end
        
        function remove_duplicated_sites(obj)
            %DESCRIPTION: remove the duplicated sites from each variable
            %
            %INPUT
            %obj        - [stem_data object]   (1x1) the stem_data object
            %
            %OUTPUT
            %
            %none: the matrices Y, X_bp, X_beta, X_z and X_p with are updated
            
            disp('Look for duplicated sites...')
            for i=1:obj.stem_varset_p.nvar
                idx=obj.stem_gridlist_p.grid{i}.duplicated_sites;
                if not(isempty(idx))
                    disp(['Removing ',num2str(length(idx)),' replicated sites for point variable ',obj.stem_varset_p.Y_name{i}]);
                    obj.site_crop('point',obj.stem_varset_p.Y_name{i},idx);
                    obj.stem_gridlist_p.grid{i}.duplicated_sites=[];
                end
            end
            
            if not(isempty(obj.stem_varset_b))
                for i=1:obj.stem_varset_b.nvar
                    idx=obj.stem_gridlist_b.grid{i}.duplicated_sites;
                    if not(isempty(idx))
                        disp(['Removing ',num2str(length(idx)),' replicated sites for pixel variable ',obj.stem_varset_b.Y_name{i}]);
                        obj.site_crop('pixel',obj.stem_varset_b.Y_name{i},idx);
                        obj.stem_gridlist_b.grid{i}.duplicated_sites=[];
                    end
                end
            end
            disp('Operation ended.')
        end        
        
        function h_fig = plot(obj,variable_name,variable_type,time_step,loading_name,loading_type)
            %DESCRIPTION: plot the variable data or the loading coefficients related to a variable
            %
            %INPUT
            %obj             - [stem_data object]   (1x1) the stem_data object
            %variable_name   - [string]             (1x1) the name of the variable
            %variable_type   - [string]             (1x1) the type of the variable, either 'point' or 'pixel'
            %time_step       - [integer >=0|string] (1x1) the time step to plot. If time_step=0 the temporal average is plotted. If time_step is a string it must be in the format dd-mm-yyyy
            %<loading_name>  - [string]             (1x1) (defalut: []) the name of the loading coefficients to plot (related to the variable specified)
            %<loading_type>  - [string]             (1x1) (default: []) the type of the loading coefficients. Can be 'beta', 'z', 'w_b' or 'w_p'
            %
            %OUTPUT
            %
            %h_fig           - [integer]            (1x1) the handle of the figure
            
            if nargin<4
                error('Not enough input arguments');
            end
            if nargin<5
                loading_name=[];
            end
            if nargin==5
                error('You must also provide the loading_type');
            end
            if ischar(time_step)
                date_num=datenum(time_step,'dd-mm-yyyy');
                time_step=find(obj.stem_datestamp.stamp==date_num);
                if isempty(time_step)
                    error(['The date stamp ',time_step,' cannot be found']);
                end
            end
            if time_step<0||time_step>obj.T
                error(['time_step out of bound. It must be between 0 and ',num2str(obj.T),' included']);
            end
            if not(strcmp(variable_type,'point')||strcmp(variable_type,'pixel'))
                error('variable_type can be either ''point'' or ''pixel''');
            end
            if nargin>=6
                if not(strcmp(loading_type,'X_beta')||strcmp(loading_type,'X_z')||strcmp(loading_type,'X_wp')||strcmp(loading_type,'X_wb'))
                    error('loading_type must be either ''X_beta'', ''X_z'', ''X_wb'' or ''X_wp''');
                end
                if strcmp(variable_type,'pixel')&&strcmp(loading_type,'X_wp')
                    error('The loading_type X_wp is not supported for pixel type data');
                end
            end
            
            if strcmp(variable_type,'point')
                index_var=obj.stem_varset_p.get_Y_index(variable_name);
                if isempty(index_var)
                    error(['The variable ',variable_name,' cannot be found as ',variable_type,' data']);
                end
                if not(isempty(loading_name))
                    if strcmp(loading_type,'X_beta')
                        indexl=obj.stem_varset_p.get_X_beta_index(loading_name,index_var);
                        if isempty(indexl)
                            error(['The loading coefficient ',loading_name,' cannot be found in X_beta for the variable ',variable_name]);
                        end
                        if not(obj.stem_varset_p.X_beta_tv)&&time_step>1
                            time_step=1;
                            disp('WARNING: the loading coefficient is time-invariant. t=1 is plotted');
                        end
                        if time_step==0
                            data=squeeze(nanmean(obj.stem_varset_p.X_beta{index_var}(:,indexl,:),3));
                        else
                            data=obj.stem_varset_p.X_beta{index_var}(:,indexl,time_step);
                        end
                    end
                    if strcmp(loading_type,'X_z')
                        indexl=obj.stem_varset_p.get_X_z_index(loading_name,index_var);
                        if isempty(indexl)
                            error(['The loading coefficient ',loading_name,' cannot be found in X_z for the variable ',variable_name]);
                        end   
                        if not(obj.stem_varset_p.X_z_tv)&&time_step>1
                            time_step=1;
                            disp('WARNING: the loading coefficient is time-invariant. t=1 is plotted');
                        end
                        if time_step==0
                            data=squeeze(nanmean(obj.stem_varset_p.X_z{index_var}(:,indexl,:),3));
                        else
                            data=obj.stem_varset_p.X_z{index_var}(:,indexl,time_step);
                        end
                    end
                    if strcmp(loading_type,'X_wb')
                        indexl=obj.stem_varset_p.get_X_bp_index(loading_name,index_var);
                        if isempty(indexl)
                            error(['The loading coefficient ',loading_name,' cannot be found in X_bp for the variable ',variable_name]);
                        end      
                        if not(obj.stem_varset_p.X_bp_tv)&&time_step>1
                            time_step=1;
                            disp('WARNING: the loading coefficient is time-invariant. t=1 is plotted');
                        end
                        if time_step==0
                            data=squeeze(nanmean(obj.stem_varset_p.X_bp{index_var}(:,indexl,:),3));
                        else
                            data=obj.stem_varset_p.X_bp{index_var}(:,indexl,time_step);
                        end
                    end
                    if strcmp(loading_type,'X_wp')
                        indexl=obj.stem_varset_p.get_X_p_index(loading_name,index_var);
                        if isempty(indexl)
                            error(['The loading coefficient ',loading_name,' cannot be found in X_p for the variable ',variable_name]);
                        end     
                        if not(obj.stem_varset_p.X_p_tv)&&time_step>1
                            time_step=1;
                            disp('WARNING: the loading coefficient is time-invariant. t=1 is plotted');
                        end
                        if time_step==0
                            data=squeeze(nanmean(obj.stem_varset_p.X_p{index_var}(:,:,:,indexl),3));
                        else
                            data=obj.stem_varset_p.X_p{index_var}(:,:,time_step,indexl);
                        end
                    end
                else
                    if time_step==0
                        data=nanmean(obj.stem_varset_p.Y{index_var},2);
                    else
                        data=obj.stem_varset_p.Y{index_var}(:,time_step);
                    end
                end
                lat=obj.stem_gridlist_p.grid{index_var}.coordinate(:,1);
                lon=obj.stem_gridlist_p.grid{index_var}.coordinate(:,2);
                if isempty(loading_name)
                    if time_step==0
                        tit=['Average ',variable_name, ' from ',datestr(obj.stem_datestamp.stamp(1)),' to ',datestr(obj.stem_datestamp.stamp(end))];
                    else
                        tit=[variable_name,' on ',datestr(obj.stem_datestamp.stamp(time_step))];
                    end
                else
                    if time_step==0
                        tit=['Average loading coefficient ',loading_name,' in ',loading_type,' for variable ',variable_name, ' from ',datestr(obj.stem_datestamp.stamp(1)),' to ',datestr(obj.stem_datestamp.stamp(end))];
                    else
                        tit=['Loading coefficient ',loading_name,' in ',loading_type,' for variable ',variable_name, ' on ',datestr(obj.stem_datestamp.stamp(time_step))];
                    end
                end
                if strcmp(obj.stem_gridlist_p.grid{index_var}.unit,'deg')
                    xlab='Longitude';
                    ylab='Latitude';
                else
                    xlab=obj.stem_gridlist_p.grid{index_var}.unit;
                    ylab=obj.stem_gridlist_p.grid{index_var}.unit;
                end
                h=stem_misc.plot_map(lat,lon,data,obj.shape,tit,xlab,ylab);
            else
                index_var=obj.stem_varset_b.get_Y_index(variable_name);
                if isempty(index_var)
                    error(['The variable ',variable_name,' cannot be found as ',variable_type,' data']);
                end
                if not(isempty(loading_name))
                    if strcmp(loading_type,'X_beta')
                        indexl=obj.stem_varset_b.get_X_beta_index(loading_name,index_var);
                        if isempty(indexl)
                            error(['The loading coefficient ',loading_name,' cannot be found in X_beta for the variable ',variable_name]);
                        end
                        if not(obj.stem_varset_b.X_beta_tv)&&time_step>1
                            time_step=1;
                            disp('WARNING: the loading coefficient is time-invariant. t=1 is plotted');
                        end
                        if time_step==0
                            data=squeeze(nanmean(obj.stem_varset_b.X_beta{index_var}(:,indexl,:),3));
                        else
                            data=obj.stem_varset_b.X_beta{index_var}(:,indexl,time_step);
                        end                        
                    end
                    if strcmp(loading_type,'X_z')
                        indexl=obj.stem_varset_b.get_X_z_index(loading_name,index_var);
                        if isempty(indexl)
                            error(['The loading coefficient ',loading_name,' cannot be found in X_z for the variable ',variable_name]);
                        end
                        if not(obj.stem_varset_b.X_z_tv)&&time_step>1
                            time_step=1;
                            disp('WARNING: the loading coefficient is time-invariant. t=1 is plotted');
                        end
                        if time_step==0
                            data=squeeze(nanmean(obj.stem_varset_b.X_z{index_var}(:,indexl,:),3));
                        else
                            data=obj.stem_varset_b.X_z{index_var}(:,indexl,time_step);
                        end                        
                    end
                    if strcmp(loading_type,'X_wb')
                        indexl=obj.stem_varset_b.get_X_bp_index(loading_name,index_var);
                        if isempty(indexl)
                            error(['The loading coefficient ',loading_name,' cannot be found in X_bp for the variable ',variable_name]);
                        end
                        if not(obj.stem_varset_b.X_bp_tv)&&time_step>1
                            time_step=1;
                            disp('WARNING: the loading coefficient is time-invariant. t=1 is plotted');
                        end
                        if time_step==0
                            data=squeeze(nanmean(obj.stem_varset_b.X_bp{index_var}(:,indexl,:),3));
                        else
                            data=obj.stem_varset_b.X_bp{index_var}(:,indexl,time_step);
                        end                        
                    end
                else
                    if time_step==0
                        data=nanmean(obj.stem_varset_b.Y{index_var},2);
                    else
                        data=obj.stem_varset_b.Y{index_var}(:,time_step);
                    end
                end
                lat=obj.stem_gridlist_b.grid{index_var}.coordinate(:,1);
                lon=obj.stem_gridlist_b.grid{index_var}.coordinate(:,2);
                lat=reshape(lat,obj.stem_gridlist_b.grid{index_var}.grid_size);
                lon=reshape(lon,obj.stem_gridlist_b.grid{index_var}.grid_size);
                data=reshape(data,obj.stem_gridlist_b.grid{index_var}.grid_size);
                if isempty(loading_name)
                    if time_step==0
                        tit=['Average ',variable_name, ' from ',datestr(obj.stem_datestamp.stamp(1)),' to ',datestr(obj.stem_datestamp.stamp(end))];
                    else
                        tit=[variable_name,' on ',datestr(obj.stem_datestamp.stamp(time_step))];
                    end
                else
                    if time_step==0
                        tit=['Average loading coefficient ',loading_name,' in ',loading_type,' for variable ',variable_name, ' from ',datestr(obj.stem_datestamp.stamp(1)),' to ',datestr(obj.stem_datestamp.stamp(end))];
                    else
                        tit=['Loading coefficient ',loading_name,' in ',loading_type,' for variable ',variable_name, ' on ',datestr(obj.stem_datestamp.stamp(time_step))];
                    end
                end                
                if strcmp(obj.stem_gridlist_b.grid{index_var}.unit,'deg')
                    xlab='Longitude';
                    ylab='Latitude';
                else
                    xlab=obj.stem_gridlist_p.grid{index_var}.unit;
                    ylab=obj.stem_gridlist_p.grid{index_var}.unit;
                end                
                h=stem_misc.plot_map(lat,lon,data,obj.shape,tit,xlab,ylab);
            end
            if nargout>0
                h_fig=h;
            end
        end        
        
        function print(obj)
            %DESCRIPTION: print the information on data and their structure
            %
            %INPUT
            %obj  - [stem_data object]   (1x1) the stem_data object
            %
            %OUTPUT
            %
            %none: the information is printed in the command window
            
            stem_misc.disp_star('Data description');
            disp(['Number of variables: ',num2str(obj.nvar)]);
            disp(' ');
            disp('Point variables:');
            for i=1:length(obj.stem_varset_p.Y_name)
                disp(['  (',num2str(i),') - ',obj.stem_varset_p.Y_name{i}]);
            end
            disp(' '); 
            if not(isempty(obj.stem_varset_b))
                disp('Pixel variables:');
                for i=1:length(obj.stem_varset_b.Y_name)
                    disp(['  (',num2str(i),') - ',obj.stem_varset_b.Y_name{i}]);
                end
                disp(' ');
            end
            disp('Date and time steps:')
            disp(['  Stating date  : ',datestr(obj.stem_datestamp.date_start)]);
            disp(['  Ending date   : ',datestr(obj.stem_datestamp.date_end)]);
            disp(['  Temporal steps: ',num2str(obj.T)]);
            disp(' ');    
            
            disp('Bounding box of the point data:');
            if strcmp(obj.stem_gridlist_p.grid{1}.unit,'deg')
                prefix_x='  Longitude ';
                prefix_y='  Latitude ';
                postfix='°';
            else
                prefix_x='X ';
                prefix_y='Y ';
                postfix=[' ',obj.stem_gridlist_p.grid{1}.unit];
            end
            disp([prefix_y,'min : ',num2str(obj.stem_gridlist_p.box(1),'%05.2f'),postfix])
            disp([prefix_y,'max : ',num2str(obj.stem_gridlist_p.box(2),'%05.2f'),postfix])
            disp([prefix_x,'min: ',num2str(obj.stem_gridlist_p.box(3),'%05.2f'),postfix])
            disp([prefix_x,'max: ',num2str(obj.stem_gridlist_p.box(4),'%05.2f'),postfix])
            disp(' ');
            if not(isempty(obj.stem_varset_b))
                disp('Bounding box of the pixel data:');
                disp([prefix_y,'min : ',num2str(obj.stem_gridlist_b.box(1),'%05.2f'),postfix])
                disp([prefix_y,'max : ',num2str(obj.stem_gridlist_b.box(2),'%05.2f'),postfix])
                disp([prefix_x,'min: ',num2str(obj.stem_gridlist_b.box(3),'%05.2f'),postfix])
                disp([prefix_x,'max: ',num2str(obj.stem_gridlist_b.box(4),'%05.2f'),postfix])
                disp(' ');
            end

            disp('Variable description - Point data');
            output{1,1}='Name';
            output{1,2}='#sites';
            output{1,3}='Mean';            
            output{1,4}='Std';            
            output{1,5}='Min';            
            output{1,6}='Max';            
            output{1,7}='Missing';            
            for i=1:obj.stem_varset_p.nvar
                output{i+1,1}=obj.stem_varset_p.Y_name{i};
                output{i+1,2}=num2str(obj.stem_varset_p.dim(i),'%03.0f');
                output{i+1,3}=num2str(nanmean(stem_misc.vec(obj.stem_varset_p.Y{i})),'%+05.2f');
                output{i+1,4}=num2str(nanstd(stem_misc.vec(obj.stem_varset_p.Y{i})),'%05.2f');
                output{i+1,5}=num2str(nanmin(stem_misc.vec(obj.stem_varset_p.Y{i})),'%+05.2f');
                output{i+1,6}=num2str(nanmax(stem_misc.vec(obj.stem_varset_p.Y{i})),'%+05.2f');
                output{i+1,7}=[num2str(sum(isnan(stem_misc.vec(obj.stem_varset_p.Y{i})))/length(stem_misc.vec(obj.stem_varset_p.Y{i}))*100,'%05.2f'),'%'];
            end
            disp(output);
            output=[];
            if not(isempty(obj.stem_varset_b))
                disp('Variable description - Pixel data');
                output{1,1}='Name';
                output{1,2}='#sites';
                output{1,3}='Mean';
                output{1,4}='Std';
                output{1,5}='Min';
                output{1,6}='Max';
                output{1,7}='Missing';
                for i=1:obj.stem_varset_b.nvar
                    output{i+1,1}=obj.stem_varset_b.Y_name{i};
                    output{i+1,2}=num2str(obj.stem_varset_b.dim(i),'%03.0f');
                    output{i+1,3}=num2str(nanmean(stem_misc.vec(obj.stem_varset_b.Y{i})),'%+05.2f');
                    output{i+1,4}=num2str(nanstd(stem_misc.vec(obj.stem_varset_b.Y{i})),'%05.2f');
                    output{i+1,5}=num2str(nanmin(stem_misc.vec(obj.stem_varset_b.Y{i})),'%+05.2f');
                    output{i+1,6}=num2str(nanmax(stem_misc.vec(obj.stem_varset_b.Y{i})),'%+05.2f');
                    output{i+1,7}=[num2str(sum(isnan(stem_misc.vec(obj.stem_varset_b.Y{i})))/length(stem_misc.vec(obj.stem_varset_b.Y{i}))*100,'%05.2f'),'%'];
                end
            end
            disp(output);
            output=[];
            for i=1:obj.stem_varset_p.nvar
                stem_misc.disp_star(['Loading coefficients of the point variable ',obj.stem_varset_p.Y_name{i}]);
                if not(isempty(obj.stem_varset_p.X_beta_name))
                    disp('Loading coefficients related to beta:');
                    output{1,1}='Name';
                    output{1,2}='Mean';
                    output{1,3}='Std';
                    output{1,4}='Min';
                    output{1,5}='Max';
                    for j=1:length(obj.stem_varset_p.X_beta_name{i})
                        output{j+1,1}=obj.stem_varset_p.X_beta_name{i}{j};
                        output{j+1,2}=num2str(nanmean(stem_misc.vec(obj.stem_varset_p.X_beta{i}(:,j,:))),'%+05.2f');
                        output{j+1,3}=num2str(nanstd(stem_misc.vec(obj.stem_varset_p.X_beta{i}(:,j,:))),'%05.2f');
                        output{j+1,4}=num2str(nanmin(stem_misc.vec(obj.stem_varset_p.X_beta{i}(:,j,:))),'%+05.2f');
                        output{j+1,5}=num2str(nanmax(stem_misc.vec(obj.stem_varset_p.X_beta{i}(:,j,:))),'%+05.2f');
                    end
                    disp(output);
                end
                output=[];
                if not(isempty(obj.stem_varset_p.X_z_name))
                    disp('Loading coefficients related to the latent variable z(t):');
                    output{1,1}='Name';
                    output{1,2}='Mean';
                    output{1,3}='Std';
                    output{1,4}='Min';
                    output{1,5}='Max';                    
                    for j=1:length(obj.stem_varset_p.X_z_name{i})
                        output{j+1,1}=obj.stem_varset_p.X_z_name{i}{j};
                        output{j+1,2}=num2str(nanmean(stem_misc.vec(obj.stem_varset_p.X_z{i}(:,j,:))),'%+05.2f');
                        output{j+1,3}=num2str(nanstd(stem_misc.vec(obj.stem_varset_p.X_z{i}(:,j,:))),'%05.2f');
                        output{j+1,4}=num2str(nanmin(stem_misc.vec(obj.stem_varset_p.X_z{i}(:,j,:))),'%+05.2f');
                        output{j+1,5}=num2str(nanmax(stem_misc.vec(obj.stem_varset_p.X_z{i}(:,j,:))),'%+05.2f');
                    end
                    disp(output);
                end
                output=[];
                if not(isempty(obj.stem_varset_p.X_p_name))
                    disp('Loading coefficients related to the latent variable w_p(s,t):');
                    output{1,1}='Name';
                    output{1,2}='Mean';
                    output{1,3}='Std';
                    output{1,4}='Min';
                    output{1,5}='Max';                       
                    for j=1:length(obj.stem_varset_p.X_p_name{i})
                        output{j+1,1}=obj.stem_varset_p.X_p_name{i}{j};
                        output{j+1,2}=num2str(nanmean(stem_misc.vec(obj.stem_varset_p.X_p{i}(:,:,:,j))),'%+05.2f');
                        output{j+1,3}=num2str(nanstd(stem_misc.vec(obj.stem_varset_p.X_p{i}(:,:,:,j))),'%05.2f');
                        output{j+1,4}=num2str(nanmin(stem_misc.vec(obj.stem_varset_p.X_p{i}(:,:,:,j))),'%+05.2f');
                        output{j+1,5}=num2str(nanmax(stem_misc.vec(obj.stem_varset_p.X_p{i}(:,:,:,j))),'%+05.2f');
                    end
                    disp(output);
                end
                output=[];
                if not(isempty(obj.stem_varset_p.X_bp_name))
                    disp('Loading coefficients related to the latent variable w_b(s,t):');
                    output{1,1}='Name';
                    output{1,2}='Mean';
                    output{1,3}='Std';
                    output{1,4}='Min';
                    output{1,5}='Max';                      
                    for j=1:length(obj.stem_varset_p.X_bp_name{i})
                        output{j+1,1}=obj.stem_varset_p.X_bp_name{i}{j};
                        output{j+1,2}=num2str(nanmean(stem_misc.vec(obj.stem_varset_p.X_bp{i}(:,j,:))),'%+05.2f');
                        output{j+1,3}=num2str(nanstd(stem_misc.vec(obj.stem_varset_p.X_bp{i}(:,j,:))),'%05.2f');
                        output{j+1,4}=num2str(nanmin(stem_misc.vec(obj.stem_varset_p.X_bp{i}(:,j,:))),'%+05.2f');
                        output{j+1,5}=num2str(nanmax(stem_misc.vec(obj.stem_varset_p.X_bp{i}(:,j,:))),'%+05.2f');
                    end
                    disp(output);
                end
                output=[];
            end
                
            if not(isempty(obj.stem_varset_b))
                for i=1:obj.stem_varset_b.nvar
                    stem_misc.disp_star(['Loading coefficients of the pixel variable ',obj.stem_varset_b.Y_name{i}]);
                    if not(isempty(obj.stem_varset_b.X_beta_name))
                        disp('Loading coefficients related to beta:');
                        output{1,1}='Name';
                        output{1,2}='Mean';
                        output{1,3}='Std';
                        output{1,4}='Min';
                        output{1,5}='Max';
                        for j=1:length(obj.stem_varset_b.X_beta_name{i})
                            output{j+1,1}=obj.stem_varset_b.X_beta_name{i}{j};
                            output{j+1,2}=num2str(nanmean(stem_misc.vec(obj.stem_varset_b.X_beta{i}(:,j,:))),'%+05.2f');
                            output{j+1,3}=num2str(nanstd(stem_misc.vec(obj.stem_varset_b.X_beta{i}(:,j,:))),'%05.2f');
                            output{j+1,4}=num2str(nanmin(stem_misc.vec(obj.stem_varset_b.X_beta{i}(:,j,:))),'%+05.2f');
                            output{j+1,5}=num2str(nanmax(stem_misc.vec(obj.stem_varset_b.X_beta{i}(:,j,:))),'%+05.2f');
                        end
                        disp(output);
                    end
                    output=[];
                    if not(isempty(obj.stem_varset_b.X_z_name))
                        disp('Loading coefficients related to the latent variable z(t):');
                        output{1,1}='Name';
                        output{1,2}='Mean';
                        output{1,3}='Std';
                        output{1,4}='Min';
                        output{1,5}='Max';
                        for j=1:length(obj.stem_varset_b.X_z_name{i})
                            output{j+1,1}=obj.stem_varset_b.X_z_name{i}{j};
                            output{j+1,2}=num2str(nanmean(stem_misc.vec(obj.stem_varset_b.X_z{i}(:,j,:))),'%+05.2f');
                            output{j+1,3}=num2str(nanstd(stem_misc.vec(obj.stem_varset_b.X_z{i}(:,j,:))),'%05.2f');
                            output{j+1,4}=num2str(nanmin(stem_misc.vec(obj.stem_varset_b.X_z{i}(:,j,:))),'%+05.2f');
                            output{j+1,5}=num2str(nanmax(stem_misc.vec(obj.stem_varset_b.X_z{i}(:,j,:))),'%+05.2f');
                        end
                        disp(output);
                    end
                    output=[];
                    if not(isempty(obj.stem_varset_b.X_bp_name))
                        disp('Loading coefficients related to the latent variable w_b(s,t):');
                        output{1,1}='Name';
                        output{1,2}='Mean';
                        output{1,3}='Std';
                        output{1,4}='Min';
                        output{1,5}='Max';
                        for j=1:length(obj.stem_varset_b.X_bp_name{i})
                            output{j+1,1}=obj.stem_varset_b.X_bp_name{i}{j};
                            output{j+1,2}=num2str(nanmean(stem_misc.vec(obj.stem_varset_b.X_bp{i}(:,j,:))),'%+05.2f');
                            output{j+1,3}=num2str(nanstd(stem_misc.vec(obj.stem_varset_b.X_bp{i}(:,j,:))),'%05.2f');
                            output{j+1,4}=num2str(nanmin(stem_misc.vec(obj.stem_varset_b.X_bp{i}(:,j,:))),'%+05.2f');
                            output{j+1,5}=num2str(nanmax(stem_misc.vec(obj.stem_varset_b.X_bp{i}(:,j,:))),'%+05.2f');
                        end
                        disp(output);
                    end
                end
            end
        end
        
        %Export methods
        function N = N(obj)
            N=obj.stem_varset_p.N;
            if not(isempty(obj.stem_varset_b))
                N=N+obj.stem_varset_b.N;
            end
        end
        
        function Nb = Nb(obj)
            if not(isempty(obj.stem_varset_b))
                Nb = obj.stem_varset_b.N;
            else
                Nb=0;
            end
        end
        
        function Np = Np(obj)
            Np = obj.stem_varset_p.N;
        end
        
        function T = T(obj)
            T=obj.stem_varset_p.T;
        end
        
        function nvar = nvar(obj)
            nvar=obj.stem_varset_p.nvar;
            if not(isempty(obj.stem_varset_b))
                nvar=nvar+obj.stem_varset_b.nvar;
            end
        end
        
        function dim = dim(obj)
            dim=obj.stem_varset_p.dim;
            if not(isempty(obj.stem_varset_b))
                dim=[dim obj.stem_varset_b.dim];
            end
        end
        
        %Class set methods
        function set.stem_varset_p(obj,stem_varset_p)
           if not(isa(stem_varset_p,'stem_varset'))
               error('stem_varset must be of class stem_varset');
           end
           obj.stem_varset_p=stem_varset_p;
        end
        
        function set.stem_varset_b(obj,stem_varset_b)
           if not(isa(stem_varset_b,'stem_varset'))
               error('stem_varset must be of class stem_varset');
           end
           
           if not(length(stem_varset_b.dim)==length(obj.stem_varset_p.dim))
               error('stem_varset_b must contain the same number of variables of stem_varset_p');
           end
           
           if not(size(stem_varset_b.Y{1},2)==size(obj.stem_varset_p.Y{1},2))
               error('The number of temporal steps cannot differ between stem_varset_b and stem_varset_p');
           end
           
           if not(isempty(stem_varset_b.X_p))
               error('X_p must be empty in stem_varset_b');
           end
          
           obj.stem_varset_b=stem_varset_b;
        end        
        
        function set.stem_gridlist_p(obj,stem_gridlist_p)
            if not(isa(stem_gridlist_p,'stem_gridlist'))
                error('stem_gridlist must be of class stem_gridlist');
            end
            if not(length(stem_gridlist_p.grid)==length(obj.stem_varset_p.Y))
                error('The number of stem_grids must be equal to the q');
            end
            for i=1:length(stem_gridlist_p.grid)
                if not(size(stem_gridlist_p.grid{i}.coordinate,1)==size(obj.stem_varset_p.Y{i},1))
                    error('The number of coordinates in the grid{i} must be equal to the number of rows of Y{i}');
                end
                if not(strcmp(stem_gridlist_p.grid{i}.site_type,'point'))
                    error('Only point data are supported in stem_gridlist_p');
                end
            end
            obj.stem_gridlist_p=stem_gridlist_p;
        end
        
        function set.stem_gridlist_b(obj,stem_gridlist_b)
            if not(isa(stem_gridlist_b,'stem_gridlist'))
                error('stem_gridlist must be of class stem_gridlist');
            end
            if not(length(stem_gridlist_b.grid)==length(obj.stem_varset_b.Y))
                error('The number of stem_grids must be equal to the q');
            end            
            for i=1:length(stem_gridlist_b.grid)
                if not(size(stem_gridlist_b.grid{i}.coordinate,1)==size(obj.stem_varset_b.Y{i},1))
                    error('The number of coordinates in the grid{i} must be equal to the number of rows of Y{i}');
                end
                if not(strcmp(stem_gridlist_b.grid{i}.site_type,'pixel'))
                    error('The grids of stem_gridlist_b must be grids of pixels');
                end
                if not(strcmp(stem_gridlist_b.grid{i}.pixel_shape,'square'))
                    error('Only square pixels are supported. Check pixel shape');
                end
                if not(stem_gridlist_b.grid{i}.pixel_side_w==stem_gridlist_b.grid{i}.pixel_side_h)
                    error('Only square pixels are supported. Check pixel_side_w and pixel_side_h');
                end
            end 
            if not(strcmp(stem_gridlist_b.grid{1}.unit,obj.stem_gridlist_p.grid{1}.unit))
                error('Both the stem_gridlist objects must contain grids with the same unit');
            end
            obj.stem_gridlist_b=stem_gridlist_b;
        end  

        function set.stem_crossval(obj,stem_crossval)
            if not(isa(stem_crossval,'stem_crossval'))
                error('stem_crossval must be of class stem_crossval');
            end
            obj.stem_crossval=stem_crossval;
        end
        
        function set.stem_datestamp(obj,stem_datestamp)
            if not(isa(stem_datestamp,'stem_datestamp'))
                error('stem_datestamp must be a stem_datestamp object');
            end
            if not(obj.stem_varset_p.T==length(stem_datestamp.stamp))
                error('The number of datestamps differs from T');
            end
            obj.stem_datestamp=stem_datestamp;
        end
        
        function set.shape(obj,shape)
            obj.shape=shape;
        end
        
    end
end

%         function google_map(obj,name,type)
%             if sum(strcmp(type,{'point','pixel'}))==0
%                 error('type must be either point or pixel');
%             end
%             if strcmp(type,'point')
%                 %not supported yet
%             else
%                 if isempty(obj.stem_varset_b)
%                     disp('No pixel variables in this model');
%                 else
%                     idx=find(strcmp(name,obj.stem_varset_b.Y_name));
%                     T=size(obj.stem_varset_b.Y{idx},2);
%                     datafile=[1,T,obj.stem_gridlist_b.grid{idx}.pixel_side_w];
%                     csvwrite('..\Data\google_bridge\parameters_kriging.csv',datafile);
%                     for t=1:T
%                         t
%                         min_value=nanmin(nanmin(obj.stem_varset_b.Y{idx}(:,t)));
%                         max_value=nanmax(nanmax(obj.stem_varset_b.Y{idx}(:,t)));
%                         data=obj.stem_varset_b.Y{idx}(:,t);
%                         value=(data-min_value)/(max_value-min_value);
%                         %color=stem_krig_result.toColor(value);
%                         datafile=[obj.stem_gridlist_b.grid{idx}.coordinate(:,1),obj.stem_gridlist_b.grid{idx}.coordinate(:,2),value,data];
%                         csvwrite(['..\Data\google_bridge\kriging',num2str(t),'.csv'],datafile);
%                     end
%                     winopen('..\Data\google_bridge\open_kriging.bat');
%                 end
%             end
%         end
       
        
%         function reset(obj)
%             % restore original data
%             if obj.can_reset
%                 disp('Restoring original data...');
% 
%                 disp('Original data restored.');
%             else
%                 warning('Reset is not allow. Data must be reloaded by using the class constructor');
%             end
%         end
        
%         function [v,theta,sigma_eps] = variogram(obj,n_steps,time_steps,graph)
%             if nargin<2
%                 n_steps=50;
%             end
%             if nargin<3
%                 time_steps=1:obj.stem_varset.T;
%             end
%             if nargin<4
%                 graph=0;
%             end
%             if isempty(time_steps)
%                 time_steps=1:obj.stem_varset.T;
%             end
%             options = optimset('Display','off');
%             
%             v=zeros(length(obj.stem_varset.dim));
%             theta=zeros(length(obj.stem_varset.dim));
%             sigma_eps=zeros(length(obj.stem_varset.dim));
%             
%             blocks=[0 cumsum(obj.stem_varset.dim)];
%             for k=1:obj.stem_varset.nvar
%                 for h=k:obj.stem_varset.nvar 
%                     y_sub1=obj.stem_varset.Y{h};
%                     y_sub2=obj.stem_varset.Y{k};
%                     DistMat_sub=obj.DistMat(blocks(h)+1:blocks(h+1),blocks(k)+1:blocks(k+1));
%                     max_distance=max(DistMat_sub(:));
%                     step=max_distance/n_steps;
%                     groups=0:step:max_distance;
%                     group=zeros(length(groups)-1,1);
%                     count=zeros(length(groups)-1,1);
%                     
%                     if h~=k
%                         tot=obj.stem_varset.dim(k)*obj.stem_varset.dim(h)*length(time_steps);
%                     else
%                         tot=((obj.stem_varset.dim(k)*(obj.stem_varset.dim(k)+1))/2)*length(time_steps);
%                     end
%                     vx=zeros(tot,1);
%                     vy=zeros(tot,1);
%                     index=1;
%                     for t=time_steps
%                         if h~=k
%                             y1=y_sub1(:,t);
%                             y2=y_sub2(:,t);
%                             for i=1:length(y1)
%                                 d=DistMat_sub(i,:);
%                                 sq=0.5*(y1(i)-y2).^2;
%                                 vx(index:index+length(d)-1)=d;
%                                 vy(index:index+length(d)-1)=sq;
%                                 index=index+length(d);
%                             end
%                         else
%                             y1=y_sub1(:,t);
%                             y2=y_sub2(:,t);
%                             for i=1:length(y1)
%                                 d=DistMat_sub(i,i:end);
%                                 sq=0.5*(y1(i)-y2(i:end)).^2;
%                                 vx(index:index+length(d)-1)=d;
%                                 vy(index:index+length(d)-1)=sq;
%                                 index=index+length(d);
%                             end
%                         end
%                     end
%                             
%                     for i=2:length(groups)
%                         data=vy(vx>=groups(i-1)&vx<groups(i));
%                         group(i-1)=nanmean(data);
%                         count(i-1)=length(data(isnotnan(data)));
%                     end
%                     xv=(groups(1:10)+(step/2))';
%                     yv=group(1:10);
%                     xv=xv(not(isnan(yv)));
%                     yv=yv(not(isnan(yv)));
%                     b=regress(yv,[ones(length(xv),1) xv]);
%                     v0=b(1);
%                     if v0<0
%                         v0=0.7;
%                     end
%                             
%                     L=isnotnan(group);
%                     L(1)=0;
%                     g=groups(1:end-1)+(step/2);
%                     x1=g(L)';
%                     y1=group(L);
%                     weights=count(L);
%                                                 
%                     if 0
%                         x0=[1 20];
%                         f = @(x,xdata) x(1)*(1-exp(-xdata/x(2)))+v0;'x';'xdata';
%                         c = lsqcurvefit(f,x0,x1,y1,[],[],options);
%                         theta0(index_m)=c(2);
%                         vv(index_m)=c(1);
%                         index_m=index_m+1;
%                     else
%                         o = fitoptions('Method','NonlinearLeastSquares');
%                         o.Weights=weights;
%                         o.StartPoint=[1 max_distance/10 v0];
%                         f = fittype('a*(1-exp(-x/b))+c');
%                         try
%                             %fit1 = fit(x1,y1-v0,f,o);
%                             fit1 = fit(x1,y1,f,o);
%                             c=coeffvalues(fit1);
%                         catch
%                             disp('Error');
%                         end
%                     end
%                             
%                     if graph
%                         figure
%                         xx=0:0.1:max(groups);
%                         yy=c(1)*(1-exp(-xx./c(2)))+c(3);
%                         plot(groups(1:end-1)+(step/2),group,'.');
%                         hold on
%                         plot(xx,yy,'-r');
%                         xlabel('Distance');
%                         if h==k
%                             title(['Variogram of ',obj.stem_varset.name{k}]);
%                         else
%                             title(['Cross-variogram between ',obj.stem_varset.name{k},' and ',obj.stem_varset.name{h}]);
%                         end
%                     end
%                 end
%             end
%         end
        

%         function temporal_cross_plot(obj,variable_name1,variable_name2,site,estimated)
%             if nargin<3
%                 error('Not enough input parameters');
%             end
%             var1_index=obj.stem_varset.get_index(variable_name1);
%             var2_index=obj.stem_varset.get_index(variable_name2);
%             if isempty(var1_index)
%                 error('The first variable_name is incorrect');
%             end
%             if isempty(var2_index)
%                 error('The second variable_name is incorrect');
%             end
%             if (nargin<4)||(isempty(site))
%                 site=1:length(obj.stem_gridlist.grid{var1_index}.coordinate);
%             end
%             if nargin<5
%                 estimated=0;
%             end
%             if estimated&&(isempty(obj.Y_hat))
%                 error('The model has not been estimated yet');
%             end
%             blocks=[0 cumsum(obj.stem_varset.dim)];
%             figure
%             done=0;
%             for j=site
%                 latlon1=obj.stem_gridlist.grid{var1_index}.coordinate(j,:);
%                 for i=1:length(obj.stem_gridlist.grid{var2_index}.coordinate)
%                     latlon2=obj.stem_gridlist.grid{var2_index}.coordinate(i,:);
%                     if sum(latlon1-latlon2)==0
%                         done=done+1;
%                         index1=blocks(var1_index)+j;
%                         index2=blocks(var2_index)+i;
%                         if estimated==0
%                             data1=obj.Y(index1,:);
%                             data2=obj.Y(index2,:);
%                         else
%                             data1=obj.Y_hat(index1,:);
%                             data2=obj.Y_hat(index2,:);
%                         end
%                         for t=1:length(data1)
%                             plot(data1(t),data2(t),'o','MarkerFaceColor',[t/length(data1) t/length(data1) t/length(data1)]);
%                             hold on
%                         end
%                     end
%                 end
%             end
%             if length(site)==1
%                 title(['Cross plot between ',variable_name1,' and ',variable_name2,' at site number',num2str(site)]);
%             else
%                 title(['Cross plot between ',variable_name1,' and ',variable_name2,' at ',num2str(done),' sites']);
%             end
%             xlabel(variable_name1);
%             ylabel(variable_name2);
%             if not(done)
%                 error('The two variables are not observed at the same site');
%             end
%         end
        
%         function v = temporal_cross_correlation(obj,type,graph)
%             % temporal cross-correlation between variables at sites. The variable with the lower number of sites is the reference variable
%             % graph: 0: no graph is plotted; 1: time-series graph and correlation graph are plotted
%             if nargin<2
%                 type='nearest';
%             end
%             if nargin<3
%                 graph=0;
%             end
%             if not(strcmp(type,'nearest')||strcmp(type,'colocated'))
%                 error('type must be either nearest or colocated');
%             end
%             v_matrix=eye(obj.stem_varset.nvar);
%             for j=1:obj.stem_varset.nvar-1
%                 for i=j+1:obj.stem_varset.nvar
%                     if obj.stem_varset.dim(i)<obj.stem_varset.dim(j)
%                         h=i;
%                         k=j;
%                         swtc=false;
%                     else
%                         h=j;
%                         k=i;
%                         swtc=true;
%                     end
%                     if graph
%                         figure
%                         l=obj.stem_varset.dim(h)^0.5;
%                         if round(l)^2==obj.stem_varset.dim(h)
%                             rows=l;
%                             cols=l;
%                         else
%                             rows=ceil(l);
%                             cols=round(l);
%                         end
%                     end
%                     ok_count=0;
%                     d=[];
%                     correlation=[];
%                     weight=[];
%                     for z=1:obj.stem_varset.dim(h)
%                         d(:,1)=abs(obj.stem_gridlist.grid{k}.coordinate(:,1)-obj.stem_gridlist.grid{h}.coordinate(z,1));
%                         d(:,2)=abs(obj.stem_gridlist.grid{k}.coordinate(:,2)-obj.stem_gridlist.grid{h}.coordinate(z,2));
%                         d=sum(d,2);
%                         [m,index]=min(d);
%                         if (strcmp(type,'colocated'))&&(sum(m)~=0)
%                             ok=false;
%                         else
%                             ok=true;
%                         end
%                         ok=true;
%                         if ok
%                             ok_count=ok_count+1;
%                             temp1=obj.stem_varset.Y{h}(z,:);
%                             temp2=obj.stem_varset.Y{k}(index,:);
%                             c=nancov(temp1,temp2);
%                             correlation(ok_count)=c(1,2)/(c(1,1)*c(2,2))^0.5;
%                             weight(ok_count)=sum(isnotnan(temp1)&isnotnan(temp2))/obj.stem_varset.T;
%                             if graph
%                                 subplot(rows,cols,z);
%                                 plot((temp1-nanmean(temp1))/nanstd(temp1),'b*');
%                                 hold on
%                                 plot((temp2-nanmean(temp2))/nanstd(temp2),'r*');
%                                 title(['Correlation: ',num2str(correlation(z))]);
%                             end
%                         end
%                     end
% %                     if graph
% %                         legend(obj.variable_name{h},obj.variable_name{k},'Location','EastOutside');
% %                     end
%                     L=isnotnan(correlation);
%                     correlation=correlation(L);
%                     weight=weight(L);
%                     if graph
%                         figure
%                         subplot(1,2,1);
%                         hist(correlation);
%                         xlabel('Correlation');
%                         subplot(1,2,2);
%                         plot((1-weight)*100,correlation,'*');
%                         xlabel('% full time-period missing');
%                         ylabel('Correlation');
%                     end
%                     mean_correlation=(correlation*weight')/sum(weight);
%                     v_matrix(h,k)=mean_correlation;
%                     v_matrix(k,h)=mean_correlation;
%                     if swtc
%                         disp(['Temporal mean weighted correlation between ',obj.stem_varset.name{h},' and ',obj.stem_varset.name{k},': ',num2str(mean_correlation), '  (',num2str(ok_count),' sites)']);
%                     else
%                         disp(['Temporal mean weighted correlation between ',obj.stem_varset.name{k},' and ',obj.stem_varset.name{h},': ',num2str(mean_correlation), '  (',num2str(ok_count),' sites)']);
%                     end
%                 end
%             end
%             if nargout>0
%                 v=v_matrix;
%                 if min(eig(v))<0
%                     warning('The output correlation matrix is not semi-definite positive!')
%                 end
%             end
%         end
        
%         function [matrix,average] = temporal_correlation(obj,variable_name)
%             var_index=obj.stem_varset.get_index(variable_name);    
%             if isempty(var_index)
%                 error('The variable name is incorrect')
%             end
%             data=obj.stem_varset.Y{var_index};
%             matrix=eye(size(data,1));
%             average=0;
%             counter=0;
%             for j=1:size(data,1)
%                 for i=j+1:size(data,1)
%                     a=data(i,:);
%                     b=data(j,:);
%                     L1=isnotnan(a);
%                     L2=isnotnan(b);
%                     L=L1&L2;
%                     a=a(L);
%                     b=b(L);
%                     temp=cov(a,b);
%                     matrix(i,j)=temp(1,2)/(var(a)*var(b))^0.5;
%                     matrix(j,i)=matrix(i,j);
%                     if isnotnan(matrix(i,j))
%                         average=average+matrix(i,j)*sum(L);
%                         counter=counter+sum(L);
%                     end
%                 end
%             end
%             average=average/counter;
%         end
        
        


