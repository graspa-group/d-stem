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

classdef stem_model < handle
    
    %CONSTANTS
    %N  = n1_p+...+nq_p+n1_b+...+nq_b - total number of observation sites
    %N_p = n1_p+...+nq_p - total number of point sites
    %N_b = n1_b+...+nq_b - total number of pixel sites
    %S = 1 if only point data are considered and S=2 if both point and pixel data are considered
    
    properties
        stem_data=[];           %[stem_data object] (1x1) object containing all the data used to estimated the model
        stem_par=[];            %[stem_par object]  (1x1) parameter set updated at each iteration of the EM algorithm
        stem_par_initial=[];    %[stem_par object]  (1x1) starting parameter set
        stem_par_sim=[];        %[stem_par object]  (1x1) parameter set used to simulate data (if data are simulated)
        estimated=0;            %[boolean] (1x1) 0: the model has not been estimated; 1: the model has been estimated
    end
    
    properties (SetAccess = private)
        stem_EM_result=[];      %[stem_EM_result object] (1x1) object containing all the results of the EM estimation
        cross_validation=0;     %[boolean] (1x1) 0: the model has been estimated considering all the data; 1: the model has bee estimated excluding the cross-validation data.
        system_size=100;        %[integer] (1x1) if N > system_size than only the diagonal is computed in matrix multiply operations
        tapering=[];            %[boolean] (1x1) 0:tapering is not enabled; 1:tapering is enabled on point sites or pixel sites
        tapering_b=[];          %[boolean] (1x1) 0:tapering is not enabled on pixel sites; 1:tapering is enabled on pixel sites
        tapering_p=[];          %[boolean] (1x1) 0:tapering is not enabled on point sites; 1:tapering is enabled on point sites
    end
    
    methods
        
        function obj = stem_model(stem_data,stem_par)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %stem_data - [stem_data object]  (1x1)
            %stem_par  - [stem_par object]   (1x1)
            %
            %OUTPUT
            %obj       - [stem_model object] (1x1)
            if nargin<2
                error('Not enough input arguments');
            end
            obj.stem_data=stem_data;
            obj.stem_par=stem_par;
            
            if not(isempty(stem_data.stem_gridlist_b))
                if not(isempty(stem_data.stem_gridlist_b.tap))
                    obj.tapering_b=1;
                else
                    obj.tapering_b=0;
                end
            end
            if not(isempty(stem_data.stem_gridlist_p.tap))
                obj.tapering_p=1;
            else
                obj.tapering_p=0;
            end
            if not(isempty(stem_data.stem_gridlist_b))
                obj.tapering=obj.tapering_p|obj.tapering_b;
            else
                obj.tapering=obj.tapering_p;
            end
            
            if stem_par.theta_clustering>0
                stem_data.update_distance('both',1);
            end
        end
        
        function [aj_bp,aj_p] = get_aj(obj)
            %DESCRIPTION: provides the vector aj_bp and aj_p used in the EM estimation
            %
            %INPUT
            %obj   - [stem_model object] (1x1)
            %
            %OUTPUT
            %aj_bp - [double]            (Nx1) is the diagonal of the NxN diagonal matrix a*J_bp
            %aj_p  - [double]            (Nx1) is the diagonal of the NxN diagonal matrix a*J_p. 
            
            %NOTE
            %The elements of aj_p from Np+1 to N are all zeros. This allows the
            %use of stem_misc.D_apply both for the pixel data and the point
            %level data avoiding the use of J_bp and J_p
            if not(isempty(obj.stem_data.stem_varset_b))
                aj_bp=zeros(obj.stem_data.N,1);
                blocks=[0 cumsum(obj.stem_data.dim)];
                for i=1:obj.stem_data.nvar
                    aj_bp(blocks(i)+1:blocks(i+1))=obj.stem_par.alpha_bp(i);
                end
            else
                aj_bp=[];
            end
            if obj.stem_par.k>0
                aj_p=zeros(obj.stem_data.N,obj.stem_par.k);
                blocks=[0 cumsum(obj.stem_data.stem_varset_p.dim)];
                for k=1:obj.stem_par.k
                    for i=1:obj.stem_data.stem_varset_p.nvar
                        aj_p(blocks(i)+1:blocks(i+1),k)=obj.stem_par.alpha_p(i,k);    
                    end
                end
            else
                aj_p=[];
            end
        end
        
        function [aj_bp_b,j_b] = get_jbp(obj,r)
            %DESCRIPTION: provides the vectors aj_bp_b and j_b used in the EM estimation
            %
            %INPUT
            %obj     - [stem_model object] (1x1)
            %r       - [integer]           (1x1) is the index between 1 and q
            %
            %OUTPUT
            %aj_bp_b - [double]            (Nx1) is the vector with elements equal to alpha_bp(r) only for the sites of the r-th variable
            %j_b     - [double]            (Nx1) is the vector with elements equal to 1 only for the sites of the r-th variable
            
            aj_bp_b=zeros(obj.stem_data.N,1);
            j_b=zeros(obj.stem_data.N,1);
            blocks=[0 cumsum(obj.stem_data.dim)];
            for i=1:obj.stem_data.nvar
                if i==r
                    j_b(blocks(i)+1:blocks(i+1))=1;
                    aj_bp_b(blocks(i)+1:blocks(i+1))=obj.stem_par.alpha_bp(i);
                end
            end
        end
        
        function [aj_p_bs,j_b] = get_jg(obj,r,s)
            %DESCRIPTION: provides the vectors aj_p_bs and j_b used in the EM estimation
            %
            %INPUT
            %obj     - [stem_model object] (1x1)
            %r       - [integer]           (1x1) is the index between 1 and q
            %s       - [integer]           (1x1) is the index between 1 and q
            %
            %OUTPUT
            %aj_p_bs - [double]            (Nx1) is the vector with elements equal to alpha_p(i,s) only for the sites of the r-th variable
            %j_b     - [double]            (Nx1) is the vector with elements equal to 1 only for the sites of the r-th variable
            aj_p_bs=zeros(obj.stem_data.N,1);
            j_b=zeros(obj.stem_data.N,1);
            blocks=[0 cumsum(obj.stem_data.dim)];
            for i=1:obj.stem_data.stem_varset_p.nvar
                if i==r
                    j_b(blocks(i)+1:blocks(i+1))=1;
                    aj_p_bs(blocks(i)+1:blocks(i+1))=obj.stem_par.alpha_p(i,s); 
                end
            end
        end        
        
        function [sigma_eps,sigma_W_b,sigma_W_p,sigma_geo,sigma_Z,aj_bp,aj_p,M] = get_sigma(obj,sigma_W_b)
            %DESCRIPTION: provides the variance-covariance matrices and some vectors that are used in the EM algorithm
            %
            %INPUT
            %
            %obj         - [stem_model object] (1x1)
            %<sigma_W_b> - [double]            (NbxNb) (default: []) variance-covariance matrix of W_b. It is provided as input argument during kriging when sigma_W_b does not change across blocks  
            %
            %OUTPUT
            %
            %sigma_eps   - [double]            (NxN) variance-covariance matrix of epsilon
            %sigma_W_b   - [double]            (N_bxN_b) variance-covariance matrix of W_b
            %sigma_W_p   - [double]            {K}(N_pxN_p) variance-covariance matrices of the K W_p_i
            %sigma_geo   - [double]            (NxN) variance-covariance matrix of the sum of all the geostatistical components (Z excluded and epsilon included)
            %sigma_Z     - [double]            (pxp) variance-covariance of Z
            %aj_bp       - [double]            (Nx1) see the details of the method get_aj;
            %aj_p        - [double]            (Nx1) see the details of the method get_aj;
            %M           - [integer >0]        (N_px1) see the details of the method update_M of the class stem_data
            %
            %NOTE
            %sigma_geo is provided only if it is time-invariant otherwise it is evaluated at each step of the EM algorithm
            disp('    Marginal variance-covariance matrices evaluation started...');
            ct1=clock;
            if nargin==1
                sigma_W_b=[];
            end
            
            nvar=obj.stem_data.nvar;
            N=obj.stem_data.N;
            M=obj.stem_data.M;
            dim=obj.stem_data.dim;
            
            %sigma_eps
            d=[];
            for i=1:nvar
                d=[d;repmat(obj.stem_par.sigma_eps(i,i),dim(i),1)];
            end
            
            I=1:length(d);
            sigma_eps=sparse(I,I,d);
            
            %sigma_W_b
            if not(isempty(obj.stem_data.stem_varset_b))
                if not(isempty(obj.stem_data.stem_varset_b.X_bp))
                    if (nargin==1)||isempty(sigma_W_b)
                        if not(obj.tapering_b)
                            sigma_W_b=zeros(obj.stem_data.stem_varset_b.N);
                        end
                        blocks=[0 cumsum(obj.stem_data.stem_varset_b.dim)];
                        if obj.stem_par.pixel_correlated
                            if obj.tapering_b
                                I=zeros(nnz(obj.stem_data.DistMat_b),1);
                                J=zeros(nnz(obj.stem_data.DistMat_b),1);
                                elements=zeros(nnz(obj.stem_data.DistMat_b),1);
                            end
                            idx=0;
                            for j=1:obj.stem_data.stem_varset_b.nvar
                                for i=j:obj.stem_data.stem_varset_b.nvar
                                    idx_r=blocks(i)+1:blocks(i+1);
                                    idx_c=blocks(j)+1:blocks(j+1);
                                    if not(obj.tapering_b)
                                        sigma_W_b(idx_r,idx_c)=obj.stem_par.v_b(i,j)*stem_misc.correlation_function(...
                                            obj.stem_par.theta_b,obj.stem_data.DistMat_b(idx_r,idx_c),obj.stem_par.correlation_type);
                                        if not(i==j)
                                            sigma_W_b(idx_c,idx_r)=sigma_W_b(idx_r,idx_c)';
                                        end
                                    else
                                        corr_result=stem_misc.correlation_function(obj.stem_par.theta_b,obj.stem_data.DistMat_b(idx_r,idx_c),obj.stem_par.correlation_type);
                                        weights=stem_misc.wendland(obj.stem_data.DistMat_b(idx_r,idx_c),obj.stem_data.stem_gridlist_b.tap);
                                        corr_result.correlation=obj.stem_par.v_b(i,j)*corr_result.correlation.*weights;
                                        siz=length(corr_result.I);
                                        I(idx+1:idx+siz)=corr_result.I+blocks(i);
                                        J(idx+1:idx+siz)=corr_result.J+blocks(j);
                                        elements(idx+1:idx+siz)=corr_result.correlation;
                                        idx=idx+siz;
                                        if not(i==j)
                                            I(idx+1:idx+siz)=corr_result.J+blocks(j);
                                            J(idx+1:idx+siz)=corr_result.I+blocks(i);
                                            elements(idx+1:idx+siz)=corr_result.correlation;
                                            idx=idx+siz;
                                        end
                                    end
                                end
                            end
                            if obj.tapering_b
                                sigma_W_b=sparse(I,J,elements);
                            end
                        else
                            if obj.tapering_b
                                idx=0;
                                nonzeros=0;
                                for i=1:obj.stem_data.stem_varset_b.nvar
                                    nonzeros=nonzeros+nnz(obj.stem_data.DistMat_b(blocks(i)+1:blocks(i+1),blocks(i)+1:blocks(i+1)));
                                end
                                I=zeros(nonzeros,1);
                                elements=zeros(nonzeros,1);
                            end
                            for i=1:obj.stem_data.stem_varset_b.nvar
                                idx_rc=blocks(i)+1:blocks(i+1);
                                if not(obj.tapering_b)
                                    sigma_W_b(idx_rc,idx_rc)=stem_misc.correlation_function(obj.stem_par.theta_b(:,i),obj.stem_data.DistMat_b(idx_rc,idx_rc),obj.stem_par.correlation_type);
                                else
                                    corr_result=stem_misc.correlation_function(obj.stem_par.theta_b(:,i),obj.stem_data.DistMat_b(idx_rc,idx_rc),obj.stem_par.correlation_type);
                                    weights=stem_misc.wendland(obj.stem_data.DistMat_b(idx_rc,idx_rc),obj.stem_data.stem_gridlist_b.tap);
                                    corr_result.correlation=obj.stem_par.v_b(i,i)*corr_result.correlation.*weights;
                                    siz=length(corr_result.I);
                                    I(idx+1:idx+siz)=corr_result.I+blocks(i);
                                    J(idx+1:idx+siz)=corr_result.J+blocks(i);
                                    elements(idx+1:idx+siz)=corr_result.correlation;
                                    idx=idx+siz;
                                end
                            end
                            if obj.tapering_b
                                sigma_W_b=sparse(I,J,elements);
                            end
                        end
                    end
                end
            end
            clear I
            clear J
            clear elements
            clear weights
            clear corr_result
            
           %sigma_W_p
           if obj.stem_par.k>0
               blocks=[0 cumsum(obj.stem_data.stem_varset_p.dim)];
               for k=1:obj.stem_par.k
                   if obj.tapering_p
                       I=zeros(nnz(obj.stem_data.DistMat_p),1);
                       J=zeros(nnz(obj.stem_data.DistMat_p),1);
                       elements=zeros(nnz(obj.stem_data.DistMat_p),1);
                   else
                       sigma_W_p{k}=obj.stem_data.DistMat_p;   
                   end
                   idx=0;
                   for j=1:obj.stem_data.stem_varset_p.nvar
                       for i=j:obj.stem_data.stem_varset_p.nvar
                           idx_r=blocks(i)+1:blocks(i+1);
                           idx_c=blocks(j)+1:blocks(j+1);
                           if not(obj.tapering_p)
                               sigma_W_p{k}(idx_r,idx_c)=obj.stem_par.v_p(i,j,k)*stem_misc.correlation_function(...
                                   obj.stem_par.theta_p(:,k),obj.stem_data.DistMat_p(idx_r,idx_c),obj.stem_par.correlation_type);
                               if not(i==j)
                                   sigma_W_p{k}(idx_c,idx_r)=sigma_W_p{k}(idx_r,idx_c)';
                               end
                           else
                               corr_result=stem_misc.correlation_function(obj.stem_par.theta_p(:,k),obj.stem_data.DistMat_p(idx_r,idx_c),obj.stem_par.correlation_type);
                               weights=stem_misc.wendland(obj.stem_data.DistMat_p(idx_r,idx_c),obj.stem_data.stem_gridlist_p.tap);
                               corr_result.correlation=obj.stem_par.v_p(i,j,k)*corr_result.correlation.*weights;
                               siz=length(corr_result.I);
                               I(idx+1:idx+siz)=corr_result.I+blocks(i);
                               J(idx+1:idx+siz)=corr_result.J+blocks(j);
                               elements(idx+1:idx+siz)=corr_result.correlation;
                               idx=idx+siz;
                               if not(i==j)
                                   I(idx+1:idx+siz)=corr_result.J+blocks(j);
                                   J(idx+1:idx+siz)=corr_result.I+blocks(i);
                                   elements(idx+1:idx+siz)=corr_result.correlation;
                                   idx=idx+siz;
                               end
                           end
                       end
                   end
                   if obj.tapering_p
                       sigma_W_p{k}=sparse(I,J,elements);
                   end
               end
           else
               sigma_W_p=[];
           end
           clear I
           clear J
           clear elements
           clear weights
           clear corr_result
           
           %sigma_geo
           [aj_bp,aj_p]=obj.get_aj;
           if not(obj.stem_data.X_tv)
               %time invariant case
               if not(isempty(obj.stem_data.stem_varset_b))
                   sigma_geo=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),obj.stem_data.X_bp(:,1,1),'b'),aj_bp,'b');
               end
               if obj.stem_par.k>0
                   if isempty(obj.stem_data.stem_varset_b)
                       %se manca il remoto allora sigma_geo non � stata
                       %ancora allocata
                       if obj.tapering
                           sigma_geo=spalloc(size(sigma_W_p{1},1),size(sigma_W_p{1},1),nnz(sigma_W_p{1}));
                       else
                           sigma_geo=zeros(N);
                       end
                   end
                   for k=1:obj.stem_par.k
                       sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},obj.stem_data.X_p(:,1,1,k),'b'),aj_p(:,k),'b');
                   end
               end
               if (obj.stem_par.k==0)&&(isempty(obj.stem_data.stem_varset_b))
                   sigma_geo=sigma_eps;
               else
                   sigma_geo=sigma_geo+sigma_eps;
               end
           else
               sigma_geo=[]; %it means sigma_geo is time-variant and has to be evaluated at each time step
           end
           
           if obj.stem_par.p>0
               sigma_Z=(eye(obj.stem_par.p^2)-kron(obj.stem_par.G,obj.stem_par.G))\obj.stem_par.sigma_eta(:);%variance of a VAR (Lutkepohl pag.27)
               sigma_Z=reshape(sigma_Z,obj.stem_par.p,obj.stem_par.p);
           else
               sigma_Z=[];
           end   
           ct2=clock;
           disp(['    Marginal variance-covariance matrices evaluation ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
        end
        
        function simulate(obj,nan_rate,nan_pattern_par)
            %DESCRIPTION: is the front-end method to simulate data from this stem_model object
            %
            %INPUT
            %
            %obj               - [stem_model object] (1x1)
            %<nan_rate>        - [double [0,1]]      (Sx1) (default: []) missing rate of the simulated data
            %<nan_pattern_par> - [double >0]         (Sx1) (default: []) (UoM: km) parameter of the exponential spatial correlation function used to define the spatial pattern of the missing data
            %
            %OUTPUT
            %
            %none: an object of class stem_sim is created. The method stem_sim.simulate actually simulates the data and it fills the stem_data object
            st_sim=stem_sim(obj);
            if nargin<2
                nan_rate=[];
                nan_pattern_par=[];
            end
            if nargin==2
                error('nan_pattern_par must also be provided');
            end
            st_sim.simulate(nan_rate,nan_pattern_par);
            obj.stem_data.simulated=1;
        end        
        
        function EM_estimate(obj,stem_EM_options)
            %DESCRIPTION: front-end method for EM estimation of this model
            %
            %INPUT
            %
            %obj                - [stem_model object]      (1x1)
            %stem_EM_options    - [stem_EM_options object] (1x)
            %
            %OUTPUT
            %
            %none: the stem_EM object is created. The method stem_EM.estimate is used to estimate the model and it updates the stem_par object

            standardized=1;
            if not(obj.stem_data.stem_varset_p.standardized)
               standardized=0;
            end
            if not(isempty(obj.stem_data.stem_varset_b))
                if not(obj.stem_data.stem_varset_b.standardized)
                    standardized=0;
                end
            end
            if not(standardized)
                error('Data must be standardized in order to avoid numerical stability problems. Use the method standardize of class stem_data.');
            end
            if not(isempty(obj.stem_data.stem_crossval))
                disp('Data modification for cross-validation started...');
                idx_var=obj.stem_data.stem_varset_p.get_Y_index(obj.stem_data.stem_crossval.variable_name);
                if isempty(idx_var)
                    error('Variable not found. Has it been deleted?');
                end
                
                %recover the indices of the cross-validation sites 
                indices=obj.stem_data.stem_crossval.indices;
                Y={obj.stem_data.stem_varset_p.Y{idx_var}(indices,:)};
                Y_name={obj.stem_data.stem_varset_p.Y_name{idx_var}};
                if not(isempty(obj.stem_data.stem_varset_p.X_bp))
                    X_bp={obj.stem_data.stem_varset_p.X_bp{idx_var}(indices,:,:)};
                    X_bp_name={obj.stem_data.stem_varset_p.X_bp{idx_var}};
                else
                    X_bp={};
                    X_bp_name={};
                end
                if not(isempty(obj.stem_data.stem_varset_p.X_beta))
                    X_beta={obj.stem_data.stem_varset_p.X_beta{idx_var}(indices,:,:)};
                    X_beta_name={obj.stem_data.stem_varset_p.X_beta_name{idx_var}};
                else
                    X_beta={};
                    X_beta_name={};
                end 
                if not(isempty(obj.stem_data.stem_varset_p.X_z))
                    X_z={obj.stem_data.stem_varset_p.X_z{idx_var}(indices,:,:)};
                    X_z_name={obj.stem_data.stem_varset_p.X_z_name{idx_var}};
                else
                    X_z={};
                    X_z_name={};
                end          
                if not(isempty(obj.stem_data.stem_varset_p.X_p))
                    X_p={obj.stem_data.stem_varset_p.X_p{idx_var}(indices,:,:,:)};
                    X_p_name={obj.stem_data.stem_varset_p.X_p_name{idx_var}};
                else
                    X_p={};
                    X_p_name={};
                end      
                
                %set the cross_mindistance vector
                if not(isempty(obj.stem_data.DistMat_p))
                    dim=obj.stem_data.dim;
                    blocks=[0 cumsum(dim)];
                    temp_dist=obj.stem_data.DistMat_p(blocks(idx_var)+1:blocks(idx_var+1),blocks(idx_var)+1:blocks(idx_var+1));
                    temp_dist=temp_dist(indices,:);
                    temp_dist(:,indices)=[];
                    obj.stem_data.stem_crossval.min_distance=min(temp_dist');
                    clear temp_dist
                end
                
                obj.stem_data.stem_crossval.stem_varset=stem_varset(Y,Y_name,X_bp,X_bp_name,X_beta,X_beta_name,X_z,X_z_name,X_p,X_p_name);
                obj.stem_data.stem_crossval.stem_gridlist=stem_gridlist();
                coordinate=obj.stem_data.stem_gridlist_p.grid{idx_var}.coordinate;
                st_grid=stem_grid(coordinate(indices,:),'deg','sparse','point');
                obj.stem_data.stem_crossval.stem_gridlist.add(st_grid);
                %remove the cross-validation data from the estimation dataset
                obj.stem_data.site_crop(obj.stem_data.stem_crossval.type,obj.stem_data.stem_crossval.variable_name,indices);
                obj.cross_validation=1;
                disp('Data modification ended.');
            else
                obj.cross_validation=0;
            end
            
            st_EM=stem_EM(obj,stem_EM_options);
            %set the current parameter value with the estimated initial value
            obj.stem_par=obj.stem_par_initial;
            if isempty(stem_EM_options.path_distributed_computing)
                obj.stem_EM_result=st_EM.estimate();
            else
                obj.stem_EM_result=st_EM.estimate_parallel(stem_EM_options.path_distributed_computing);
            end
            obj.estimated=1;
            if obj.cross_validation
                st_krig=stem_krig(obj);
                block_size=1000;
                back_transform=0;
                no_varcov=1;
                crossval=1;
                obj.stem_data.stem_crossval.stem_krig_result=st_krig.kriging(obj.stem_data.stem_crossval.variable_name,obj.stem_data.stem_crossval.stem_gridlist.grid{1},block_size,[],[],back_transform,no_varcov,crossval);
                obj.stem_data.stem_crossval.res=obj.stem_data.stem_crossval.stem_varset.Y{idx_var}-obj.stem_data.stem_crossval.stem_krig_result.y_hat;
                obj.stem_data.stem_crossval.mse=nanvar(obj.stem_data.stem_crossval.res');
                obj.stem_data.stem_crossval.relative_mse=obj.stem_data.stem_crossval.mse./nanvar(obj.stem_data.stem_crossval.stem_varset.Y{idx_var}');

                s=obj.stem_data.stem_varset_p.Y_stds{idx_var};
                m=obj.stem_data.stem_varset_p.Y_means{idx_var};
                if (obj.stem_data.stem_varset_p.standardized)&&not(obj.stem_data.stem_varset_p.log_transformed)
                    y_hat_back=obj.stem_data.stem_crossval.stem_krig_result.y_hat*s+m;
                    y=obj.stem_data.stem_crossval.stem_varset.Y{idx_var}*s+m;
                end
                
                if (obj.stem_data.stem_varset_p.standardized)&&(obj.stem_data.stem_varset_p.log_transformed)
                    y_hat_back=obj.stem_data.stem_crossval.stem_krig_result.y_hat;
                    st=nanstd(obj.stem_data.stem_varset_p.Y{idx_var});
                    st=repmat(st,[size(y_hat_back,1),1]);
                    st=st.^2*s;
                    y_hat_back=exp(y_hat_back*s+m+st/2);
                    y=exp(obj.stem_data.stem_crossval.stem_varset.Y{idx_var}*s+m);
                end
                obj.stem_data.stem_crossval.res_back=y-y_hat_back;
                obj.stem_data.stem_crossval.y_back=y;
                obj.stem_data.stem_crossval.y_hat_back=y_hat_back;
            end
        end
          
        function print(obj)
            %DESCRIPTION: print the information on the estimation result
            %
            %INPUT
            %obj  - [stem_model object]   (1x1) the stem_data object
            %
            %OUTPUT
            %
            %none: the information is printed in the command window          
            
            stem_misc.disp_star('Model estimation results')   
            if obj.estimated
                if obj.tapering
                    disp('* Tapering is enabled.');    
                    if not(isempty(obj.tapering_p))
                        disp(['  Point data tapering: ',num2str(obj.stem_data.stem_gridlist_p.tap),' km']);
                    else
                        disp('   Tapering is NOT enabled on point data');
                    end
                    if not(isempty(obj.tapering_b))
                        disp(['  Pixel data tapering: ',num2str(obj.stem_data.stem_gridlist_b.tap),' km']);
                    else
                        disp('  Tapering is NOT enabled on pixel data');
                    end
                else
                    disp('* Tapering is not enabled');
                end
                disp(' ');
                if not(isempty(obj.stem_EM_result.logL))
                    disp(['* Observed data log-likelihood: ',num2str(obj.stem_EM_result.logL,'%05.3f')]);
                else
                    disp('* Observed data log-likelihood: not computed. Use the method set_logL of the class stem_model.');
                end
                disp(' ');
                counter=1;
                if not(isempty(obj.stem_par.beta))
                    for i=1:obj.stem_data.stem_varset_p.nvar
                        disp(['* Beta coefficients related to the point variable ',obj.stem_data.stem_varset_p.Y_name{i}]);
                        output=[];
                        output{1,1}='Loading coefficient';
                        output{1,2}='Value';
                        output{1,3}='Std';
                        for j=1:length(obj.stem_data.stem_varset_p.X_beta_name{i})
                            output{j+1,1}=obj.stem_data.stem_varset_p.X_beta_name{i}{j};
                            output{j+1,2}=num2str(obj.stem_par.beta(counter),'%+05.3f');
                            if not(isempty(obj.stem_EM_result.varcov))
                                output{j+1,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                            else
                                output{j+1,3}='Not computed';
                            end
                            counter=counter+1;
                        end
                        disp(output);
                    end
                    if not(isempty(obj.stem_data.stem_varset_b))
                        if not(isempty(obj.stem_data.stem_varset_b.X_beta))
                            for i=1:obj.stem_data.stem_varset_b.nvar
                                disp(['* Beta coefficients related to the pixel variable ',obj.stem_data.stem_varset_b.Y_name{i}]);
                                output=[];
                                output{1,1}='Loading coefficient';
                                output{1,2}='Value';
                                output{1,3}='Std';
                                for j=1:length(obj.stem_data.stem_varset_b.X_beta_name{i})
                                    output{j+1,1}=obj.stem_data.stem_varset_b.X_beta_name{i}{j};
                                    output{j+1,2}=num2str(obj.stem_par.beta(counter),'%+05.3f');
                                    if not(isempty(obj.stem_EM_result.varcov))
                                        output{j+1,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                                    else
                                        output{j+1,3}='Not computed';
                                    end
                                    counter=counter+1;
                                end
                                disp(output);
                            end
                        end
                    end
                end
                output=[];
                disp('* Sigma_eps diagonal elements (Variance)')
                output{1,1}='Variable';
                output{1,2}='Value';
                output{1,3}='Std';
                for i=1:obj.stem_data.stem_varset_p.nvar
                    output{i+1,1}=obj.stem_data.stem_varset_p.Y_name{i};
                    output{i+1,2}=num2str(obj.stem_par.sigma_eps(i,i),'%05.3f');
                    if not(isempty(obj.stem_EM_result.varcov))
                        output{i+1,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                    else
                        output{i+1,3}='Not computed';
                    end
                    counter=counter+1;
                end
                if not(isempty(obj.stem_data.stem_varset_b))
                    delta=obj.stem_data.stem_varset_p.nvar;
                    for i=1:obj.stem_data.stem_varset_b.nvar
                        output{i+1+delta,1}=obj.stem_data.stem_varset_b.Y_name{i};
                        output{i+1+delta,2}=num2str(obj.stem_par.sigma_eps(i+delta,i+delta),'%05.3f');
                        if not(isempty(obj.stem_EM_result.varcov))
                            output{i+1+delta,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                        else
                            output{i+1+delta,3}='Not computed';
                        end
                        counter=counter+1;
                    end
                end
                disp(output);
                output=[];
                if not(isempty(obj.stem_data.stem_varset_b))
                    disp('* alpha_bp elements')
                    output{1,1}='Variable';
                    output{1,2}='Value';
                    output{1,3}='Std';
                    for i=1:obj.stem_data.stem_varset_p.nvar
                        output{i+1,1}=obj.stem_data.stem_varset_p.Y_name{i};
                        output{i+1,2}=num2str(obj.stem_par.alpha_bp(i),'%+05.3f');
                        if not(isempty(obj.stem_EM_result.varcov))
                            output{i+1,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                        else
                            output{i+1,3}='Not computed';
                        end
                        counter=counter+1;
                    end
                    delta=obj.stem_data.stem_varset_p.nvar;
                    for i=1:obj.stem_data.stem_varset_b.nvar
                        output{i+1+delta,1}=obj.stem_data.stem_varset_b.Y_name{i};
                        output{i+1+delta,2}=num2str(obj.stem_par.alpha_bp(i+delta),'%+05.3f');
                        if not(isempty(obj.stem_EM_result.varcov))
                            output{i+1+delta,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                        else
                            output{i+1+delta,3}='Not computed';
                        end
                        counter=counter+1;
                    end
                    disp(output);
                    output=[];
                    if obj.stem_par.pixel_correlated
                        disp('* Pixel data are cross-correlated.');
                        disp(' ');
                        output{1,2}='Value [km]';
                        output{1,3}='Std [km]';
                        output{2,1}='Theta_b';
                        output{2,2}=num2str(obj.stem_par.theta_b(1),'%05.3f');
                        if not(isempty(obj.stem_EM_result.varcov))
                            output{2,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                        else
                            output{2,3}='Not computed';
                        end
                        counter=counter+1;
                        disp(output);
                        output=[];
                        disp('* v_b matrix:');
                        for i=1:obj.stem_data.stem_varset_b.nvar
                            output{1,i+1}=obj.stem_data.stem_varset_b.Y_name{i};
                            output{i+1,1}=obj.stem_data.stem_varset_b.Y_name{i};
                            output{i+1,i+1}=num2str(1,'%+5.2f');
                        end
                        for i=1:obj.stem_data.stem_varset_b.nvar
                            for j=i+1:obj.stem_data.stem_varset_b.nvar
                                if not(isempty(obj.stem_EM_result.varcov))
                                    output{i+1,j+1}=[num2str(obj.stem_par.v_b(i,j),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.2f'),')'];
                                    output{j+1,i+1}=output{i+1,j+1};
                                else
                                    output{i+1,j+1}=num2str(obj.stem_par.v_b(i,j),'%+05.2f');
                                    output{j+1,i+1}=output{i+1,j+1};
                                end
                                counter=counter+1;
                            end
                        end
                        disp(output);
                        output=[];
                    else
                        disp('* Pixel data are NOT cross-correlated.');
                        disp(' ');
                        disp('* Theta_b elements:');
                        output{1,1}='Variable';
                        output{1,2}='Value [km]';
                        output{1,3}='Std [km]';
                        for i=1:obj.stem_data.stem_varset_b.nvar
                            output{i+1,1}=obj.stem_data.stem_varset_b.Y_name{i};
                            output{i+1,2}=num2str(obj.stem_par.theta_b(i),'%05.3f');
                            if not(isempty(obj.stem_EM_result.varcov))
                                output{i+1,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                            else
                                output{i+1,3}='Not computed';
                            end
                            counter=counter+1;
                        end
                        disp(output);
                        output=[];
                    end
                end
                if obj.stem_par.k>0
                    disp(['* ',num2str(obj.stem_par.k),' fine-scale coregionalization components w_p']);
                    disp(' ');
                    disp(['* alpha_p elements:'])
                    for i=1:obj.stem_data.stem_varset_p.nvar
                        output{i*2,1}=obj.stem_data.stem_varset_p.Y_name{i};
                        for k=1:obj.stem_par.k
                            output{i*2-1,k+1}=obj.stem_data.stem_varset_p.X_p_name{i}{k};
                            if not(isempty(obj.stem_EM_result.varcov))
                                output{i*2,k+1}=[num2str(obj.stem_par.alpha_p(i,k),'%+05.3f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f'),')'];
                            else
                                output{i*2,k+1}=num2str(obj.stem_par.alpha_p(i,k),'%+05.3f');
                            end
                            counter=counter+1;
                        end
                    end
                    disp(output);
                    output=[];
                    disp(['* theta_p elements:']);
                    output{1,1}='Coreg. component';
                    output{1,2}='Value [km]';
                    output{1,3}='Std [km]';
                    for k=1:obj.stem_par.k
                        if k==1
                            postfix='st';
                        end
                        if k==2
                            postfix='nd';
                        end
                        if k==3
                            postfix='rd';
                        end
                        if k>3
                            postfix='th';
                        end
                        output{k+1,1}=[num2str(k),postfix];
                        output{k+1,2}=num2str(obj.stem_par.theta_p(:,k),'%06.2f');
                        if not(isempty(obj.stem_EM_result.varcov))
                            output{k+1,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.2f');
                        else
                            output{k+1,3}='Not computed';
                        end
                        counter=counter+1;
                    end
                    disp(output);
                    output=[];
                    for k=1:obj.stem_par.k
                        if k==1
                            postfix='st';
                        end
                        if k==2
                            postfix='nd';
                        end
                        if k==3
                            postfix='rd';
                        end
                        if k>3
                            postfix='th';
                        end                        
                        disp(['* v_p matrix for the ',num2str(k),postfix,' coreg. component:']);
                        for i=1:obj.stem_data.stem_varset_p.nvar
                            output{1,i+1}=obj.stem_data.stem_varset_p.Y_name{i};
                            output{i+1,1}=obj.stem_data.stem_varset_p.Y_name{i};
                            output{i+1,i+1}=num2str(1,'%+5.2f');
                        end
                        for i=1:obj.stem_data.stem_varset_p.nvar
                            for j=i+1:obj.stem_data.stem_varset_p.nvar
                                if not(isempty(obj.stem_EM_result.varcov))
                                    output{i+1,j+1}=[num2str(obj.stem_par.v_p(i,j,k),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%03.2f'),')'];
                                    output{j+1,i+1}=output{i+1,j+1};
                                else
                                    output{i+1,j+1}=num2str(obj.stem_par.v_p(i,j,k),'%+05.2f');
                                    output{j+1,i+1}=output{i+1,j+1};
                                end
                                counter=counter+1;
                            end
                        end                        
                        disp(output);
                        output=[];
                    end
                end
                if obj.stem_par.p>0
                    disp('* Transition matrix G:');
                    c=1;
                    for i=1:obj.stem_data.stem_varset_p.nvar
                        for j=1:size(obj.stem_data.stem_varset_p.X_z{i},2)
                            output{1,c+1}=[obj.stem_data.stem_varset_p.Y_name{i},' - ',obj.stem_data.stem_varset_p.X_z_name{i}{j}];
                            output{c+1,1}=output{1,c+1};
                            c=c+1;
                        end
                    end
                    if not(isempty(obj.stem_data.stem_varset_b))
                        if not(isempty(obj.stem_data.stem_varset_b.X_z))
                            for i=1:obj.stem_data.stem_varset_b.nvar
                                for j=1:size(obj.stem_data.stem_varset_b.X_z{i},2)
                                    output{1,c+1}=[obj.stem_data.stem_varset_b.Y_name{i},' - ',obj.stem_data.stem_varset_b.X_z_name{i}{j}];
                                    output{c+1,1}=output{1,c+1};
                                    c=c+1;
                                end
                            end
                        end
                    end
                    if obj.stem_par.time_diagonal
                        for j=1:obj.stem_par.p
                            for i=1:obj.stem_par.p
                                if i==j
                                    if not(isempty(obj.stem_EM_result.varcov))
                                        output{i+1,i+1}=[num2str(obj.stem_par.G(i,i),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%03.2f'),')'];
                                    else
                                        output{i+1,i+1}=num2str(obj.stem_par.G(i,i),'%+05.2f');
                                    end
                                    counter=counter+1;
                                else
                                    output{i+1,j+1}='0';
                                end
                            end
                        end
                    else
                        for j=1:obj.stem_par.p
                            for i=1:obj.stem_par.p
                                if not(isempty(obj.stem_EM_result.varcov))
                                    output{i+1,j+1}=[num2str(obj.stem_par.G(i,j),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%03.2f'),')'];
                                else
                                    output{i+1,j+1}=num2str(obj.stem_par.G(i,j),'%+05.2f');
                                end
                                counter=counter+1;
                            end
                        end
                    end
                    disp(output);
                    output=[];
                    disp('* Sigma_eta matrix:');
                    c=1;
                    for i=1:obj.stem_data.stem_varset_p.nvar
                        for j=1:size(obj.stem_data.stem_varset_p.X_z{i},2)
                            output{1,c+1}=[obj.stem_data.stem_varset_p.Y_name{i},' - ',obj.stem_data.stem_varset_p.X_z_name{i}{j}];
                            output{c+1,1}=output{1,c+1};
                            c=c+1;
                        end
                    end
                    if not(isempty(obj.stem_data.stem_varset_b))
                        if not(isempty(obj.stem_data.stem_varset_b.X_z))
                            for i=1:obj.stem_data.stem_varset_b.nvar
                                for j=1:size(obj.stem_data.stem_varset_b.X_z{i},2)
                                    output{1,c+1}=[obj.stem_data.stem_varset_b.Y_name{i},' - ',obj.stem_data.stem_varset_b.X_z_name{i}{j}];
                                    output{c+1,1}=output{1,c+1};
                                    c=c+1;
                                end
                            end
                        end
                    end
                    if obj.stem_par.time_diagonal
                        for j=1:obj.stem_par.p
                            for i=1:obj.stem_par.p
                                if i==j
                                    if not(isempty(obj.stem_EM_result.varcov))
                                        output{i+1,i+1}=[num2str(obj.stem_par.sigma_eta(i,i),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%03.2f'),')'];
                                    else
                                        output{i+1,i+1}=num2str(obj.stem_par.sigma_eta(i,i),'%+05.2f');
                                    end
                                    counter=counter+1;
                                else
                                    output{i+1,j+1}='0';
                                end
                            end
                        end
                    else
                        for j=1:obj.stem_par.p
                            for i=j:obj.stem_par.p
                                if not(isempty(obj.stem_EM_result.varcov))
                                    output{i+1,j+1}=[num2str(obj.stem_par.sigma_eta(i,j),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%03.2f'),')'];
                                else
                                    output{i+1,j+1}=num2str(obj.stem_par.sigma_eta(i,j),'%+05.2f');
                                end
                                output{j+1,i+1}=output{i+1,j+1};
                                counter=counter+1;
                            end
                        end
                    end
                    disp(output);
                  
                end
            else
                disp('The model has not been estimated yet. Use the method print of the class stem_data to print data information.');
            end            
        end
        
        function set_system_size(obj,dim)
            %DESCRIPTION: evaluate the minimum N after which it is faster to evaluate only the diagonal of a matrix product instead of the full matrix
            %
            %INPUT
            %
            %obj - [stem_model object] (1x1)
            %
            %OUTPUT
            %
            %none: the system_size property is updated
            if nargin<2
                dim=100;
                t1=0;
                t2=100;
                while t1<t2
                    m1=randn(dim);
                    m2=randn(dim);
                    tic
                    m=diag(m1*m2);
                    t1=toc;
                    tic
                    for i=1:size(m1)
                        h=m1(i,:)*m2(:,i);
                    end
                    t2=toc;
                    dim=dim+50;
                end
                obj.system_size=dim;
            else
                obj.system_size=dim;
            end
        end
        
        function set_logL(obj)
            %DESCRIPTION: evaluate the observed data log-likelihood
            %
            %INPUT
            %
            %obj - [stem_model object] (1x1)
            %
            %OUTPUT
            %
            %none: the logL property of the object stem_EM_result is updated
            disp('Log-Likelihood computation...');
            st_kalman=stem_kalman(obj);
            [st_kalmansmoother_result,~,~,~,~,~,~,~,~] = st_kalman.smoother(1);
            obj.stem_EM_result.logL=st_kalmansmoother_result.logL;
            disp('Log-Likelihood computation ended.');
        end
        
        function set_varcov(obj)
            %DESCRIPTION: evaluate the variance-covariance matrix of the model parameters
            %
            %INPUT
            %
            %obj - [stem_model object] (1x1)
            %
            %OUTPUT
            %
            %none: the varcov property of the object stem_EM_result is updated    
            
            %parameter order: beta,sigma_eps,alpha_bp,theta_b,v_b,alpha_p,theta_p,v_p,G,sigma_eta
            
            if obj.estimated==0
                error('The model has not been estimated yet');
            end
            
            N=obj.N;
            Np=obj.Np;
            T=obj.T;
            dim=obj.dim;            
            
            data=obj.stem_data;
            par=obj.stem_par;
            q=par.q;
            p=par.p;
            k=par.k;

            %evaluate the number of parameters
            n_beta=0;
            n_bp_alpha=0;
            n_bp_theta=0;
            n_bp_v=0;
            n_p_alpha=0;
            n_p_theta=0;
            n_p_v=0;
            n_time_G=0;
            n_time_s2e=0;
                       
            if not(isempty(data.X_bp))
                n_bp_alpha=2*q;
                if par.pixel_correlated
                    n_bp_v=q*(q-1)/2; %v_bp extra-diagonal and univoc
                    n_bp_theta=1;
                else
                    %nothing to do for V as the V matrix is the identity
                    n_bp_theta=q;
                end
            end
            
            if not(isempty(data.X_beta))
                n_beta=length(par.beta);
            end
            
            if p>0
                if par.time_diagonal
                    n_time_G=p;
                    n_time_s2e=p;
                else
                    n_time_G=p^2;
                    n_time_s2e=(p+1)*p/2;
                end
                n_time=n_time_G+n_time_s2e;
            else
                n_time=0;
            end
            
            if not(isempty(data.X_p))
                n_p_alpha=q*k; %alpha_p
                n_p_v=q*(q-1)/2*k;
                n_p_theta=k;
            end
            n_eps=size(par.sigma_eps,1);
            
            n_psi=n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+n_p_v+n_time_G+n_time_s2e;
            
            if p>0
                %kalman filter
                st_kalman=stem_kalman(obj);
                compute_logL=0;
                enable_varcov_computation=1;
                [st_kalmanfilter_result,sigma_eps,sigma_W_b,sigma_W_p,~,aj_bp,aj_p,M,sigma_geo] = st_kalman.filter(compute_logL,enable_varcov_computation);
            else
                [sigma_eps,sigma_W_b,sigma_W_p,sigma_geo,~,aj_bp,aj_p,M] = obj.get_sigma();
                st_kalmanfilter_result=stem_kalmanfilter_result([],[],[],[],[],[],[]);
            end            
            J=st_kalmanfilter_result.J(:,:,2:end); %J for t=0 is deleted
            st_kalmanfilter_result.J=st_kalmanfilter_result.J(:,:,2:end);
            st_kalmanfilter_result.zk_f=st_kalmanfilter_result.zk_f(:,2:end);
            st_kalmanfilter_result.Pk_f=st_kalmanfilter_result.Pk_f(:,:,2:end);            
           
            disp('Derivatives allocation...');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  allocation time variant derivatives  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            d_P=zeros(p,p,n_psi); %the first n_beta derivatives are zero
            d_P_lag=zeros(p,p,n_psi);
            
            %d_J=zeros(p,N,n_psi); %the first n_beta derivatives are zero
            %d_J_lag=zeros(p,N,n_psi);
            for i=1:n_psi
                d_J_Lt{i}=sparse(p,N);
            end
            
            d_Z=zeros(p,n_psi); %the first n_beta derivatives are zero
            d_Z_lag=zeros(p,n_psi);
            
            d_e=zeros(N,n_psi);
            d_e_lag=zeros(N,n_psi);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  allocation of time invariant derivatives  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            d_beta=zeros(n_beta,1,n_psi);
            d_G=zeros(p,p,n_psi);
            d_s2e=zeros(p,p,n_psi);
            
            disp('Derivatives initialization...');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  initialization of time invariant derivatives  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            blocks=[0 cumsum(dim)];
            %beta
            if n_beta>0
                for i=1:n_beta
                    temp=zeros(n_beta,1);
                    temp(i)=1;
                    d_beta(:,1,i)=temp;
                end
            end
            %G
            if par.time_diagonal
                for i=1:p
                    temp=zeros(p);
                    temp(i,i)=1;
                    d_G(:,:,n_psi-n_time+i)=temp;
                end
            else
                for i=1:p^2
                    temp=zeros(p^2,1);
                    temp(i)=1;
                    d_G(:,:,n_psi-n_time+i)=reshape(temp,p,p);
                end
            end
            %sigma2_eta
            if par.time_diagonal
                for i=1:p
                    temp=zeros(p);
                    temp(i,i)=1;
                    d_s2e(:,:,n_psi-n_time_s2e+i)=temp;
                end
            else
                j=1;
                for i=1:(p*(p+1))/2
                    temp=zeros(p);
                    l=1;
                    for h=1:p
                        for k=h:p
                            if i==l
                                temp(h,k)=1;
                                temp(k,h)=1;
                            end
                            l=l+1;
                        end
                    end
                    d_s2e(:,:,n_psi-n_time_s2e+i)=temp;
                    j=j+1;
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %      compute d_Sgeo_prel     %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %beta
            for i=1:n_beta
                d_Sgeo_prel{i}=sparse(N,N);
            end
            %sigma_eps
            for i=n_beta+1:n_beta+n_eps
                Id=blocks(i-n_beta)+1:blocks(i-n_beta+1);
                d_Sgeo_prel{i}=sparse(Id,Id,ones(length(Id),1),N,N);
            end
            %alpha_bp
            if n_bp_alpha>0
                result=stem_misc.M_apply(sigma_W_b,M,'b');
                for j=1:n_bp_alpha
                    Id=[];
                    Jd=[];
                    elements=[];
                    for i=1:n_bp_alpha
                        temp=result(blocks(i)+1:blocks(i+1),blocks(j)+1:blocks(j+1));
                        if i==j
                            temp2=2*par.alpha_bp(j)*temp;
                            L=find(temp2);
                            [idx_I,idx_J]=ind2sub(size(temp2),L);
                            Id=[Id; idx_I+blocks(i)];
                            Jd=[Jd; idx_J+blocks(j)];
                            elements=[elements; temp2(L)];
                        else
                            temp2=par.alpha_bp(i)*temp;
                            L=find(temp2);
                            [idx_I,idx_J]=ind2sub(size(temp2),L);
                            Id=[Id; idx_I+blocks(i)];
                            Jd=[Jd; idx_J+blocks(j)];
                            elements=[elements; temp2(L)];
                            Id=[Id; idx_J+blocks(j)];
                            Jd=[Jd; idx_I+blocks(i)];
                            elements=[elements; temp2(L)];
                        end
                    end
                    zero_density=(1-length(Id)/(N^2))*100;
                    if obj.tapering||zero_density>60
                        d_Sgeo_prel{n_beta+n_eps+j}=sparse(Id,Jd,elements,N,N);
                    else
                        d_Sgeo_prel{n_beta+n_eps+j}=zeros(N);
                        temp=sub2ind([N,N],Id,Jd);
                        d_Sgeo_prel{n_beta+n_eps+j}(temp)=elements;
                    end
                end
            end
            %theta_bp
            
            if n_bp_theta>0
                d=stem_misc.M_apply(obj.stem_data.DistMat_b,M,'b');
                if not(par.pixel_correlated)&&(q>1)
                    for j=1:q
                        Id=[];
                        Jd=[];
                        elements=[];
                        for i=1:2 
                            idx=blocks(j+(i-1)*2)+1:blocks(j+(i-1)*2+1);
                            result1=d(idx,idx);
                            result2=result(idx,idx);
                            result2=stem_misc.D_apply(result2,aj_bp(idx),'b');
                            if strcmp(par.correlation_type,'exponential')
                                result3=result1/(par.theta_b(j)^2).*result2;
                            elseif strcmp(par.correleation_type,'matern32')
                                result3=-sqrt(3)*result1/(par.theta_b(j)^2).*result2./(1+sqrt(3).*result1./par.theta_b(j))+result2.*sqrt(3).*result1./(par.theta_b(j)^2);
                            else
                                result3=(-sqrt(5)*result1/(par.theta_b(j)^2)-10/3*result1.^2/par.theta_b(j)^3).*result2./(1+sqrt(5).*result1./par.theta_b(j)+5/3*result1.^2/par.theta_b(j)^2)+result2.*sqrt(5).*result1./(par.theta_b(j)^2);
                            end
                            
                            L=find(result3);
                            [idx_I,idx_J]=ind2sub(size(result3),L);
                            Id=[Id; idx_I+blocks(j+(i-1)*2)];
                            Jd=[Jd; idx_J+blocks(j+(i-1)*2)];
                            elements=[elements; result3(L)];
                        end
                        zero_density=(1-length(Id)/(N^2))*100;
                        if obj.tapering||zero_density>60
                            d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+j}=sparse(Id,Jd,elements,N,N);
                        else
                            d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+j}=zeros(N);
                            temp=sub2ind([N,N],Id,Jd);
                            d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+j}(temp)=elements;
                        end
                    end
                else
                    d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+1}=stem_misc.D_apply(result,aj_bp,'b').*d/(par.theta_b^2);
                end
            end
            %v_bp
            if n_bp_v>0
                if par.pixel_correlated
                    error('The variance-covariance matrix for cross-correlated pixel data is not available yet');
                    z=1;
                    for j=1:q
                        for i=j+1:q
                            Id=[];
                            Jd=[];
                            elements=[];
                            for h=1:2
                                idx=blocks(j+(h-1)*2)+1:blocks(j+(h-1)*2);
                                result1=result(idx,idx);
                                result2=stem_misc.D_apply(result1,aj_bp(idx),'b')/par.v_b(i,j); 
                                L=find(result2);
                                [idx_I,idx_J]=ind2sub(size(result2),L);
                                Id=[Id; idx_I+blocks(j+(h-1)*2)];
                                Jd=[Jd; idx_J+blocks(j+(h-1)*2)];
                                elements=[elements; result2(L)];
                            end
                            zero_density=(1-length(Id)/(N^2))*100;
                            if obj.tapering||zero_density>60
                                d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+z}=sparse(Id,Jd,elements,N,N);
                            else
                                d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+z}=zeros(N);
                                temp=sub2ind([N,N],Id,Jd);
                                d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+z}(temp)=elements;
                            end
                            z=z+1;
                        end
                    end
                else
                    %nothing as the matrix V_b is not estimated
                end
            end
            
            %alpha_p
            if n_p_alpha>0
                z=1;
                for k=1:par.k
                    for j=1:q
                        Id=[];
                        Jd=[];
                        elements=[];
                        for i=1:q
                            temp=sigma_W_p{k}(blocks(i)+1:blocks(i+1),blocks(j)+1:blocks(j+1));
                            if i==j
                                temp2=2*par.alpha_p(i,k)*temp;
                                L=find(temp2);
                                [idx_I,idx_J]=ind2sub(size(temp2),L);
                                Id=[Id; idx_I+blocks(i)];
                                Jd=[Jd; idx_J+blocks(j)];
                                elements=[elements; temp2(L)];
                            else
                                temp2=par.alpha_p(i,k)*temp;
                                L=find(temp2);
                                [idx_I,idx_J]=ind2sub(size(temp2),L);
                                Id=[Id; idx_I+blocks(i)];
                                Jd=[Jd; idx_J+blocks(j)];
                                elements=[elements; temp2(L)];
                                Id=[Id; idx_J+blocks(j)];
                                Jd=[Jd; idx_I+blocks(i)];
                                elements=[elements; temp2(L)];
                            end
                        end
                        zero_density=(1-length(Id)/(N^2))*100;
                        if obj.tapering||zero_density>60
                            d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+z}=sparse(Id,Jd,elements,Np,Np); %note that N,N is not used
                        else
                            d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+z}=zeros(Np);
                            temp=sub2ind([Np,Np],Id,Jd);
                            d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+z}(temp)=elements;
                        end
                        z=z+1;
                    end
                end
            end
            
            %theta_p
            if n_p_theta>0
                for k=1:par.k
                    if strcmp(par.correlation_type,'exponential')
                        d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+k}=stem_misc.D_apply(sigma_W_p{k}.*obj.stem_data.DistMat_p/(par.theta_p(k)^2),aj_p(:,k),'b');
                    elseif strcmp(par.correlation_type,'matern32')
                        d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+k}=stem_misc.D_apply(-sqrt(3)*obj.stem_data.DistMat_p/(par.theta_p(k)^2).*exp(-sqrt(3)*obj.stem_data.DistMat_p/par.theta_p(k))+sigma_W_p{k}.*sqrt(3).*obj.stem_data.DistMat_p/(par.theta_p(k)^2),aj_p(:,k),'b');
                    else
                        d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+k}=stem_misc.D_apply((-sqrt(5)*obj.stem_data.DistMat_p/par.theta_p(k)^2-10/3.*obj.stem_data.DistMat_p.^2/par.theta_p(k)^3).*exp(-sqrt(5)*obj.stem_data.DistMat_p/par.theta_p(k))+sigma_W_p{k}.*sqrt(5).*obj.stem_data.DistMat_p/(par.theta_p(k)^2),aj_p(:,k),'b');
                    end
                    if stem_misc.zero_density(d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+k})>60
                        d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+k}=sparse(d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+k});
                    end
                end
            end
            %v_p
            if n_p_v>0
                z=1;
                for k=1:par.k
                    for j=1:q
                        for i=j+1:q
                            Id=[];
                            Jd=[];
                            elements=[];
                            %since the block is extra-diagonal it has two separated D_apply!
                            result1=sigma_W_p{k}(blocks(j)+1:blocks(j+1),blocks(i)+1:blocks(i+1));
                            result2=stem_misc.D_apply(result1,aj_p(blocks(j)+1:blocks(j+1),k),'l');
                            result2=stem_misc.D_apply(result2,aj_p(blocks(i)+1:blocks(i+1),k),'r');
                            result2=result2/par.v_p(i,j,k);
                            L=find(result2);
                            [idx_I,idx_J]=ind2sub(size(result2),L);
                            Id=[Id; idx_I+blocks(j)];
                            Jd=[Jd; idx_J+blocks(i)];
                            elements=[elements; result2(L)];
                            Id=[Id; idx_J+blocks(i)];
                            Jd=[Jd; idx_I+blocks(j)];
                            elements=[elements; result2(L)];
                            zero_density=(1-length(Id)/(N^2))*100;
                            if obj.tapering||zero_density>60
                                d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+z}=sparse(Id,Jd,elements,Np,Np);
                            else
                                d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+z}=zeros(Np);
                                temp=sub2ind([Np,Np],Id,Jd);
                                d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+z}(temp)=elements;
                            end
                            z=z+1;
                        end
                    end
                end
            end
            
            %G and s2e
            for i=n_psi-n_time+1:n_psi
                d_Sgeo_prel{i}=sparse(N,N);
            end
                
            if not(data.X_tv) 
                %compute d_Sgeo in the time-invariant case
                for i=1:n_beta+n_eps
                    d_Sgeo{i}=d_Sgeo_prel{i};
                end
                if n_bp_alpha>0
                    for j=1:n_bp_alpha
                        d_Sgeo{n_beta+n_eps+j}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+j},obj.stem_data.X_bp(:,1,1),'b');
                    end
                end
                if n_bp_theta>0
                    for j=1:n_bp_theta
                        d_Sgeo{n_beta+n_eps+n_bp_alpha+j}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+j},obj.stem_data.X_bp(:,1,1),'b');   
                    end
                end
                if n_bp_v>0
                    for j=1:n_bp_v
                        d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+j}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+j},obj.stem_data.X_bp(:,1,1),'b');                           
                    end
                end
                if n_p_alpha>0
                    z=1;
                    for k=1:par.k
                        for j=1:q
                            d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+z}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+z},obj.stem_data.X_p(:,1,1,k),'b');
                            L=find(d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+z});
                            [Id,Jd]=ind2sub(size(d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+z}),L);
                            elements=d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+z}(L);
                            d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+z}=sparse(Id,Jd,elements,N,N);
                            z=z+1;
                        end
                    end
                end
                if n_p_theta>0
                    for k=1:par.k
                        X_p_orlated=[obj.stem_data.X_p(:,1,1,k);zeros(N-size(obj.stem_data.X_p(:,1,1,k),1),1)];
                        d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+k}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+k},X_p_orlated,'b');
                    end
                end
                if n_p_v>0
                    z=1;
                    for k=1:par.k
                        for j=1:q*(q-1)/2
                            d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+z}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+z},obj.stem_data.X_p(:,1,1,k),'b');
                            L=find(d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+z});
                            [Id,Jd]=ind2sub(size(d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+z}),L);
                            elements=d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+z}(L);
                            d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+z}=sparse(Id,Jd,elements,N,N);
                            z=z+1;
                        end
                    end
                end
                for i=n_psi-n_time+1:n_psi
                    d_Sgeo{i}=d_Sgeo_prel{i};   
                end
            end
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Information matrix evaluation  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            disp('Hessian evaluation...');
            c0 = 0;
            IM=zeros(n_psi);
            tot=n_psi*(n_psi+1)/2*T;
            counter=1;
            
            for t=1:T
                if data.X_bp_tv
                    tBP=t;
                else
                    tBP=1;
                end
                if data.X_z_tv
                    tT=t;
                else
                    tT=1;
                end
                if data.X_beta_tv
                    tbeta=t;
                else
                    tbeta=1;
                end      
                if data.X_p_tv
                    tP=t;
                else
                    tP=1;
                end                  
                
                if data.X_tv
                    %compute sigma_geo in the time-variant case
                    if not(isempty(data.X_bp))
                        sigma_geo=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),data.X_bp(:,1,tBP),'b'),aj_bp,'b');
                    end
                    if not(isempty(data.X_p))
                        if isempty(data.X_bp)
                            if (obj.tapering)&&(p==0)
                                sigma_geo=spalloc(size(sigma_W_p{1},1),size(sigma_W_p{1},1),nnz(sigma_W_p{1}));
                            else
                                sigma_geo=zeros(N);
                            end
                        end
                        for k=1:size(data.X_p,4)
                            sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},data.X_p(:,1,tP,k),'b'),aj_p(:,k),'b');
                        end
                    end
                    if isempty(data.X_p)&&isempty(data.X_bp)
                        sigma_geo=sigma_eps;
                    else
                        sigma_geo=sigma_geo+sigma_eps;
                    end
                    
                    %compute d_Sgeo in the time-variant case
                    for i=1:n_beta+n_eps
                        d_Sgeo{i}=d_Sgeo_prel{i};
                    end
                    if n_bp_alpha>0
                        for j=1:n_bp_alpha
                            d_Sgeo{n_beta+n_eps+j}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+j},obj.stem_data.X_bp(:,1,tBP),'b');
                        end
                    end
                    if n_bp_theta>0
                        for j=1:n_bp_theta
                            d_Sgeo{n_beta+n_eps+n_bp_alpha+j}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+j},obj.stem_data.X_bp(:,1,tBP),'b');
                        end
                    end
                    if n_bp_v>0
                        for j=1:n_bp_v
                            d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+j}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+j},obj.stem_data.X_bp(:,1,tBP),'b');
                        end
                    end
                    if n_p_alpha>0
                        z=1;
                        for k=1:par.k
                            for j=1:q
                                d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+z}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+z},obj.stem_data.X_p(:,1,tP,k),'b');
                                L=find(d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+z});
                                [Id,Jd]=ind2sub(size(d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+z}),L);
                                elements=d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+z}(L);
                                d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+z}=sparse(Id,Jd,elements,N,N);
                                z=z+1;
                            end
                        end
                    end
                    if n_p_theta>0
                        for k=1:par.k
                            X_p_orlated=[obj.stem_data.X_p(:,1,tP,k);zeros(N-size(obj.stem_data.X_p(:,1,tP,k),1),1)];
                            d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+k}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+k},X_p_orlated,'b');
                        end
                    end
                    if n_p_v>0
                        z=1;
                        for k=1:par.k
                            for j=1:q*(q-1)/2
                                d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+z}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+z},obj.stem_data.X_p(:,1,tP,k),'b');
                                L=find(d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+z});
                                [Id,Jd]=ind2sub(size(d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+z}),L);
                                elements=d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+z}(L);
                                d_Sgeo{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+z}=sparse(Id,Jd,elements,N,N);
                                z=z+1;
                            end
                        end
                    end
                    for i=n_psi-n_time+1:n_psi
                        d_Sgeo{i}=d_Sgeo_prel{i};
                    end
                end
                
                Lt=not(isnan(data.Y(:,t)));
                if n_time>0
                    X_z_orlated=data.X_z(:,:,tT);
                    X_z_orlated=[X_z_orlated;zeros(N-size(X_z_orlated,1),size(X_z_orlated,2))];
                    sigma_t_Lt=X_z_orlated(Lt,:)*st_kalmanfilter_result.Pk_f(:,:,t)*X_z_orlated(Lt,:)'+sigma_geo(Lt,Lt);
                else
                    sigma_t_Lt=sigma_geo(Lt,Lt);
                end
                
                et=zeros(N,1);
                if n_time>0
                    if n_beta>0
                        X_beta_orlated=data.X_beta(:,:,tbeta);
                        X_beta_orlated=[X_beta_orlated;zeros(N-size(X_beta_orlated,1),size(X_beta_orlated,2))];
                        et(Lt)=data.Y(Lt,t)-X_beta_orlated(Lt,:)*par.beta-X_z_orlated(Lt,:)*st_kalmanfilter_result.zk_f(:,t);
                    else
                        et(Lt)=data.Y(Lt,t)-X_z_orlated(Lt,:)*st_kalmanfilter_result.zk_f(:,t);
                    end
                else
                    if n_beta>0
                        X_beta_orlated=data.X_beta(:,:,tbeta);
                        X_beta_orlated=[X_beta_orlated;zeros(N-size(X_beta_orlated,1),size(X_beta_orlated,2))];
                        et(Lt)=data.Y(Lt,t)-X_beta_orlated(Lt,:)*par.beta;
                    else
                        et(Lt)=data.Y(Lt,t);
                    end
                end
                
                clear sigma_eps
                clear sigma_W_b
                clear sigma_W_p
                clear result
                clear result1
                clear result2
                clear result3
                clear temp
                clear temp2
                clear Id
                clear Jd
                clear L
                clear X_p_orlated
                clear aj_p
                clear aj_bp
                clear d
                clear elements
                clear idx
                clear idx_I
                clear idx_J
                
                if t==1
                    %d_P
                    for i=n_psi-n_time_s2e+1:n_psi
                        d_P(:,:,i)=d_s2e(:,:,i);
                    end
                    
                    %d_Z is zero for each parameter at t=1
                    
                    %d_St
                    for i=1:n_psi-n_time
                        %with respect to sigma_eps
                        d_St_Lt{i}=d_Sgeo{i}(Lt,Lt);
                    end
                    
                    %with respect to G are zero
                    for i=n_psi-n_time+1:n_psi-n_time_s2e
                        d_St_Lt{i}=sparse(sum(Lt),sum(Lt));
                    end
                    
                    for i=n_psi-n_time_s2e+1:n_psi
                        %with respect to sigma_eta
                        d_St_Lt{i}=X_z_orlated(Lt,:)*d_s2e(:,:,i)*X_z_orlated(Lt,:)';
                        if stem_misc.zero_density(d_St_Lt{i})>60
                            d_St_Lt{i}=sparse(d_St_Lt{i});
                        end
                    end
                    
                    %d_J
                    if n_time>0
                        for i=1:n_beta
                            d_J_Lt{i}=sparse(p,sum(Lt));
                        end
                        for i=n_beta+1:n_psi-n_time
                            d_J_Lt{i}=(-st_kalmanfilter_result.J(:,Lt,t)*d_Sgeo{i}(Lt,Lt))/sigma_t_Lt;
                        end
                        
                        for i=n_psi-n_time+1:+n_psi-n_time+n_time_G
                            %with respect to G
                            d_J_Lt{i}=(d_G(:,:,i)*st_kalmanfilter_result.Pk_f(:,:,t)*X_z_orlated(Lt,:)')/sigma_t_Lt;
                        end
                        
                        for i=n_psi-n_time_s2e+1:n_psi
                            %with respect to sigma_eta
                            d_J_Lt{i}=(par.G*d_P(:,:,i)*X_z_orlated(Lt,:)')/sigma_t_Lt;
                        end
                    end
                    
                    %d_e
                    for i=1:n_beta
                        d_e(:,i)=-X_beta_orlated*d_beta(:,1,i);
                    end
                else
                    Lt1=not(isnan(data.Y(:,t-1)));
                    %d_e
                    for i=1:n_psi
                        if (i<=n_beta)
                            if n_time>0
                                d_e(:,i)=-X_beta_orlated*d_beta(:,1,i)-X_z_orlated*d_Z_lag(:,i);
                            else
                                d_e(:,i)=-X_beta_orlated*d_beta(:,1,i);
                            end
                        else
                            if n_time>0
                                d_e(:,i)=-X_z_orlated*d_Z_lag(:,i);
                            end
                        end
                    end
                    
                    %d_Z
                    if n_time>0
                        for i=1:n_beta
                            d_Z(:,i)=par.G*d_Z_lag(:,i)+J(:,Lt1,t-1)*d_e_lag(Lt1,i);
                        end
                        for i=n_beta+1:n_psi-n_time
                            d_Z(:,i)=par.G*d_Z_lag(:,i)+d_J_lag_Lt1{i}*et_lag(Lt1)+J(:,Lt1,t-1)*d_e_lag(Lt1,i);
                        end
                        for i=n_psi-n_time+1:n_psi-n_time+n_time_G
                            d_Z(:,i)=d_G(:,:,i)*st_kalmanfilter_result.zk_f(:,t-1)+par.G*d_Z_lag(:,i)+d_J_lag_Lt1{i}*et_lag(Lt1)+J(:,Lt1,t-1)*d_e_lag(Lt1,i);
                        end
                        for i=n_psi-n_time_s2e+1:n_psi
                            d_Z(:,i)=par.G*d_Z_lag(:,i)+d_J_lag_Lt1{i}*et_lag(Lt1)+J(:,Lt1,t-1)*d_e_lag(Lt1,i);
                        end
                        
                        %d_P
                        for i=1:n_psi
                            d_P(:,:,i)=d_G(:,:,i)*(eye(p)-J(:,Lt1,t-1)*X_z_orlated(Lt1,:))*st_kalmanfilter_result.Pk_f(:,:,t-1)*par.G'-...
                                par.G*d_J_lag_Lt1{i}*X_z_orlated(Lt1,:)*st_kalmanfilter_result.Pk_f(:,:,t-1)*par.G'+...
                                par.G*(eye(p)-J(:,Lt1,t-1)*X_z_orlated(Lt1,:))*d_P_lag(:,:,i)*par.G'+...
                                par.G*(eye(p)-J(:,Lt1,t-1)*X_z_orlated(Lt1,:))*st_kalmanfilter_result.Pk_f(:,:,t-1)*d_G(:,:,i)'+d_s2e(:,:,i);
                        end
                    end
                    
                    %d_St
                    for i=1:n_psi
                        if n_time>0
                            d_St_Lt{i}=X_z_orlated(Lt,:)*d_P(:,:,i)*X_z_orlated(Lt,:)'+d_Sgeo{i}(Lt,Lt);
                        else
                            d_St_Lt{i}=d_Sgeo{i}(Lt,Lt);
                        end
                    end
                    
                    %d_J
                    if n_time>0
                        for i=1:n_psi
                            d_J_Lt{i}=(d_G(:,:,i)*st_kalmanfilter_result.Pk_f(:,:,t)*X_z_orlated(Lt,:)'+par.G*d_P(:,:,i)*X_z_orlated(Lt,:)'-J(:,Lt,t)*d_St_Lt{i})/sigma_t_Lt;
                        end
                    end
                end
                
                if n_time>0
                    temp1=X_z_orlated(Lt,:)*st_kalmanfilter_result.Pk_f(:,:,t)*X_z_orlated(Lt,:)';
                    sigma_t_Lt=temp1+sigma_geo(Lt,Lt);
                else
                    sigma_t_Lt=sigma_geo(Lt,Lt);
                end
                
                clear temp0
                clear temp1
                clear temp2
                d_e_Lt=d_e(Lt,:);
                if issparse(sigma_t_Lt)
                    r = symamd(sigma_t_Lt);
                    c=chol(sigma_t_Lt(r,r));
                    temp0(r,:)=stem_misc.chol_solve(c,d_e_Lt(r,:));
                else
                    c=chol(sigma_t_Lt);
                    temp0=stem_misc.chol_solve(c,d_e_Lt);
                end
                
                for i=1:n_psi
                    temp1{i}=d_e_Lt(:,i)'*temp0;
                    d_St_i_Lt=d_St_Lt{i};
                    if nnz(d_St_i_Lt)>0
                        if issparse(sigma_t_Lt)
                            temp2{i}(r,:)=full(stem_misc.chol_solve(c,d_St_i_Lt(r,:)));
                        else
                            temp2{i}=full(stem_misc.chol_solve(c,d_St_i_Lt));
                        end
                    else
                        temp2{i}=spalloc(size(sigma_t_Lt,1),size(sigma_t_Lt,2),0);
                    end
                end
                for i=1:n_psi
                    blocks=0:50:size(temp2{i},1);
                    if not(blocks(end)==size(temp2{i},1))
                        blocks=[blocks size(temp2{i},1)];
                    end
                    nnz_temp2_i=nnz(temp2{i});
                    for j=i:n_psi
                        IM(i,j)=IM(i,j)+temp1{i}(j);
                        if (nnz_temp2_i>0)&&(nnz(temp2{j})>0)
                            sumtrace=0;
                            for z=1:length(blocks)-1
                                idx=blocks(z)+1:blocks(z+1);
                                sumtrace=sumtrace+trace(temp2{i}(idx,:)*temp2{j}(:,idx));
                            end
                            IM(i,j)=IM(i,j)+0.5*sumtrace+0.25*trace(temp2{i})*trace(temp2{j});
                        end
                        counter=counter+1;
                        if (mod(counter,100)==0)||(counter==tot)
                            if mod(round(counter/tot*100),1)==0&&c0~=round(counter/tot*100)
                                c0 = round(counter/tot*100);
                                disp(['Hessian evaluation: ',num2str(round(counter/tot*100)),'% completed']);
                            end
                        end
                    end
                end
                d_P_lag=d_P;
                d_J_lag_Lt1=d_J_Lt;
                d_Z_lag=d_Z;
                d_e_lag=d_e;
                et_lag=et;
            end
            IM=IM+triu(IM,1)';
            obj.stem_EM_result.varcov=inv(IM);
        end      
        
        function set_initial_values(obj,stem_par)
            %DESCRIPTION: set the initial values of the model parameters
            %
            %INPUT
            %
            %obj - [stem_model object] (1x1)
            %
            %OUTPUT
            %
            %none: the stem_par_initial property is updated                
            if not(isa(stem_par,'stem_par'))
                error('The input argument must be of class stem_par');
            end
            obj.stem_par_initial=stem_par;
        end
        
        %Export functions. Useful to avoid access to the properties of nested objects
        function N = N(obj)
            N=obj.stem_data.N();
        end
        
        function Nb = Nb(obj)
            Nb=obj.stem_data.Nb();
        end
        
        function Np = Np(obj)
            Np=obj.stem_data.Np();
        end
        
        function T = T(obj)
            T=obj.stem_data.T();
        end
        
        function nvar=nvar(obj)
            nvar=obj.stem_data.nvar();
        end
        
        function dim=dim(obj)
            dim=obj.stem_data.dim();
        end
        
        %Initial values estimation functions (only beta at the moment)
        function [beta0] = get_beta0(obj)
            if obj.stem_par.n_beta>0
                N = size(obj.stem_data.X_beta,1);
                y = obj.stem_data.Y(1:N,:);
                y=y(:);
                T = obj.T;
                x = zeros(N*T,size(obj.stem_data.X_beta,2));
                
                for t=1:T
                    if size(obj.stem_data.X_beta,3)==T
                        tT=t;
                    else
                        tT=1;
                    end
                    x((t-1)*N+1:t*N,:)=obj.stem_data.X_beta(:,:,tT);
                end
                
                L=not(isnan(y));
                beta0 = (x(L,:)'*x(L,:))\x(L,:)'*y(L);
            else
                disp('WARNING: the model does not include data to estimate beta');
                beta0=[];
            end
        end
      
        %Class set functions 
        function set.stem_data(obj,stem_data)
            if strcmp(class(stem_data),'stem_data')
                obj.stem_data=stem_data;
            else
                error('The argument must be of class stem_data');
            end
        end
        
        function set.stem_par(obj,stem_par)
            if not(isa(stem_par,'stem_par'))
                error('stem_par_initial must be of class stem_par');
            end
            obj.stem_par=stem_par;
            if not(length(stem_par.beta)==size(obj.stem_data.X_beta,2))
                error(['The length of beta in stem_par must be equal to ',num2str(size(obj.stem_data.X_beta,2))]);
            end
        end        
       
        function set.stem_par_initial(obj,stem_par_initial)
            if not(isa(stem_par_initial,'stem_par'))
                error('stem_par_initial must be of class stem_par');
            else
                if obj.stem_par.time_diagonal
                    stem_par_initial.G=diag(diag(stem_par_initial.G));
                    stem_par_initial.sigma_eta=diag(diag(stem_par_initial.sigma_eta));
                end
                obj.stem_par_initial=stem_par_initial;
            end
        end
    end

end
